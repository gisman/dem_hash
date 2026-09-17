from __future__ import annotations

import argparse
import csv
import json
import math
import resource
import time
from contextlib import nullcontext
from pathlib import Path
from typing import Any, Iterator

import geohash2
import numpy as np
from numba import njit

from .dem import iter_dem_points
from .geohash_index import encode_geohash, iter_geohash_records
from .leveldb_store import LevelDBStore

try:
    import psycopg
except ImportError:  # pragma: no cover
    psycopg = None


POSTGIS_9CELL_SQL = """
WITH point_geom AS (
    SELECT ST_Transform(
        ST_SetSRID(ST_MakePoint(%s, %s), 4326),
        (SELECT ST_SRID(rast) FROM dem LIMIT 1)
    ) AS geom
),
center_pixel AS (
    SELECT
        d.rast,
        ST_WorldToRasterCoordX(d.rast, p.geom) AS cx,
        ST_WorldToRasterCoordY(d.rast, p.geom) AS cy
    FROM dem d
    JOIN point_geom p ON ST_Intersects(d.rast, p.geom)
),
offsets AS (
    SELECT dx, dy
    FROM generate_series(-1, 1) AS dx
    CROSS JOIN generate_series(-1, 1) AS dy
)
SELECT
    o.dx,
    o.dy,
    ST_Value(cp.rast, 1, cp.cx + o.dx, cp.cy + o.dy, true) AS elevation,
    ST_RasterToWorldCoordX(cp.rast, cp.cx + o.dx, cp.cy + o.dy) AS pixel_x,
    ST_RasterToWorldCoordY(cp.rast, cp.cx + o.dx, cp.cy + o.dy) AS pixel_y
FROM center_pixel cp
CROSS JOIN offsets o
ORDER BY o.dy, o.dx;
"""


POSTGIS_9CELL_BATCH_COUNT_SQL = """
WITH input_points AS (
    SELECT
        row_number() OVER () AS pid,
        lon,
        lat
    FROM unnest(%s::double precision[], %s::double precision[]) AS t(lon, lat)
),
point_geom AS (
    SELECT
        ip.pid,
        ST_Transform(
            ST_SetSRID(ST_MakePoint(ip.lon, ip.lat), 4326),
            (SELECT ST_SRID(rast) FROM dem LIMIT 1)
        ) AS geom
    FROM input_points ip
),
center_pixel AS (
    SELECT
        pg.pid,
        d.rast,
        ST_WorldToRasterCoordX(d.rast, pg.geom) AS cx,
        ST_WorldToRasterCoordY(d.rast, pg.geom) AS cy
    FROM dem d
    JOIN point_geom pg ON ST_Intersects(d.rast, pg.geom)
),
offsets AS (
    SELECT dx, dy
    FROM generate_series(-1, 1) AS dx
    CROSS JOIN generate_series(-1, 1) AS dy
)
SELECT COUNT(*)
FROM center_pixel cp
CROSS JOIN offsets o;
"""


@njit(cache=True)
def _solve_linear_system_gaussian(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    n = a.shape[0]
    aug = np.empty((n, n + 1), dtype=np.float64)

    for i in range(n):
        for j in range(n):
            aug[i, j] = a[i, j]
        aug[i, n] = b[i]

    for col in range(n):
        pivot_row = col
        pivot_abs = abs(aug[col, col])
        for row in range(col + 1, n):
            cand = abs(aug[row, col])
            if cand > pivot_abs:
                pivot_abs = cand
                pivot_row = row

        if pivot_abs < 1e-12:
            aug[col, col] += 1e-8
            pivot_abs = abs(aug[col, col])

        if pivot_row != col:
            for k in range(col, n + 1):
                tmp = aug[col, k]
                aug[col, k] = aug[pivot_row, k]
                aug[pivot_row, k] = tmp

        pivot = aug[col, col]
        if abs(pivot) < 1e-14:
            pivot = 1e-14
            aug[col, col] = pivot

        for row in range(col + 1, n):
            factor = aug[row, col] / pivot
            if factor == 0.0:
                continue
            for k in range(col, n + 1):
                aug[row, k] -= factor * aug[col, k]

    x = np.zeros(n, dtype=np.float64)
    for i in range(n - 1, -1, -1):
        rhs = aug[i, n]
        for j in range(i + 1, n):
            rhs -= aug[i, j] * x[j]
        denom = aug[i, i]
        if abs(denom) < 1e-14:
            denom = 1e-14
        x[i] = rhs / denom

    return x


@njit(cache=True)
def _uk_predict_local_xy(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> float:
    """주변 샘플 점들로부터 현 위치의 고도값을 추정하는 로컬 UK(Universal Kriging) 예측기.

    이 함수는 입력 좌표계에서 각 샘플 지점이 (x, y, z) = (동서거리, 남북거리, 고도)
    를 가진다고 보고, 기준 위치에서의 고도를 예측한다. 핵심 아이디어는 다음과 같다.

    1) 샘플 값 z는 단순한 평면식이 아니라 공간적 상관관계(spatial correlation)를 가진
       랜덤 변수로 가정한다.
    2) z를 "전역 평균 + 선형 drift(기울기) + 공간적 오차" 형태로 분해한다.
       여기서는 x, y의 선형 항을 drift로 두고, 그 나머지를 지점 간 거리 기반의
       Gaussian covariance로 모델링한다.
    3) 각 샘플 점의 값이 알려져 있을 때, 예측 지점과의 거리별 covariance를 계산하고,
       Kriging 정상방정식을 풀어 예측 계수(lambda)를 구한다.
    4) 예측값은 lambda와 샘플 값의 가중합으로 계산한다.

    코드의 구현은 "universal kriging with linear drift"의 간단한 버전으로,
    (1) 평균/기울기 조건식, (2) covariance matrix, (3) 선형 시스템 풀이를 직접 구성한다.
    """
    n = x.shape[0]
    if n == 0:
        # 입력 포인트가 하나도 없으면 예측 불가
        return np.nan

    # 모든 샘플 점이 (0, 0)에 몰려 있으면, 그 위치의 실제 고도값을 그대로 반환한다.
    # 이 경우에는 거리 기반 covariances가 의미가 없기 때문에 특별 처리한다.
    zero_mask = (x == 0.0) & (y == 0.0)
    if np.any(zero_mask):
        return float(z[zero_mask][0])

    # 1) z값의 평균을 계산한다.
    #    이 값은 추정 과정에서 "전역 수준(global mean)" 또는 baseline의 근사치로 사용된다.
    z_mean = float(np.mean(z))

    # 2) z의 분산을 추정한다.
    #    이것은 Gaussian covariance의 sill(최대 분산) 값에 해당한다.
    #    즉, 지점 간 거리별 상관관계가 완전히 끊어지는 수준의 분산을 가정한다.
    z_var = float(np.mean((z - z_mean) ** 2))
    if z_var < 1e-8:
        # 너무 작은 분산은 수치적으로 불안정하므로 최소값을 보장한다.
        z_var = 1e-8

    # 3) 가장 먼 샘플 점 사이의 거리를 계산한다.
    #    이것은 지점 간 상관관계가 얼마나 멀리까지 유지되는지를 나타내는 range 파라미터의
    #    근사값으로 사용된다. 보통 range는 cov(d)=exp(-d/range) 가 0에 가까워지는 거리 척도이다.
    dx_mat = x[:, None] - x[None, :]
    dy_mat = y[:, None] - y[None, :]
    dist_mat = np.sqrt(dx_mat * dx_mat + dy_mat * dy_mat)
    max_dist = float(np.max(dist_mat))

    if max_dist < 1.0:
        # 너무 작은 유효 거리는 의미가 없으므로 최소 길이를 지정한다.
        max_dist = 1.0

    # 4) 가우시안 커널 기반 covariance를 설정한다.
    #    cov(d) = sill * exp(-d / range)
    #    여기서 sill는 전역 분산, range는 거리에서 상관관계가 얼마나 빠르게 떨어지는지를 조절한다.
    sill = z_var
    range_param = max_dist

    # nugget(미세 노이즈) 항: 동일 지점에서의 자기 상관이 전부 정확히 같은 값이 아니라
    # 아주 작은 추가 분산을 가진다고 가정하기 위해 넣는다. 수치 안정성을 위해 tiny 값도 추가한다.
    reg = sill * 1e-6 + 1e-8

    # 5) Universal Kriging의 정상방정식 크기 설정
    #    n개의 관측값에 대해 기본 선형 drift 항을 3개(x, y, 상수)로 두므로,
    #    총 크기는 n + 3이다. 실제로 마지막 3개 행/열은 다음의 제약조건을 위해 사용된다.
    #       - sum(lambda) = 1       (상수 drift 조건)
    #       - sum(lambda * x_i) = x0 (x drift 조건)
    #       - sum(lambda * y_i) = y0 (y drift 조건)
    #    이는 예측 위치가 포함된 선형 drift를 보존하도록 만들기 위한 조건이다.
    p = 3
    size = n + p
    mat = np.zeros((size, size), dtype=np.float64)
    rhs = np.zeros(size, dtype=np.float64)

    # 6) 관측치 간 covariance 행렬을 채운다.
    #    K_ij = cov(distance between sample i and sample j)
    #    대각선 요소는 자기 상관을 포함하고, nugget(reg)를 더해 수치적으로 안정화한다.
    cov_mat = sill * np.exp(-dist_mat / range_param)
    for i in range(n):
        cov_mat[i, i] += reg
    mat[:n, :n] = cov_mat

    # 각 관측점에 대해 선형 drift 항의 기여를 추가한다.
    #    [1, x_i, y_i] 를 row의 마지막 3개 열에 적재한다.
    mat[:n, n] = 1.0
    mat[:n, n + 1] = x
    mat[:n, n + 2] = y

    # 예측 지점과 관측값 사이의 covariances를 우변 rhs에 채운다.
    # 이 함수는 입력에 "예측 위치가 아닌, 관측 위치들의 값만 알려진 상태"에서
    # 예측 위치의 z 값을 구해야 하므로, 각 관측점 i에 대해 아래의 상관성을 rhs[i]에 담는다.
    #   cov(distance(sample_i, prediction_point))
    d0 = np.sqrt(x * x + y * y)
    rhs[:n] = sill * np.exp(-d0 / range_param)

    # 7) 마지막 3개 행/열은 drift 조건식을 표현한다.
    #    관측점에 대한 [1, x_i, y_i] 계수와 예측점에서의 [1, x0, y0] 조건을 맞춘다.
    #    우변 rhs[n], rhs[n+1], rhs[n+2]는 예측 지점의 선형 drift 값으로 1,0,0 을 넣는다.
    #    이 조건에 의해 kriging 가중치들이 합이 1이고, x/y 기울기까지 보존되도록 강제한다.
    mat[n, :n] = 1.0
    mat[n + 1, :n] = x
    mat[n + 2, :n] = y

    rhs[n] = 1.0
    rhs[n + 1] = 0.0
    rhs[n + 2] = 0.0

    # 8) 정상방정식을 풀어 계수 sol을 구한다.
    #    sol의 앞 n개 원소는 각 샘플의 kriging weight lambda_i 이고,
    #    마지막 3개는 Lagrange multiplier이다.
    sol = _solve_linear_system_gaussian(mat, rhs)

    # 9) 예측값은 가중합으로 계산한다.
    #    pred = Σ λ_i * z_i
    #    여기서 λ_i는 샘플 값과 지점 간 covariances를 기반으로 산출된 최적 weight이다.
    pred = float(np.dot(sol[:n], z))
    return pred


def _estimate_elevation_from_store_kriging(
    store: LevelDBStore,
    lat: float,
    lon: float,
    precision: int,
    per_hash_limit: int,
) -> dict:
    center_hash, search_hashes = _search_hashes_for_point(lat, lon, precision=precision)

    candidates: list[dict] = []
    for geohash_prefix in search_hashes:
        rows = store.search_by_prefix(geohash_prefix, limit=per_hash_limit)
        candidates.extend(rows)

    if not candidates:
        return {
            "query": {"latitude": lat, "longitude": lon},
            "precision": precision,
            "center_geohash": center_hash,
            "searched_geohashes": search_hashes,
            "candidate_count": 0,
            "estimated_elevation": None,
            "message": "no candidate points found in level-7 3x3 geohash cells",
        }

    lats = np.array([float(item["latitude"]) for item in candidates], dtype=np.float64)
    lons = np.array([float(item["longitude"]) for item in candidates], dtype=np.float64)
    elev = np.array([float(item["elevation"]) for item in candidates], dtype=np.float64)

    radius = 6_371_000.0
    deg_to_rad = math.pi / 180.0
    lat0_rad = lat * deg_to_rad
    x = (lons - lon) * deg_to_rad * radius * math.cos(lat0_rad)
    y = (lats - lat) * deg_to_rad * radius

    estimate = _uk_predict_local_xy(x, y, elev)
    estimated_elevation = None if np.isnan(estimate) else float(estimate)

    return {
        "query": {"latitude": lat, "longitude": lon},
        "precision": precision,
        "center_geohash": center_hash,
        "searched_geohashes": search_hashes,
        "candidate_count": len(candidates),
        "points_used_for_kriging": len(candidates),
        "kriging_model": "universal-linear-drift",
        "estimated_elevation": estimated_elevation,
    }


def _records_from_dem(dem_path: str, precision: int, stride: int, skip_rows: int = 0):
    points = iter_dem_points(dem_path, stride=stride, skip_rows=skip_rows)
    return iter_geohash_records(points, precision=precision)


def cmd_build_leveldb(args: argparse.Namespace) -> None:
    records = _records_from_dem(
        args.dem_path, args.precision, args.stride, skip_rows=args.skip_rows
    )
    with LevelDBStore(args.db_path) as store:
        count = store.put_many(records)
    print(f"saved {count} records to LevelDB: {args.db_path}")


def _haversine_m(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    radius = 6_371_000.0
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)

    a = (
        math.sin(dphi / 2.0) ** 2
        + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2.0) ** 2
    )
    c = 2.0 * math.atan2(math.sqrt(a), math.sqrt(1.0 - a))
    return radius * c


def _search_hashes_for_point(
    lat: float, lon: float, precision: int
) -> tuple[str, list[str]]:
    center_hash = encode_geohash(lat, lon, precision=precision)
    center_lat, center_lon, lat_err, lon_err = geohash2.decode_exactly(center_hash)
    step_lat = lat_err * 2.0
    step_lon = lon_err * 2.0

    search_hashes: list[str] = []
    for row_offset in (1, 0, -1):
        for col_offset in (-1, 0, 1):
            neighbor_hash = encode_geohash(
                center_lat + (step_lat * row_offset),
                center_lon + (step_lon * col_offset),
                precision=precision,
            )
            if neighbor_hash not in search_hashes:
                search_hashes.append(neighbor_hash)

    return center_hash, search_hashes


def _estimate_elevation_from_store(
    store: LevelDBStore,
    lat: float,
    lon: float,
    precision: int,
    per_hash_limit: int,
) -> dict:
    center_hash, search_hashes = _search_hashes_for_point(lat, lon, precision=precision)

    candidates: list[dict] = []
    for geohash_prefix in search_hashes:
        rows = store.search_by_prefix(geohash_prefix, limit=per_hash_limit)
        candidates.extend(rows)

    if not candidates:
        return {
            "query": {"latitude": lat, "longitude": lon},
            "precision": precision,
            "center_geohash": center_hash,
            "searched_geohashes": search_hashes,
            "candidate_count": 0,
            "estimated_elevation": None,
            "message": "no candidate points found in level-7 3x3 geohash cells",
        }

    for point in candidates:
        point["distance_m"] = _haversine_m(
            lat, lon, point["latitude"], point["longitude"]
        )

    points_for_idw = sorted(candidates, key=lambda item: item["distance_m"])[:9]
    power = 2.0

    if points_for_idw[0]["distance_m"] == 0.0:
        estimated_elevation = points_for_idw[0]["elevation"]
        for point in points_for_idw:
            point["idw_weight"] = 1.0 if point["distance_m"] == 0.0 else 0.0
    else:
        weighted_sum = 0.0
        weight_sum = 0.0
        for point in points_for_idw:
            weight = 1.0 / (point["distance_m"] ** power)
            point["idw_weight"] = weight
            weighted_sum += weight * point["elevation"]
            weight_sum += weight
        estimated_elevation = weighted_sum / weight_sum if weight_sum > 0.0 else None

    return {
        "query": {"latitude": lat, "longitude": lon},
        "precision": precision,
        "center_geohash": center_hash,
        "searched_geohashes": search_hashes,
        "candidate_count": len(candidates),
        "points_used_for_idw": len(points_for_idw),
        "idw_power": power,
        "estimated_elevation": estimated_elevation,
        "neighbors": points_for_idw,
    }


def _iter_csv_coordinates(csv_path: str | Path) -> Iterator[tuple[float, float]]:
    path = Path(csv_path)

    with path.open("r", encoding="utf-8", newline="") as fp:
        # Header mode: supports latitude/longitude and lat/lon naming.
        dict_reader = csv.DictReader(fp)
        if dict_reader.fieldnames:
            normalized = {name.strip().lower(): name for name in dict_reader.fieldnames}
            pnu = normalized.get("pnu")
            lat_name = normalized.get("latitude") or normalized.get("lat")
            lon_name = normalized.get("longitude") or normalized.get("lon")

            if lat_name and lon_name:
                for row_no, row in enumerate(dict_reader, start=2):
                    lat_raw = row.get(lat_name)
                    lon_raw = row.get(lon_name)
                    if lat_raw is None or lon_raw is None:
                        raise ValueError(f"Missing lat/lon value at line {row_no}")
                    yield pnu, float(lat_raw), float(lon_raw)
                return

    # with path.open("r", encoding="utf-8", newline="") as fp:
    #     reader = csv.reader(fp)
    #     # The first row is treated as a title/header and always ignored.
    #     next(reader, None)
    #     for row_no, row in enumerate(reader, start=2):
    #         if not row or all(not cell.strip() for cell in row):
    #             continue
    #         if len(row) < 2:
    #             raise ValueError(f"Expected at least 2 columns at line {row_no}")
    #         yield float(row[0]), float(row[1])


def cmd_search_geohash(args: argparse.Namespace) -> None:
    if args.per_hash_limit < 1:
        raise ValueError("--per-hash-limit must be >= 1")

    precision = 7
    with LevelDBStore(args.db_path) as store:
        result = _estimate_elevation_from_store(
            store,
            lat=args.lat,
            lon=args.lon,
            precision=precision,
            per_hash_limit=args.per_hash_limit,
        )

    print(json.dumps(result, ensure_ascii=False, indent=2))


def cmd_search_csv(args: argparse.Namespace) -> None:
    if args.limit < -1:
        raise ValueError("--limit must be >= 0 or -1 (unlimited)")

    if args.per_hash_limit < 1:
        raise ValueError("--per-hash-limit must be >= 1")

    precision = 7
    processed_count = 0
    peak_memory_mb = 0.0
    output_path = Path(f"{args.csv_path}.out") if args.write_output else None

    with LevelDBStore(args.db_path) as store:
        # 제한이 없는 실행에서만 워밍업을 수행한다.
        if args.limit == -1:
            # 정확한 성능 측정을 위해 100건 정도 읽어본다.
            warmup_count = 0
            for lon, lat in _iter_csv_coordinates(args.csv_path):
                _estimate_elevation_from_store(
                    store,
                    lat=lat,
                    lon=lon,
                    precision=precision,
                    per_hash_limit=args.per_hash_limit,
                )
                warmup_count += 1
                if warmup_count >= 100:
                    break

        started_at = time.perf_counter()
        if output_path is None:
            # Fast path: no output file, only parse coordinates and estimate.
            for lon, lat in _iter_csv_coordinates(args.csv_path):
                if args.limit >= 0 and processed_count >= args.limit:
                    break

                result = _estimate_elevation_from_store(
                    store,
                    lat=lat,
                    lon=lon,
                    precision=precision,
                    per_hash_limit=args.per_hash_limit,
                )
                processed_count += 1
                if processed_count <= 200:
                    print(
                        f'lat={lat}, lon={lon}, elevation={result["estimated_elevation"]}'
                    )
        else:
            with Path(args.csv_path).open("r", encoding="utf-8", newline="") as in_fp:
                out_context = (
                    output_path.open("w", encoding="utf-8", newline="")
                    if output_path is not None
                    else nullcontext(None)
                )
                with out_context as out_fp:
                    reader = csv.reader(in_fp)
                    writer = csv.writer(out_fp) if out_fp is not None else None

                    first_row = next(reader, None)
                    if first_row is not None:
                        normalized = {
                            name.strip().lower(): idx
                            for idx, name in enumerate(first_row)
                        }
                        lat_idx = (
                            normalized["latitude"]
                            if "latitude" in normalized
                            else normalized.get("lat")
                        )
                        lon_idx = (
                            normalized["longitude"]
                            if "longitude" in normalized
                            else normalized.get("lon", normalized.get("lng"))
                        )

                        if writer is not None:
                            writer.writerow(first_row + ["elevation"])

                        for row_no, row in enumerate(reader, start=2):
                            if args.limit >= 0 and processed_count >= args.limit:
                                break

                            if lat_idx is not None and lon_idx is not None:
                                if lat_idx >= len(row) or lon_idx >= len(row):
                                    raise ValueError(
                                        f"Missing lat/lon value at line {row_no}"
                                    )
                                lon = float(row[lat_idx])
                                lat = float(row[lon_idx])
                            else:
                                if not row or all(not cell.strip() for cell in row):
                                    continue
                                if len(row) < 2:
                                    raise ValueError(
                                        f"Expected at least 2 columns at line {row_no}"
                                    )
                                lon = float(row[1])
                                lat = float(row[2])

                            result = _estimate_elevation_from_store(
                                store,
                                lat=lat,
                                lon=lon,
                                precision=precision,
                                per_hash_limit=args.per_hash_limit,
                            )
                            if writer is not None:
                                writer.writerow(row + [result["estimated_elevation"]])
                            processed_count += 1
                            if processed_count <= 200:
                                print(
                                    f'lat={lat}, lon={lon}, elevation={result["estimated_elevation"]}'
                                )

    elapsed_seconds = time.perf_counter() - started_at
    # ru_maxrss is KB on Linux; convert to MiB for readable output.
    peak_memory_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    print(
        f"elapsed_seconds={elapsed_seconds:.3f} count={processed_count} "
        f"peak_memory_mb={peak_memory_mb:.2f}"
    )
    if output_path is not None:
        print(f"output_path={output_path}")
    else:
        print("output_path=disabled")


def cmd_search_csv_kriging(args: argparse.Namespace) -> None:
    if args.limit < -1:
        raise ValueError("--limit must be >= 0 or -1 (unlimited)")

    if args.per_hash_limit < 1:
        raise ValueError("--per-hash-limit must be >= 1")

    precision = 7
    processed_count = 0
    output_path = Path(f"{args.csv_path}.out") if args.write_output else None

    with LevelDBStore(args.db_path) as store:
        if args.limit == -1:
            warmup_count = 0
            for lon, lat in _iter_csv_coordinates(args.csv_path):
                _estimate_elevation_from_store_kriging(
                    store,
                    lat=lat,
                    lon=lon,
                    precision=precision,
                    per_hash_limit=args.per_hash_limit,
                )
                warmup_count += 1
                if warmup_count >= 100:
                    break

        started_at = time.perf_counter()
        if output_path is None:
            for lon, lat in _iter_csv_coordinates(args.csv_path):
                if args.limit >= 0 and processed_count >= args.limit:
                    break

                result = _estimate_elevation_from_store_kriging(
                    store,
                    lat=lat,
                    lon=lon,
                    precision=precision,
                    per_hash_limit=args.per_hash_limit,
                )
                processed_count += 1
                if processed_count <= 200:
                    print(
                        f'lat={lat}, lon={lon}, elevation={result["estimated_elevation"]}'
                    )
        else:
            with Path(args.csv_path).open("r", encoding="utf-8", newline="") as in_fp:
                out_context = (
                    output_path.open("w", encoding="utf-8", newline="")
                    if output_path is not None
                    else nullcontext(None)
                )
                with out_context as out_fp:
                    reader = csv.reader(in_fp)
                    writer = csv.writer(out_fp) if out_fp is not None else None

                    first_row = next(reader, None)
                    if first_row is not None:
                        normalized = {
                            name.strip().lower(): idx
                            for idx, name in enumerate(first_row)
                        }
                        lat_idx = (
                            normalized["latitude"]
                            if "latitude" in normalized
                            else normalized.get("lat")
                        )
                        lon_idx = (
                            normalized["longitude"]
                            if "longitude" in normalized
                            else normalized.get("lon", normalized.get("lng"))
                        )

                        if writer is not None:
                            writer.writerow(first_row + ["elevation"])

                        for row_no, row in enumerate(reader, start=2):
                            if args.limit >= 0 and processed_count >= args.limit:
                                break

                            if lat_idx is not None and lon_idx is not None:
                                if lat_idx >= len(row) or lon_idx >= len(row):
                                    raise ValueError(
                                        f"Missing lat/lon value at line {row_no}"
                                    )
                                lon = float(row[lat_idx])
                                lat = float(row[lon_idx])
                            else:
                                if not row or all(not cell.strip() for cell in row):
                                    continue
                                if len(row) < 2:
                                    raise ValueError(
                                        f"Expected at least 2 columns at line {row_no}"
                                    )
                                lon = float(row[1])
                                lat = float(row[2])

                            result = _estimate_elevation_from_store_kriging(
                                store,
                                lat=lat,
                                lon=lon,
                                precision=precision,
                                per_hash_limit=args.per_hash_limit,
                            )
                            if writer is not None:
                                writer.writerow(row + [result["estimated_elevation"]])
                            processed_count += 1
                            if processed_count <= 200:
                                print(
                                    f'lat={lat}, lon={lon}, elevation={result["estimated_elevation"]}'
                                )

                            if processed_count % 100_000 == 0:
                                print(f"Processed {processed_count:,} rows")

    elapsed_seconds = time.perf_counter() - started_at
    peak_memory_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    print(
        f"elapsed_seconds={elapsed_seconds:.3f} count={processed_count} "
        f"peak_memory_mb={peak_memory_mb:.2f}"
    )
    if output_path is not None:
        print(f"output_path={output_path}")
    else:
        print("output_path=disabled")


def _load_coordinates_for_benchmark(csv_path: str | Path) -> list[tuple[float, float]]:
    path = Path(csv_path)
    coords: list[tuple[float, float]] = []

    with path.open("r", encoding="utf-8", newline="") as fp:
        reader = csv.DictReader(fp)
        if reader.fieldnames:
            normalized = {name.strip().lower(): name for name in reader.fieldnames}
            lat_name = (
                normalized.get("latitude")
                or normalized.get("lat")
                or normalized.get("y_coord")
                or normalized.get("y")
            )
            lon_name = (
                normalized.get("longitude")
                or normalized.get("lon")
                or normalized.get("lng")
                or normalized.get("x_coord")
                or normalized.get("x")
            )
            if lat_name and lon_name:
                for row in reader:
                    lat_raw = row.get(lat_name)
                    lon_raw = row.get(lon_name)
                    if lat_raw is None or lon_raw is None:
                        continue
                    coords.append((float(lat_raw), float(lon_raw)))

    if coords:
        return coords

    with path.open("r", encoding="utf-8", newline="") as fp:
        plain_reader = csv.reader(fp)
        next(plain_reader, None)
        for row in plain_reader:
            if len(row) < 3:
                continue
            coords.append((float(row[2]), float(row[1])))

    return coords


def _iter_sample_coordinates(
    coords: list[tuple[float, float]], sample_count: int
) -> Iterator[tuple[float, float]]:
    coord_count = len(coords)
    for idx in range(sample_count):
        yield coords[idx % coord_count]


def _run_hash_benchmark(
    store: LevelDBStore,
    coords: list[tuple[float, float]],
    sample_count: int,
    precision: int,
    per_hash_limit: int,
) -> tuple[float, int]:
    started_at = time.perf_counter()
    total_rows = 0

    for lat, lon in _iter_sample_coordinates(coords, sample_count):
        _, search_hashes = _search_hashes_for_point(lat, lon, precision=precision)
        for geohash_prefix in search_hashes:
            total_rows += len(
                store.search_by_prefix(geohash_prefix, limit=per_hash_limit)
            )

    elapsed_seconds = time.perf_counter() - started_at
    return elapsed_seconds, total_rows


def _run_postgis_benchmark(
    coords: list[tuple[float, float]],
    sample_count: int,
    conn: Any,
    batch_size: int,
) -> tuple[float, int]:
    if batch_size < 1:
        raise ValueError("batch_size must be >= 1")

    started_at = time.perf_counter()
    total_rows = 0

    lons: list[float] = []
    lats: list[float] = []

    with conn.cursor() as cur:
        for lat, lon in _iter_sample_coordinates(coords, sample_count):
            lats.append(lat)
            lons.append(lon)

            if len(lats) >= batch_size:
                cur.execute(POSTGIS_9CELL_BATCH_COUNT_SQL, (lons, lats))
                row = cur.fetchone()
                total_rows += int(row[0]) if row is not None else 0
                lats.clear()
                lons.clear()

        if lats:
            cur.execute(POSTGIS_9CELL_BATCH_COUNT_SQL, (lons, lats))
            row = cur.fetchone()
            total_rows += int(row[0]) if row is not None else 0

    elapsed_seconds = time.perf_counter() - started_at
    return elapsed_seconds, total_rows


def _parse_sample_sizes(sample_sizes: str) -> list[int]:
    values = [part.strip() for part in sample_sizes.split(",") if part.strip()]
    parsed = [int(value) for value in values]
    if not parsed or any(value < 1 for value in parsed):
        raise ValueError("--sample-sizes must include one or more positive integers")
    return parsed


def cmd_compare_hash_postgis(args: argparse.Namespace) -> None:
    if args.per_hash_limit < 1:
        raise ValueError("--per-hash-limit must be >= 1")
    if psycopg is None:
        raise RuntimeError(
            "psycopg is required for PostGIS benchmark. "
            "Install it with: uv add 'psycopg[binary]'"
        )

    sample_sizes = _parse_sample_sizes(args.sample_sizes)
    coords = _load_coordinates_for_benchmark(args.csv_path)
    if not coords:
        raise ValueError("No coordinates found in --csv-path")

    print(f"coordinate_pool={len(coords)} from={args.csv_path}")
    if max(sample_sizes) > len(coords):
        print(
            "sample_mode=cyclic-repeat (reuses coordinate pool for larger sample sizes)"
        )

    with LevelDBStore(args.db_path) as store:
        with psycopg.connect(
            host=args.pg_host,
            port=args.pg_port,
            user=args.pg_user,
            password=args.pg_password,
            dbname=args.pg_database,
            autocommit=True,
        ) as conn:
            if args.warmup > 0:
                warmup_count = min(args.warmup, len(coords))
                _run_hash_benchmark(
                    store,
                    coords,
                    sample_count=warmup_count,
                    precision=args.precision,
                    per_hash_limit=args.per_hash_limit,
                )
                _run_postgis_benchmark(
                    coords,
                    sample_count=warmup_count,
                    conn=conn,
                    batch_size=args.postgis_batch_size,
                )

            for sample_count in sample_sizes:
                hash_elapsed, _ = _run_hash_benchmark(
                    store,
                    coords,
                    sample_count=sample_count,
                    precision=args.precision,
                    per_hash_limit=args.per_hash_limit,
                )
                hash_qps = (
                    sample_count / hash_elapsed if hash_elapsed > 0 else float("inf")
                )

                postgis_elapsed, _ = _run_postgis_benchmark(
                    coords,
                    sample_count=sample_count,
                    conn=conn,
                    batch_size=args.postgis_batch_size,
                )
                postgis_qps = (
                    sample_count / postgis_elapsed
                    if postgis_elapsed > 0
                    else float("inf")
                )

                speed_ratio = (
                    hash_elapsed / postgis_elapsed
                    if postgis_elapsed > 0
                    else float("inf")
                )
                print(f"sample_count={sample_count}")
                print(
                    f"hash_elapsed_seconds={hash_elapsed:.6f} hash_qps={hash_qps:.2f}"
                )
                print(
                    f"postgis_elapsed_seconds={postgis_elapsed:.6f} "
                    f"postgis_qps={postgis_qps:.2f}"
                )
                print(f"hash_vs_postgis_time_ratio={speed_ratio:.4f}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="DEM geohash index utility")
    sub = parser.add_subparsers(dest="command", required=True)

    p_build = sub.add_parser("build-leveldb", help="Index DEM and store into LevelDB")
    p_build.add_argument("--dem-path", default="data/한반도90m_GRS80.img")
    p_build.add_argument("--db-path", default="data/dem.leveldb")
    p_build.add_argument("--precision", type=int, default=7)
    p_build.add_argument("--stride", type=int, default=10)
    p_build.add_argument("--skip-rows", type=int, default=0)
    p_build.set_defaults(func=cmd_build_leveldb)

    p_search_hash = sub.add_parser(
        "search-geohash", help="Estimate elevation from lat/lon using geohash+IDW"
    )
    p_search_hash.add_argument("--db-path", default="data/dem.leveldb")
    p_search_hash.add_argument("--lat", type=float, required=True)
    p_search_hash.add_argument("--lon", type=float, required=True)
    p_search_hash.add_argument("--per-hash-limit", type=int, default=5000)
    p_search_hash.set_defaults(func=cmd_search_geohash)

    p_search_csv = sub.add_parser(
        "search-csv", help="Process all coordinates in CSV and estimate elevations"
    )
    p_search_csv.add_argument("--db-path", default="data/dem.leveldb")
    p_search_csv.add_argument("--csv-path", default="test.csv")
    p_search_csv.add_argument("--limit", type=int, default=-1)
    p_search_csv.add_argument(
        "--write-output",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    p_search_csv.add_argument("--per-hash-limit", type=int, default=5)
    p_search_csv.set_defaults(func=cmd_search_csv)

    p_search_csv_kriging = sub.add_parser(
        "search-csv-kriging",
        help="Process all coordinates in CSV and estimate elevations using Universal Kriging",
    )
    p_search_csv_kriging.add_argument("--db-path", default="data/dem.leveldb")
    p_search_csv_kriging.add_argument("--csv-path", default="test.csv")
    p_search_csv_kriging.add_argument("--limit", type=int, default=-1)
    p_search_csv_kriging.add_argument(
        "--write-output",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    p_search_csv_kriging.add_argument("--per-hash-limit", type=int, default=5)
    p_search_csv_kriging.set_defaults(func=cmd_search_csv_kriging)

    p_compare = sub.add_parser(
        "compare-hash-postgis",
        help="Compare hash vs PostGIS query speed using 9-neighbor search",
    )
    p_compare.add_argument("--db-path", default="data/dem.leveldb")
    p_compare.add_argument("--csv-path", default="data/pnu_centroid_100k.csv")
    p_compare.add_argument("--precision", type=int, default=7)
    p_compare.add_argument("--per-hash-limit", type=int, default=5)
    p_compare.add_argument("--sample-sizes", default="100000,1000000")
    p_compare.add_argument("--warmup", type=int, default=1000)
    p_compare.add_argument("--pg-host", default="localhost")
    p_compare.add_argument("--pg-port", type=int, default=5432)
    p_compare.add_argument("--pg-user", default="postgres")
    p_compare.add_argument("--pg-password", default="xhddlfd0")
    p_compare.add_argument("--pg-database", default="cube_db")
    p_compare.add_argument("--postgis-batch-size", type=int, default=5000)
    p_compare.set_defaults(func=cmd_compare_hash_postgis)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
