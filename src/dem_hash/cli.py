from __future__ import annotations

import argparse
import csv
import json
import math
import time
from pathlib import Path
from typing import Iterator

import geohash2

from .dem import iter_dem_points
from .geohash_index import encode_geohash, iter_geohash_records
from .leveldb_store import LevelDBStore


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
            lat_name = normalized.get("latitude") or normalized.get("lat")
            lon_name = normalized.get("longitude") or normalized.get("lon")

            if lat_name and lon_name:
                for row_no, row in enumerate(dict_reader, start=2):
                    lat_raw = row.get(lat_name)
                    lon_raw = row.get(lon_name)
                    if lat_raw is None or lon_raw is None:
                        raise ValueError(f"Missing lat/lon value at line {row_no}")
                    yield float(lat_raw), float(lon_raw)
                return

    with path.open("r", encoding="utf-8", newline="") as fp:
        reader = csv.reader(fp)
        # The first row is treated as a title/header and always ignored.
        next(reader, None)
        for row_no, row in enumerate(reader, start=2):
            if not row or all(not cell.strip() for cell in row):
                continue
            if len(row) < 2:
                raise ValueError(f"Expected at least 2 columns at line {row_no}")
            yield float(row[0]), float(row[1])


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
    if args.per_hash_limit < 1:
        raise ValueError("--per-hash-limit must be >= 1")

    precision = 7
    processed_count = 0

    with LevelDBStore(args.db_path) as store:
        # 정확한 성능 측정을 위해 100건 정도 읽어본다.
        for _ in range(100):
            try:
                lat, lon = next(_iter_csv_coordinates(args.csv_path))
            except StopIteration:
                break
            _estimate_elevation_from_store(
                store,
                lat=lat,
                lon=lon,
                precision=precision,
                per_hash_limit=args.per_hash_limit,
            )

        started_at = time.perf_counter()
        for lon, lat in _iter_csv_coordinates(args.csv_path):
            result =_estimate_elevation_from_store(
                store,
                lat=lat,
                lon=lon,
                precision=precision,
                per_hash_limit=args.per_hash_limit,
            )
            processed_count += 1
            if processed_count < 200:
                print(f'lat={lat}, lon={lon}, elevation={result["estimated_elevation"]}')

    elapsed_seconds = time.perf_counter() - started_at
    print(f"elapsed_seconds={elapsed_seconds:.3f} count={processed_count}")


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
    p_search_csv.add_argument("--per-hash-limit", type=int, default=5)
    p_search_csv.set_defaults(func=cmd_search_csv)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
