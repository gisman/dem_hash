# DEM Hash

`data/한반도90m_GRS80.img` 파일은 DEM 원본 데이터다.

이 프로젝트는 아래 기능을 제공한다.
- DEM를 geohash로 인덱싱하여 LevelDB에 저장
- 위도/경도 입력 기반 LevelDB 검색
- 7레벨 geohash 3x3(9셀) 기반 9점 IDW 고도 보정
- CSV 전체 좌표 일괄 고도 산출

## 설치

```bash
uv sync
```

`plyvel` 설치 시 C++ 컴파일러가 필요할 수 있다.

## 1) DEM -> LevelDB 저장

```bash
uv run dem-hash build-leveldb \
	--dem-path data/한반도90m_GRS80.img \
	--db-path data/dem.leveldb \
	--precision 7 \
	--stride 1
```

- `precision`: geohash 길이
- `stride`: 샘플링 간격 (값이 클수록 데이터 수가 줄어 속도가 빨라짐)
- `skip-rows`: 앞에서부터 건너뛸 유효 포인트 수 (중단 후 재시작 용도)

## 2) LevelDB geohash 검색 (lat/lon + IDW)

```bash
uv run dem-hash search-geohash \
	--db-path data/dem.leveldb \
	--lat 37.5665 \
	--lon 126.9780

uv run dem-hash search-geohash \
	--db-path data/dem.leveldb \
	--lat 37.5512 \
	--lon 126.9882

uv run dem-hash search-geohash \
	--db-path data/dem.leveldb \
	--lat 37.3536757 \
	--lon 128.1682600

-- 바다 (해상)
uv run dem-hash search-geohash \
	--db-path data/dem.leveldb \
	--lat  35.7762819 \
	--lon 125.4667766


```

126.9882, 37.5512

- 입력 좌표를 7레벨 geohash로 변환한 뒤, 인접 8셀 포함 총 9셀에서 후보를 조회
- 후보 중 거리 기준 최근접 9점을 사용해 IDW(power=2)로 고도 추정

## 3) CSV 전체 좌표 고도 산출

```bash
uv run dem-hash search-csv \
	--db-path data/dem.leveldb \
	--csv-path data/test.csv
```

- CSV 헤더는 `lat,lon` 또는 `latitude,longitude`를 지원
- 헤더가 없으면 각 행의 첫 두 컬럼을 `lat,lon`으로 간주
- 출력은 총 소요 시간과 처리 건수만 출력


