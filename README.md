# DEM Hash

DEM Hash는 DEM 데이터를 geohash 기반으로 인덱싱하고 LevelDB에 저장하여, 위도/경도 입력을 통해 고도를 빠르게 조회할 수 있는 툴이다.

# 기능
- DEM을 geohash로 인덱싱하여 [LevelDB](https://github.com/google/leveldb)(경량 Key-value DB)에 저장
- 위도/경도 입력 기반 LevelDB 검색
- 7레벨 geohash 3x3(9셀) 기반 Kriging 보정
- CSV 데이터 일괄 고도 산출

# 한계

* LevelDB는 Local Key-value DB이므로 다중 사용자 환경에서는 사용하기 어렵다. (Amazon DynamoDB와 같은 클라우드 기반 Key-value DB를 고려할 수 있다.)

# 설치

## 패키지 설치
```bash
uv sync
```

주의; `plyvel` 설치 시 C++ 컴파일러가 필요할 수 있다.

## DEM 파일 설치

브이월드에서 [수치표고모델](https://www.vworld.kr/dtmk/dtmk_ntads_s002.do?dsId=30206)을 다운로드 받아 `data/한반도90m_GRS80.img` 경로에 저장한다.


## DEM -> LevelDB 변환

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

# 실행

## LevelDB geohash 검색 (lat/lon + 보간)

```bash
uv run dem-hash search-geohash \
	--db-path data/dem.leveldb \
	--lat 37.545467355260556 \
	--lon 127.07256862776254

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

## CSV 파일 일괄 고도 산출

```bash
uv run dem-hash search-csv \
	--db-path data/dem.leveldb \
	--csv-path data/pnu_centroid_100k.csv \
	--limit 100000
```

- CSV 헤더는 `lat,lon` 또는 `latitude,longitude`를 지원
- 헤더가 없으면 각 행의 첫 두 컬럼을 `lat,lon`으로 간주
- 출력은 총 소요 시간과 처리 건수만 출력

