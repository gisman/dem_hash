from __future__ import annotations

from pathlib import Path
from typing import Iterable, List

import orjson

from .geohash_index import GeohashRecord

try:
    import plyvel
except ImportError as exc:  # pragma: no cover
    plyvel = None
    _PLYVEL_IMPORT_ERROR = exc
else:
    _PLYVEL_IMPORT_ERROR = None


class LevelDBStore:
    def __init__(self, db_path: str | Path):
        if plyvel is None:  # pragma: no cover
            raise RuntimeError(
                "plyvel is required for LevelDB support. " "Install dependencies first."
            ) from _PLYVEL_IMPORT_ERROR

        self.db_path = Path(db_path)
        self._db = plyvel.DB(str(self.db_path), create_if_missing=True)

    def close(self) -> None:
        self._db.close()

    def __enter__(self) -> "LevelDBStore":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:  # type: ignore[override]
        self.close()

    @staticmethod
    def _build_key(record: GeohashRecord) -> bytes:
        key = f"{record.geohash}:{record.latitude:.7f}:{record.longitude:.7f}"
        return key.encode("utf-8")

    @staticmethod
    def _build_value(record: GeohashRecord) -> bytes:
        payload = {
            "geohash": record.geohash,
            "latitude": record.latitude,
            "longitude": record.longitude,
            "elevation": record.elevation,
        }
        return orjson.dumps(payload)

    def put_many(
        self, records: Iterable[GeohashRecord], batch_size: int = 100_000
    ) -> int:
        if batch_size < 1:
            raise ValueError("batch_size must be >= 1")

        count = 0
        batch_no = 0
        iterator = iter(records)

        while True:
            wrote = 0

            with self._db.write_batch() as batch:
                for _ in range(batch_size):

                    try:
                        record = next(iterator)
                    except StopIteration:
                        break
                    # if record.elevation != 0.0:
                    batch.put(self._build_key(record), self._build_value(record))
                    wrote += 1

            if wrote == 0:
                break

            batch_no += 1
            count += wrote
            print(
                f"[LevelDB] batch {batch_no}: wrote {wrote} records " f"(total={count})"
            )

        return count

    def search_by_prefix(self, geohash_prefix: str, limit: int = 100) -> List[dict]:
        prefix = geohash_prefix.encode("utf-8")
        results: List[dict] = []

        for _, value in self._db.iterator(prefix=prefix):
            results.append(orjson.loads(value))
            if len(results) >= limit:
                break

        return results
