from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Iterator

import geohash2

from .dem import DEMPoint


@dataclass(frozen=True)
class GeohashRecord:
    geohash: str
    latitude: float
    longitude: float
    elevation: float


def encode_geohash(latitude: float, longitude: float, precision: int) -> str:
    return geohash2.encode(latitude, longitude, precision=precision)


def iter_geohash_records(
    points: Iterable[DEMPoint], precision: int
) -> Iterator[GeohashRecord]:
    for point in points:
        yield GeohashRecord(
            geohash=encode_geohash(
                point.latitude, point.longitude, precision=precision
            ),
            latitude=point.latitude,
            longitude=point.longitude,
            elevation=point.elevation,
        )
