from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import rasterio
from pyproj import Transformer
from rasterio.crs import CRS


@dataclass(frozen=True)
class DEMPoint:
    latitude: float
    longitude: float
    elevation: float


def iter_dem_points(
    dem_path: str | Path,
    stride: int = 1,
    source_epsg_if_missing: int = 5186,
    skip_rows: int = 0,
) -> Iterator[DEMPoint]:
    if stride < 1:
        raise ValueError("stride must be >= 1")
    if skip_rows < 0:
        raise ValueError("skip_rows must be >= 0")

    dem_path = Path(dem_path)
    with rasterio.open(dem_path) as dataset:
        if dataset.crs is None:
            source_crs = CRS.from_epsg(source_epsg_if_missing)
            print(
                "[DEM] CRS is missing; "
                f"assuming EPSG:{source_epsg_if_missing} for coordinate transform"
            )
        else:
            source_crs = dataset.crs

        transformer: Transformer | None
        if source_crs.to_epsg() == 4326:
            transformer = None
        else:
            transformer = Transformer.from_crs(source_crs, "EPSG:4326", always_xy=True)

        band = dataset.read(1, masked=True)
        nodata = dataset.nodata
        skipped = 0

        for row in range(0, dataset.height, stride):
            for col in range(0, dataset.width, stride):
                cell = band[row, col]
                if getattr(cell, "mask", False):
                    continue

                elevation = float(cell)
                if nodata is not None and elevation == nodata:
                    continue

                if skipped < skip_rows:
                    skipped += 1
                    continue

                x, y = dataset.xy(row, col)
                if transformer is None:
                    lon, lat = float(x), float(y)
                else:
                    lon, lat = transformer.transform(x, y)
                    lon, lat = float(lon), float(lat)

                yield DEMPoint(latitude=lat, longitude=lon, elevation=elevation)
