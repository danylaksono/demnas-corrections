#!/usr/bin/env python3
"""
demnas_dtm.py - DEMNAS (surface model) to DTM (ground) correction and smoothing.

Pipeline
  1. Mosaic DEMNAS tiles, clip to an AOI (e.g. a KHG boundary), reproject to UTM.
  2. Fetch vegetation corrections (default: cloud-native COGs, read by window):
       - land cover: ESA WorldCover 10 m (public COGs on AWS)
       - canopy height: Meta/WRI global canopy height (public COGs on AWS)
     or use local rasters instead.
  3. Vegetation correction:  DTM = DEM - k[landcover] * CHM - bias[landcover]
  4. Mask cells that should not shape the ground surface (water, built-up,
     optional vector masks such as canals or settlements), remove spikes.
  5. Smooth to a dome-scale surface (gaussian, median, polynomial trend,
     or a combination), fill gaps, write GeoTIFFs and a run report.

  6. Ground truthing (optional, on by default): fetch ICESat-2 ATL08 (and GEDI L2A)
     via SlideRule or read local points, convert to the DEM datum, fit k and bias
     per land cover class, apply them, and report errors on held-out points.

Usage
  python demnas_dtm.py --dem "demnas/*.tif" --aoi khg.gpkg --out out_khg01
  python demnas_dtm.py --config config.yaml
  python demnas_dtm.py --dem dem.tif --out out --no-download \
         --landcover lc.tif --chm chm.tif --ground local --ground-points pts.csv
  python demnas_dtm.py --dem dem.tif --out out --ground none   # skip calibration
"""
from __future__ import annotations

import argparse
import copy
import glob
import json
import logging
import math
import os
import sys
import time
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.fill import fillnodata
from rasterio.merge import merge
from rasterio.transform import from_origin
from rasterio.warp import reproject, transform_bounds, calculate_default_transform
from scipy import ndimage

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ground_truth as gt  # noqa: E402

log = logging.getLogger("demnas_dtm")

# --------------------------------------------------------------------------------------
# Default configuration (every value can be overridden by YAML or CLI)
# --------------------------------------------------------------------------------------
DEFAULTS: dict = {
    "dem": [],                    # list of paths / globs to DEMNAS tiles
    "aoi": None,                  # vector file (KHG boundary); None = full DEM extent
    "aoi_buffer_m": 500,          # buffer around AOI so smoothing has context at the edges
    "out": "out_dtm",
    "target_crs": "auto-utm",     # "auto-utm" or any CRS string, e.g. "EPSG:32749"
    "target_res_m": 8.0,          # DEMNAS native is ~8 m (0.27 arc-second)
    "dem_valid_range_m": [-50, 4000],   # values outside are treated as nodata

    "corrections": {
        "mode": "download",       # "download" | "local" | "none"
        "landcover_path": None,   # used when mode = local
        "chm_path": None,         # used when mode = local
        "worldcover_url": (
            "https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/"
            "ESA_WorldCover_10m_2021_v200_{tile}_Map.tif"),
        "chm_url": (
            "https://dataforgood-fb-data.s3.amazonaws.com/forests/v1/"
            "alsgedi_global_v6_float/chm/{quadkey}.tif"),
        "chm_quadkey_zoom": 9,
        "chm_min_m": 2.0,         # canopy below this is treated as 0 (noise)
        "chm_max_m": 60.0,        # clip implausible heights
        "chm_fallback_m": {"10": 20.0, "95": 12.0},  # used per class if no CHM tile is found
    },

    # Penetration factor k per ESA WorldCover class: share of canopy height present in the DEM.
    # PLACEHOLDERS - calibrate against ICESat-2/GEDI before relying on them.
    "k_by_class": {
        "10": 0.8,   # tree cover
        "20": 0.5,   # shrubland
        "30": 0.2,   # grassland
        "40": 0.2,   # cropland
        "50": 0.0,   # built-up (masked anyway)
        "60": 0.0,   # bare / sparse vegetation
        "70": 0.0,   # snow and ice
        "80": 0.0,   # permanent water bodies (masked anyway)
        "90": 0.4,   # herbaceous wetland
        "95": 0.8,   # mangroves
        "100": 0.0,  # moss and lichen
    },
    "k_default": 0.5,             # for classes not listed above
    "bias_by_class": {},          # constant offset per class (m), e.g. {"10": 0.5}
    "bias_default": 0.0,

    "mask": {
        "landcover_classes": [80, 50],   # water, built-up: excluded from smoothing input
        "vectors": [],                   # e.g. [{"path": "canals.gpkg", "buffer_m": 15}]
        "despike": True,
        "despike_window_m": 400,         # local reference for spike detection
        "despike_threshold_m": 2.5,      # |DTM - local median| above this is removed
    },

    "smoothing": {
        # Steps applied in order. Available methods:
        #   gaussian   (sigma_m)
        #   median     (window_m)           - decimated, robust to remaining outliers
        #   polynomial (degree, sample, huber_m) - robust dome-scale trend surface
        "steps": [
            {"method": "median", "window_m": 300},
            {"method": "gaussian", "sigma_m": 150},
        ],
        "fill_max_search_px": 200,       # gap filling after masking
    },

    "ground": copy.deepcopy(gt.GROUND_DEFAULTS),

    "outputs": {
        "write_intermediates": True,     # landcover, chm, k raster, corrected (unsmoothed) DTM
        "compress": "deflate",
    },
}


# --------------------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------------------
def deep_update(base: dict, upd: dict) -> dict:
    for k, v in (upd or {}).items():
        if isinstance(v, dict) and isinstance(base.get(k), dict):
            deep_update(base[k], v)
        else:
            base[k] = v
    return base


def load_config(args: argparse.Namespace) -> dict:
    cfg = copy.deepcopy(DEFAULTS)
    if args.config:
        import yaml
        with open(args.config, "r", encoding="utf-8") as f:
            deep_update(cfg, yaml.safe_load(f) or {})
    if args.dem:
        cfg["dem"] = args.dem
    if args.aoi:
        cfg["aoi"] = args.aoi
    if args.out:
        cfg["out"] = args.out
    if args.res:
        cfg["target_res_m"] = args.res
    if args.crs:
        cfg["target_crs"] = args.crs
    if args.no_download:
        cfg["corrections"]["mode"] = "local" if (args.landcover or args.chm) else "none"
    if args.landcover:
        cfg["corrections"]["landcover_path"] = args.landcover
    if args.chm:
        cfg["corrections"]["chm_path"] = args.chm
    if args.ground:
        cfg["ground"]["source"] = args.ground
    if args.ground_points:
        cfg["ground"]["source"] = "local"
        cfg["ground"]["local"]["path"] = args.ground_points
    if args.datasets:
        cfg["ground"]["datasets"] = args.datasets
    if isinstance(cfg["dem"], str):
        cfg["dem"] = [cfg["dem"]]
    return cfg


@dataclass
class Grid:
    crs: rasterio.crs.CRS
    transform: rasterio.Affine
    width: int
    height: int

    @property
    def res(self) -> float:
        return self.transform.a

    @property
    def bounds(self):
        left, top = self.transform.c, self.transform.f
        return (left, top - self.height * self.res, left + self.width * self.res, top)

    def profile(self, dtype="float32", nodata=np.nan, compress="deflate") -> dict:
        return dict(driver="GTiff", crs=self.crs, transform=self.transform, width=self.width,
                    height=self.height, count=1, dtype=dtype, nodata=nodata, tiled=True,
                    blockxsize=512, blockysize=512, compress=compress, BIGTIFF="IF_SAFER")


def auto_utm(lon: float, lat: float) -> str:
    zone = int((lon + 180) // 6) + 1
    return f"EPSG:{(32600 if lat >= 0 else 32700) + zone}"


def write_raster(path: Path, arr: np.ndarray, grid: Grid, cfg: dict, dtype="float32", nodata=np.nan):
    prof = grid.profile(dtype=dtype, nodata=nodata, compress=cfg["outputs"]["compress"])
    if dtype != "float32" and isinstance(nodata, float) and math.isnan(nodata):
        prof["nodata"] = 0
    with rasterio.open(path, "w", **prof) as dst:
        dst.write(arr.astype(dtype), 1)
    log.info("wrote %s", path)


def gdal_cloud_env():
    """GDAL settings for efficient, unsigned windowed reads of public COGs."""
    return rasterio.Env(
        GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
        CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif,.tiff",
        GDAL_HTTP_MULTIRANGE="YES",
        GDAL_HTTP_MERGE_CONSECUTIVE_RANGES="YES",
        GDAL_HTTP_MAX_RETRY="4",
        GDAL_HTTP_RETRY_DELAY="2",
        VSI_CACHE="TRUE",
        AWS_NO_SIGN_REQUEST="YES",
    )


# --------------------------------------------------------------------------------------
# 1. DEM: mosaic, clip, reproject
# --------------------------------------------------------------------------------------
def read_aoi(path: str | None, buffer_m: float):
    """Return AOI geometry (shapely, EPSG:4326) or None."""
    if not path:
        return None
    import geopandas as gpd
    gdf = gpd.read_file(path)
    if gdf.crs is None:
        raise ValueError(f"AOI {path} has no CRS")
    union = (lambda g: g.union_all() if hasattr(g, "union_all") else g.unary_union)
    utm = auto_utm(*union(gdf.to_crs(4326)).centroid.coords[0])
    geom = union(gdf.to_crs(utm)).buffer(buffer_m)
    return gpd.GeoSeries([geom], crs=utm).to_crs(4326).iloc[0]


def load_dem(cfg: dict) -> tuple[np.ndarray, Grid, object]:
    paths: list[str] = []
    for p in cfg["dem"]:
        hits = sorted(glob.glob(p))
        paths.extend(hits if hits else [p])
    if not paths:
        raise SystemExit("No DEM input given (use --dem or the 'dem' config key).")
    log.info("DEM tiles: %d", len(paths))

    srcs = [rasterio.open(p) for p in paths]
    src_crs = srcs[0].crs
    aoi = read_aoi(cfg["aoi"], cfg["aoi_buffer_m"])
    bounds = None
    if aoi is not None:
        bounds = transform_bounds("EPSG:4326", src_crs, *aoi.bounds)
    src_nodata = srcs[0].nodata if srcs[0].nodata is not None else -99999.0
    mosaic, mtrans = merge(srcs, bounds=bounds, nodata=src_nodata, dtype="float32")
    for s in srcs:
        s.close()
    dem = mosaic[0].astype("float32")
    dem[dem == src_nodata] = np.nan
    lo, hi = cfg["dem_valid_range_m"]
    dem[(dem < lo) | (dem > hi)] = np.nan

    # target grid
    h, w = dem.shape
    left, top = mtrans.c, mtrans.f
    right, bottom = left + w * mtrans.a, top + h * mtrans.e
    if cfg["target_crs"] == "auto-utm":
        lon0, lat0, lon1, lat1 = transform_bounds(src_crs, "EPSG:4326", left, bottom, right, top)
        dst_crs = rasterio.crs.CRS.from_string(auto_utm((lon0 + lon1) / 2, (lat0 + lat1) / 2))
    else:
        dst_crs = rasterio.crs.CRS.from_string(cfg["target_crs"])
    res = float(cfg["target_res_m"])
    t, tw, th = calculate_default_transform(src_crs, dst_crs, w, h, left, bottom, right, top, resolution=res)
    grid = Grid(dst_crs, from_origin(t.c, t.f, res, res), tw, th)

    out = np.full((grid.height, grid.width), np.nan, dtype="float32")
    reproject(dem, out, src_transform=mtrans, src_crs=src_crs, src_nodata=np.nan,
              dst_transform=grid.transform, dst_crs=grid.crs, dst_nodata=np.nan,
              resampling=Resampling.bilinear)

    aoi_mask = None
    if aoi is not None:
        from rasterio.features import geometry_mask
        import geopandas as gpd
        aoi_proj = gpd.GeoSeries([aoi], crs=4326).to_crs(grid.crs).iloc[0]
        aoi_mask = geometry_mask([aoi_proj], out_shape=out.shape, transform=grid.transform, invert=True)
        out[~aoi_mask] = np.nan
    log.info("DEM grid: %s, %d x %d px at %.1f m", grid.crs.to_string(), grid.width, grid.height, res)
    if aoi is None:
        from shapely.geometry import box
        aoi = box(*transform_bounds(grid.crs, "EPSG:4326", *grid.bounds, densify_pts=21))
    return out, grid, aoi_mask, aoi


# --------------------------------------------------------------------------------------
# 2. Corrections: land cover + canopy height
# --------------------------------------------------------------------------------------
def grid_bounds_ll(grid: Grid):
    return transform_bounds(grid.crs, "EPSG:4326", *grid.bounds, densify_pts=21)


def worldcover_tiles(bounds_ll) -> list[str]:
    lon0, lat0, lon1, lat1 = bounds_ll
    tiles = []
    for lat in range(int(math.floor(lat0 / 3) * 3), int(math.floor(lat1 / 3) * 3) + 1, 3):
        for lon in range(int(math.floor(lon0 / 3) * 3), int(math.floor(lon1 / 3) * 3) + 1, 3):
            ns = "N" if lat >= 0 else "S"
            ew = "E" if lon >= 0 else "W"
            tiles.append(f"{ns}{abs(lat):02d}{ew}{abs(lon):03d}")
    return tiles


def chm_quadkeys(bounds_ll, zoom: int) -> list[str]:
    import mercantile
    lon0, lat0, lon1, lat1 = bounds_ll
    return [mercantile.quadkey(t) for t in mercantile.tiles(lon0, lat0, lon1, lat1, zooms=[zoom])]


def pick_overview(src, target_res_src_units: float) -> int | None:
    """Largest overview whose resolution is still <= half the target resolution."""
    ovr = src.overviews(1)
    best = None
    for i, f in enumerate(ovr):
        if src.res[0] * f <= target_res_src_units / 2:
            best = i
    return best


def warp_into(path: str, dst: np.ndarray, grid: Grid, resampling: Resampling,
              remote: bool, nodata_override=None) -> bool:
    """Warp a (possibly remote) raster into dst where dst is still empty. Returns success."""
    vpath = f"/vsicurl/{path}" if remote and path.startswith("http") else path
    try:
        with rasterio.open(vpath) as probe:
            src_crs = probe.crs
            # target resolution expressed in source units (approximate for geographic CRS)
            tgt = grid.res if not src_crs.is_geographic else grid.res / 111320.0
            ovl = pick_overview(probe, tgt) if remote else None
        with rasterio.open(vpath, overview_level=ovl) if ovl is not None else rasterio.open(vpath) as src:
            nod = nodata_override if nodata_override is not None else src.nodata
            tmp = np.full(dst.shape, np.nan, dtype="float32")
            reproject(rasterio.band(src, 1), tmp, src_nodata=nod,
                      dst_transform=grid.transform, dst_crs=grid.crs, dst_nodata=np.nan,
                      resampling=resampling)
        fill = np.isnan(dst) & ~np.isnan(tmp)
        dst[fill] = tmp[fill]
        log.info("  read %s%s", path, f" (overview {ovl})" if ovl is not None else "")
        return True
    except Exception as e:  # missing tile (ocean, no data) or network error
        log.warning("  could not read %s: %s", path, e)
        return False


def get_corrections(cfg: dict, grid: Grid) -> tuple[np.ndarray | None, np.ndarray | None, dict]:
    c = cfg["corrections"]
    meta = {"mode": c["mode"], "landcover_sources": [], "chm_sources": []}
    shape = (grid.height, grid.width)
    if c["mode"] == "none":
        log.info("Corrections: none (DEM used as is)")
        return None, None, meta

    lc = np.full(shape, np.nan, dtype="float32")
    chm = np.full(shape, np.nan, dtype="float32")
    if c["mode"] == "local":
        if c.get("landcover_path"):
            warp_into(c["landcover_path"], lc, grid, Resampling.nearest, remote=False)
            meta["landcover_sources"].append(c["landcover_path"])
        if c.get("chm_path"):
            warp_into(c["chm_path"], chm, grid, Resampling.average, remote=False)
            meta["chm_sources"].append(c["chm_path"])
    elif c["mode"] == "download":
        bll = grid_bounds_ll(grid)
        with gdal_cloud_env():
            log.info("Fetching ESA WorldCover tiles ...")
            for t in worldcover_tiles(bll):
                url = c["worldcover_url"].format(tile=t)
                if warp_into(url, lc, grid, Resampling.nearest, remote=True, nodata_override=0):
                    meta["landcover_sources"].append(url)
            log.info("Fetching canopy height tiles ...")
            for qk in chm_quadkeys(bll, int(c["chm_quadkey_zoom"])):
                url = c["chm_url"].format(quadkey=qk)
                if warp_into(url, chm, grid, Resampling.average, remote=True):
                    meta["chm_sources"].append(url)
    else:
        raise ValueError(f"unknown corrections.mode: {c['mode']}")

    if not meta["landcover_sources"]:
        log.warning("No land cover available: k_default applied everywhere, no class masking.")
        lc = None
    if meta["chm_sources"]:
        chm = np.clip(chm, 0, float(c["chm_max_m"]))
        chm[chm < float(c["chm_min_m"])] = 0.0
    else:
        log.warning("No canopy height available.")
        chm = None
    # fill CHM gaps (or the whole CHM) with per-class fallback heights
    fb = {int(k): float(v) for k, v in (c.get("chm_fallback_m") or {}).items()}
    if lc is not None and fb:
        if chm is None:
            chm = np.zeros(shape, dtype="float32")
            meta["chm_fallback_used"] = "all"
        gaps = np.isnan(chm)
        for cls, h in fb.items():
            chm[gaps & (lc == cls)] = h
        chm[np.isnan(chm)] = 0.0
    elif chm is not None:
        chm[np.isnan(chm)] = 0.0
    return lc, chm, meta


def class_lookup(lc: np.ndarray | None, table: dict, default: float, shape) -> np.ndarray:
    out = np.full(shape, float(default), dtype="float32")
    if lc is None:
        return out
    for k, v in table.items():
        out[lc == int(k)] = float(v)
    return out


# --------------------------------------------------------------------------------------
# 3-5. Correction, masking, smoothing
# --------------------------------------------------------------------------------------
def nan_gaussian(z: np.ndarray, sigma_px: float) -> np.ndarray:
    valid = np.isfinite(z)
    num = ndimage.gaussian_filter(np.where(valid, z, 0.0).astype("float64"), sigma_px, mode="nearest")
    den = ndimage.gaussian_filter(valid.astype("float64"), sigma_px, mode="nearest")
    out = np.where(den > 0.05, num / np.maximum(den, 1e-12), np.nan)
    return out.astype("float32")


def coarse_median(z: np.ndarray, grid: Grid, window_m: float) -> np.ndarray:
    """Robust large-window median: block-median to a coarse grid, 3x3 median, bilinear back."""
    cell_m = max(grid.res, window_m / 3.0)
    f = max(1, int(round(cell_m / grid.res)))
    h, w = z.shape
    H, W = math.ceil(h / f), math.ceil(w / f)
    pad = np.full((H * f, W * f), np.nan, dtype="float32")
    pad[:h, :w] = z
    blocks = pad.reshape(H, f, W, f).transpose(0, 2, 1, 3).reshape(H, W, f * f)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        coarse = np.nanmedian(blocks, axis=2).astype("float32")
    cmask = np.isfinite(coarse)
    if cmask.any() and (~cmask).any():
        coarse = fillnodata(coarse, mask=cmask.astype("uint8"), max_search_distance=50)
    coarse = ndimage.median_filter(np.nan_to_num(coarse, nan=np.nanmean(coarse)), size=3, mode="nearest")
    coarse[~np.isfinite(coarse)] = np.nan
    ctrans = rasterio.Affine(grid.res * f, 0, grid.transform.c, 0, -grid.res * f, grid.transform.f)
    out = np.full(z.shape, np.nan, dtype="float32")
    reproject(coarse.astype("float32"), out, src_transform=ctrans, src_crs=grid.crs,
              dst_transform=grid.transform, dst_crs=grid.crs, resampling=Resampling.bilinear,
              src_nodata=np.nan, dst_nodata=np.nan)
    return out


def poly_trend(z: np.ndarray, degree: int = 3, sample: int = 200_000, huber_m: float = 1.0,
               iters: int = 6, seed: int = 0) -> np.ndarray:
    """Robust polynomial trend surface (IRLS with Huber weights) - dome-scale ground form."""
    h, w = z.shape
    rows, cols = np.nonzero(np.isfinite(z))
    if rows.size < 50:
        raise ValueError("too few valid cells for polynomial trend")
    rng = np.random.default_rng(seed)
    if rows.size > sample:
        sel = rng.choice(rows.size, sample, replace=False)
        rows, cols = rows[sel], cols[sel]
    x = cols / max(w - 1, 1) * 2 - 1
    y = rows / max(h - 1, 1) * 2 - 1
    terms = [(i, j) for i in range(degree + 1) for j in range(degree + 1 - i)]
    A = np.stack([x ** i * y ** j for i, j in terms], axis=1)
    b = z[rows, cols].astype("float64")
    wts = np.ones_like(b)
    coef = None
    for _ in range(iters):
        sw = np.sqrt(wts)
        coef, *_ = np.linalg.lstsq(A * sw[:, None], b * sw, rcond=None)
        r = np.abs(b - A @ coef)
        wts = np.where(r <= huber_m, 1.0, huber_m / np.maximum(r, 1e-9))
    out = np.empty((h, w), dtype="float32")
    xs = np.arange(w) / max(w - 1, 1) * 2 - 1
    for r0 in range(0, h, 512):
        r1 = min(h, r0 + 512)
        yy = (np.arange(r0, r1) / max(h - 1, 1) * 2 - 1)[:, None]
        acc = np.zeros((r1 - r0, w))
        for (i, j), c in zip(terms, coef):
            acc += c * (xs[None, :] ** i) * (yy ** j)
        out[r0:r1] = acc
    return out


def vector_mask(cfg: dict, grid: Grid) -> np.ndarray | None:
    items = cfg["mask"].get("vectors") or []
    if not items:
        return None
    import geopandas as gpd
    from rasterio.features import geometry_mask
    m = np.zeros((grid.height, grid.width), dtype=bool)
    for it in items:
        gdf = gpd.read_file(it["path"]).to_crs(grid.crs)
        buf = float(it.get("buffer_m", 0))
        geoms = gdf.geometry.buffer(buf) if buf > 0 else gdf.geometry
        geoms = [g for g in geoms if g is not None and not g.is_empty]
        if geoms:
            m |= geometry_mask(geoms, out_shape=m.shape, transform=grid.transform, invert=True)
        log.info("vector mask %s (buffer %.0f m): %d features", it["path"], buf, len(geoms))
    return m


def smooth(z: np.ndarray, grid: Grid, cfg: dict) -> np.ndarray:
    s = z.copy()
    for step in cfg["smoothing"]["steps"]:
        m = step["method"]
        t0 = time.time()
        if m == "gaussian":
            s = nan_gaussian(s, float(step["sigma_m"]) / grid.res)
        elif m == "median":
            s = coarse_median(s, grid, float(step["window_m"]))
        elif m == "polynomial":
            s = poly_trend(s, int(step.get("degree", 3)), int(step.get("sample", 200_000)),
                           float(step.get("huber_m", 1.0)))
        else:
            raise ValueError(f"unknown smoothing method: {m}")
        log.info("smoothing step %s done (%.1f s)", m, time.time() - t0)
    return s


def fill_gaps(z: np.ndarray, max_px: int) -> np.ndarray:
    valid = np.isfinite(z)
    if valid.all() or not valid.any():
        return z
    filled = fillnodata(np.where(valid, z, 0).astype("float32"), mask=valid.astype("uint8"),
                        max_search_distance=max_px, smoothing_iterations=0)
    filled[~valid & (filled == 0)] = np.nan
    return filled


# --------------------------------------------------------------------------------------
# 6. Ground truthing
# --------------------------------------------------------------------------------------
def ground_calibration(cfg, grid, dem, lc, chm, vm, aoi_ll, out):
    g = cfg["ground"]
    pts = gt.get_points(aoi_ll, g, out)
    if len(pts) == 0:
        log.warning("no ground points available; prior k values kept")
        return None, None, None
    x, y = gt.to_grid_xy(pts, grid)
    em = gt.edge_mask(lc, float(g["exclude"]["edge_buffer_m"]), grid.res)
    s = gt.sample_arrays(x, y, grid, {"dem": dem, "chm": chm, "lc": lc,
                                      "edge": None if em is None else em.astype("float32"),
                                      "vmask": None if vm is None else vm.astype("float32")})
    pts = gt.to_dem_datum(pts, g, s["lc"])
    df = pts.drop(columns="geometry").copy()
    df["x"], df["y"] = x, y
    for k_ in ("dem", "chm", "lc", "edge", "vmask"):
        df[k_] = s[k_]
    df["chm"] = df["chm"].fillna(0.0)

    reason = np.full(len(df), "", dtype=object)
    reason[~s["inside"]] = "outside grid"
    reason[(reason == "") & ~np.isfinite(df["dem"].to_numpy())] = "no DEM"
    reason[(reason == "") & ~np.isfinite(df["H"].to_numpy())] = "no height"
    if lc is not None:
        reason[(reason == "") & df["lc"].isin([int(c) for c in g["exclude"]["landcover_classes"]]).to_numpy()] = "excluded class"
    reason[(reason == "") & (df["edge"].to_numpy() == 1)] = "near land-cover edge"
    if g["exclude"]["use_vector_masks"]:
        reason[(reason == "") & (df["vmask"].to_numpy() == 1)] = "inside vector mask"
    gross = np.abs(df["dem"].to_numpy() - df["H"].to_numpy()) > float(g["exclude"]["max_abs_residual_m"])
    reason[(reason == "") & gross] = "gross outlier"
    df["exclude_reason"] = reason
    use = reason == ""
    log.info("ground points usable: %d of %d", int(use.sum()), len(df))
    log.info("exclusions: %s", df.loc[~use, "exclude_reason"].value_counts().to_dict())
    diag = gt.datum_diagnostics(df, g["diagnostics"], float(g["fit"]["bias_from_chm_below_m"]))
    with open(out / "datum_diagnostics.json", "w", encoding="utf-8") as f:
        json.dump(diag, f, indent=2)

    df["split"] = ""
    if use.any():
        val = gt.spatial_split(df["x"].to_numpy()[use], df["y"].to_numpy()[use], g["split"])
        df.loc[use, "split"] = np.where(val, "validation", "calibration")
    cal = df[df["split"] == "calibration"]
    if len(cal) < int(g["fit"]["min_points_per_class"]):
        log.warning("too few calibration points (%d); prior k values kept", len(cal))
        calib = None
    else:
        log.info("fitting k and bias on %d calibration points", len(cal))
        calib = gt.fit_k_bias(cal, cfg["k_by_class"], cfg["k_default"], g["fit"])
        calib["n_validation"] = int((df["split"] == "validation").sum())
    import geopandas as gpd
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df["x"], df["y"]), crs=grid.crs)
    if calib:
        with open(out / "calibration.json", "w", encoding="utf-8") as f:
            json.dump(calib, f, indent=2)
    return calib, gdf, diag


def ground_validation(gpts, grid, surfaces: dict, out) -> list | None:
    s = gt.sample_arrays(gpts["x"].to_numpy(), gpts["y"].to_numpy(), grid, surfaces)
    for name in surfaces:
        gpts["err_" + name] = s[name] - gpts["H"].to_numpy()
    gpts.to_file(out / "ground_points.gpkg", driver="GPKG")
    log.info("wrote %s", out / "ground_points.gpkg")
    val = gpts[gpts["split"] == "validation"]
    if len(val) == 0:
        log.warning("no validation points; skipping validation table")
        return None
    vsel = (gpts["split"] == "validation").to_numpy()
    tab = gt.error_table(val["H"].to_numpy(), {n: s[n][vsel] for n in surfaces},
                         val["lc"].to_numpy(), gt.WORLDCOVER_NAMES)
    tab.to_csv(out / "validation.csv", index=False)
    log.info("wrote %s", out / "validation.csv")
    allrows = tab[tab["class"] == "all"]
    for _, r in allrows.iterrows():
        log.info("  validation %-15s n=%5d bias=%6.2f  RMSE=%6.2f  NMAD=%5.2f",
                 r["surface"], r["n"], r["bias_m"], r["rmse_m"], r["nmad_m"])
    return allrows.to_dict(orient="records")


# --------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------
def run(cfg: dict) -> dict:
    out = Path(cfg["out"])
    out.mkdir(parents=True, exist_ok=True)
    t_start = time.time()

    dem, grid, aoi_mask, aoi_ll = load_dem(cfg)
    shape = dem.shape

    lc, chm, cmeta = get_corrections(cfg, grid)
    vm = vector_mask(cfg, grid)

    def correct(k_table, k_def, b_table, b_def):
        k_ = class_lookup(lc, k_table, k_def, shape)
        b_ = class_lookup(lc, b_table, b_def, shape)
        corr = (k_ * chm if chm is not None else 0.0) + b_
        return (dem - corr).astype("float32"), k_, corr

    dtm_prior, _, _ = correct(cfg["k_by_class"], cfg["k_default"], cfg["bias_by_class"], cfg["bias_default"])

    # ground truthing: fit k and bias, then apply
    calib, gpts, diag = None, None, None
    if cfg["ground"]["enabled"] and cfg["ground"]["source"] != "none":
        try:
            calib, gpts, diag = ground_calibration(cfg, grid, dem, lc, chm, vm, aoi_ll, out)
        except Exception as e:
            if cfg["ground"]["fail_on_error"]:
                raise
            log.error("ground truthing failed, continuing with prior k values: %s", e)
    k_table, b_table = dict(cfg["k_by_class"]), dict(cfg["bias_by_class"])
    b_def = cfg["bias_default"]
    if calib:
        for cls, e in calib["classes"].items():
            k_table[cls] = e["k"]
            b_table[cls] = e["bias"]
        b_def = calib["bias_global"]["value"]
    dtm, k, correction = correct(k_table, cfg["k_default"], b_table, b_def)

    # mask cells that should not shape the ground surface
    excl = np.zeros(shape, dtype=bool)
    if lc is not None:
        for cls in cfg["mask"]["landcover_classes"]:
            excl |= lc == int(cls)
    if vm is not None:
        excl |= vm
    base = np.where(excl, np.nan, dtm)

    spikes = np.zeros(shape, dtype=bool)
    if cfg["mask"]["despike"]:
        ref = coarse_median(base, grid, float(cfg["mask"]["despike_window_m"]))
        spikes = np.isfinite(base) & (np.abs(base - ref) > float(cfg["mask"]["despike_threshold_m"]))
        base[spikes] = np.nan
        log.info("despike: removed %.2f%% of cells", 100 * spikes.mean())

    smoothed = smooth(base, grid, cfg)
    smoothed = fill_gaps(smoothed, int(cfg["smoothing"]["fill_max_search_px"]))
    if aoi_mask is not None:
        smoothed[~aoi_mask] = np.nan
    residual = dtm - smoothed

    validation = None
    if gpts is not None and len(gpts):
        validation = ground_validation(gpts, grid, {"dem_raw": dem, "dtm_prior_k": dtm_prior,
                                                    "dtm_calibrated": dtm, "dtm_smooth": smoothed},
                                       out)

    # outputs
    write_raster(out / "dtm_smooth.tif", smoothed, grid, cfg)
    write_raster(out / "dtm_residual.tif", residual, grid, cfg)
    if cfg["outputs"]["write_intermediates"]:
        write_raster(out / "dem_utm.tif", dem, grid, cfg)
        write_raster(out / "dtm_corrected.tif", dtm, grid, cfg)
        write_raster(out / "correction_m.tif", np.where(np.isfinite(dem), correction, np.nan), grid, cfg)
        write_raster(out / "k_factor.tif", k, grid, cfg)
        write_raster(out / "excluded_mask.tif", (excl | spikes).astype("uint8"), grid, cfg, "uint8", 255)
        if lc is not None:
            write_raster(out / "landcover.tif", np.nan_to_num(lc, nan=0), grid, cfg, "uint8", 0)
        if chm is not None:
            write_raster(out / "canopy_height.tif", chm, grid, cfg)

    def stats(a):
        a = a[np.isfinite(a)]
        if a.size == 0:
            return None
        return {"min": float(a.min()), "p05": float(np.percentile(a, 5)), "median": float(np.median(a)),
                "p95": float(np.percentile(a, 95)), "max": float(a.max())}

    report = {
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),
        "runtime_s": round(time.time() - t_start, 1),
        "grid": {"crs": grid.crs.to_string(), "res_m": grid.res, "width": grid.width,
                 "height": grid.height, "bounds": grid.bounds},
        "corrections": cmeta,
        "excluded_share": float((excl | spikes).mean()),
        "spike_share": float(spikes.mean()),
        "stats": {"dem": stats(dem), "correction_m": stats(np.where(np.isfinite(dem), correction, np.nan)),
                  "dtm_smooth": stats(smoothed), "residual": stats(residual)},
        "calibration": calib,
        "datum_diagnostics": diag,
        "validation_all_classes": validation,
        "k_used": {k_: float(v) for k_, v in k_table.items()},
        "config": cfg,
        "note": ("k and bias fitted from ground points" if calib else
                 "k and bias are uncalibrated placeholders (no ground calibration applied)"),
    }
    with open(out / "run_report.json", "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, default=str)
    log.info("done in %.1f s -> %s", time.time() - t_start, out)
    return report


def main(argv=None):
    ap = argparse.ArgumentParser(description="DEMNAS to DTM vegetation correction and smoothing")
    ap.add_argument("--config", help="YAML config (overrides defaults; CLI overrides config)")
    ap.add_argument("--dem", nargs="+", help="DEMNAS tile paths or globs")
    ap.add_argument("--aoi", help="AOI vector (e.g. KHG boundary)")
    ap.add_argument("--out", help="output folder")
    ap.add_argument("--res", type=float, help="target resolution in metres")
    ap.add_argument("--crs", help="target CRS, default auto-utm")
    ap.add_argument("--no-download", action="store_true", help="do not fetch cloud corrections")
    ap.add_argument("--landcover", help="local land cover raster (ESA WorldCover classes)")
    ap.add_argument("--chm", help="local canopy height raster (m)")
    ap.add_argument("--ground", choices=["sliderule", "local", "none"], help="ground truth source")
    ap.add_argument("--ground-points", help="local ground points (CSV lon,lat,h or vector)")
    ap.add_argument("--datasets", nargs="+", choices=["icesat2_atl08", "gedi_l2a"],
                    help="SlideRule datasets to fetch")
    ap.add_argument("--print-default-config", action="store_true")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")
    if args.print_default_config:
        import yaml
        print(yaml.safe_dump(DEFAULTS, sort_keys=False, allow_unicode=True))
        return
    run(load_config(args))


if __name__ == "__main__":
    main()
