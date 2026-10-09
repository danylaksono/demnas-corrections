"""
ground_truth.py - ground reference points for calibrating the DEMNAS vegetation correction.

Steps (called from demnas_dtm.run):
  1. get_points      ICESat-2 ATL08 (PhoREAL via SlideRule) and optionally GEDI L2A,
                     or a local CSV / vector file of ground heights.
  2. to_dem_datum    ellipsoidal heights -> DEM vertical datum (PROJ geoid, geoid raster,
                     or a constant), plus an optional subsidence epoch shift.
  3. sample_points   DEM, canopy height, land cover and exclusion masks at each point.
  4. spatial_split   calibration / validation split by spatial blocks (not random points,
                     because neighbouring footprints along a track are correlated).
  5. fit_k_bias      robust fit of  DEM - ground = bias + k[class] * canopy_height.
  6. validate        error statistics per land cover class for each surface.

Notes
  - SlideRule's public service is used through its Python client (pip install sliderule).
    Field names follow the client at the time of writing (h_te_median for ATL08-PhoREAL,
    elevation_lm for GEDI L2A); if the service changes, export points to CSV and use
    ground.source: local.
  - All heights from ICESat-2 and GEDI are ellipsoidal (ITRF2014 / EPSG:7912).
"""
from __future__ import annotations

import logging
import math

import numpy as np
import pandas as pd

log = logging.getLogger("demnas_dtm.ground")

GROUND_DEFAULTS: dict = {
    "enabled": True,
    "fail_on_error": False,            # False: log the error and continue with prior k values
    "source": "sliderule",             # "sliderule" | "local" | "none"
    "datasets": ["icesat2_atl08"],     # add "gedi_l2a" for denser coverage
    "t0": "2019-01-01T00:00:00Z",
    "t1": "2025-12-31T23:59:59Z",
    "cache": True,                     # reuse <out>/ground_points_raw.gpkg if it exists
    "simplify_aoi_deg": 0.001,         # simplify AOI polygon before sending to SlideRule
    "tile_deg": 0.5,                   # split larger AOIs into tiles of this size (0 = no tiling)
    "icesat2": {
        "len": 100.0,                  # segment length (m); ATL08 standard is 100 m
        "res": 100.0,                  # step between segments (m)
        "atl08_class": ["atl08_ground", "atl08_canopy", "atl08_top_of_canopy"],
        "phoreal": {"binsize": 1.0, "geoloc": "center", "use_abs_h": False,
                    "send_waveform": False, "above_classifier": False},
        "night_only": False,           # solar_elevation < 0 only (better SNR, fewer points)
        "min_ground_photons": 10,      # applied if the service returns a photon-count column
    },
    "gedi": {
        "degrade_flag": 0,
        "l2_quality_flag": 1,
        "min_sensitivity": 0.95,       # applied if a sensitivity column is returned
    },
    "local": {
        "path": None,                  # CSV (lon, lat, h, [source], [time]) or vector file
        "lon_col": "lon", "lat_col": "lat", "h_col": "h",
        "time_col": "time", "source_col": "source",
        "height_is_ellipsoidal": True,
    },
    "datum": {
        "mode": "pyproj",              # "pyproj" | "raster" | "constant" | "none"
        "target": "EPSG:4326+3855",    # WGS84 + EGM2008 heights (check DEMNAS metadata)
        "raster_path": None,           # geoid undulation N (m) raster, e.g. INAGeoid2020
        "constant_m": 0.0,             # used when mode = constant
    },
    "epoch": {
        "dem_year": None,              # e.g. 2016; None = no subsidence shift
        "subsidence_m_per_yr_by_class": {},   # e.g. {"40": 0.04, "20": 0.03}
        "subsidence_default_m_per_yr": 0.0,
    },
    "exclude": {
        "landcover_classes": [80, 50],
        "edge_buffer_m": 50,           # drop points within this distance of a land-cover change
        "use_vector_masks": True,      # also drop points inside mask.vectors (e.g. canal buffers)
        "max_abs_residual_m": 30.0,    # gross outliers (cloud returns, misclassified photons)
    },
    "split": {"validation_fraction": 0.3, "block_size_m": 2000, "seed": 42},
    "fit": {
        "min_points_per_class": 30,
        "bias_from_chm_below_m": 2.0,  # points with canopy below this define the global bias
        "fit_bias_per_class": False,
        "joint_intercept": False,      # fit DEM - ground = a[class] + k[class] * chm on all points of a class
        "k_bounds": [0.0, 1.0],
        "huber_m": 1.5,
    },
    "diagnostics": {
        "warn_median_m": 1.0,          # warn when median DEM - ground exceeds this (datum, or points off)
        "warn_class_spread_m": 1.0,    # warn when low-canopy offsets differ this much between classes
    },
}


# --------------------------------------------------------------------------------------
# 1. Points
# --------------------------------------------------------------------------------------
def _empty():
    import geopandas as gpd
    return gpd.GeoDataFrame({"lon": [], "lat": [], "h_ell": [], "source": [], "time": []},
                            geometry=[], crs="EPSG:4326")


def _first_col(df, names):
    for n in names:
        if n in df.columns:
            return n
    return None


def _year_from_time(series) -> np.ndarray:
    try:
        t = pd.to_datetime(series, utc=True, errors="coerce")
        return (t.dt.year + (t.dt.dayofyear - 1) / 365.25).to_numpy(dtype="float64")
    except Exception:
        return np.full(len(series), np.nan)


def fetch_sliderule(aoi_ll, gcfg: dict):
    """Fetch from SlideRule, splitting a large AOI into tiles. SlideRule refuses requests that
    match more than 300 ATL03 granules, so a region larger than about 0.5 degrees fails whole."""
    import geopandas as gpd
    from shapely.geometry import box
    tile = float(gcfg.get("tile_deg") or 0)
    x0, y0, x1, y1 = aoi_ll.bounds
    if tile <= 0 or max(x1 - x0, y1 - y0) <= tile:
        return _fetch_one(aoi_ll, gcfg)
    frames, nfail = [], 0
    xs = np.arange(x0, x1, tile)
    ys = np.arange(y0, y1, tile)
    log.info("SlideRule: AOI split into up to %d tiles of %.2f deg", len(xs) * len(ys), tile)
    for x in xs:
        for y in ys:
            part = aoi_ll.intersection(box(x, y, min(x + tile, x1), min(y + tile, y1)))
            if part.is_empty or part.area < 1e-6:
                continue
            try:
                g = _fetch_one(part, gcfg)
            except Exception as e:
                nfail += 1
                log.error("SlideRule tile (%.2f, %.2f) failed: %s", x, y, e)
                continue
            if len(g):
                frames.append(g)
    if nfail:
        log.warning("SlideRule: %d tile(s) failed; ground points are incomplete", nfail)
    if not frames:
        return _empty()
    df = pd.concat(frames, ignore_index=True)
    return gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.lon, df.lat), crs="EPSG:4326")


def _fetch_one(aoi_ll, gcfg: dict):
    """Return GeoDataFrame [lon, lat, h_ell, source, time, year] from SlideRule."""
    import geopandas as gpd
    from sliderule import sliderule, icesat2, gedi

    sliderule.init(verbose=False)
    poly = aoi_ll.simplify(gcfg["simplify_aoi_deg"]) if gcfg["simplify_aoi_deg"] else aoi_ll
    region = sliderule.toregion(gpd.GeoDataFrame(geometry=[poly], crs="EPSG:4326"))
    frames = []

    if "icesat2_atl08" in gcfg["datasets"]:
        ic = gcfg["icesat2"]
        parms = {"poly": region["poly"], "t0": gcfg["t0"], "t1": gcfg["t1"],
                 "len": ic["len"], "res": ic["res"], "atl08_class": ic["atl08_class"],
                 "phoreal": ic["phoreal"]}
        log.info("SlideRule: requesting ICESat-2 ATL08 (PhoREAL) ...")
        g = icesat2.atl08p(parms)
        log.info("  ATL08 segments returned: %d", len(g))
        if len(g):
            hcol = _first_col(g, ["h_te_median", "h_te_best_fit", "h_te_mean"])
            if hcol is None:
                raise RuntimeError(f"no terrain height column in ATL08 result: {list(g.columns)}")
            keep = np.isfinite(g[hcol].to_numpy(dtype="float64"))
            if ic["night_only"] and "solar_elevation" in g:
                keep &= g["solar_elevation"].to_numpy() < 0
            pc = _first_col(g, ["ground_photon_count", "n_te_photons", "gnd_ph_count"])
            if pc and ic["min_ground_photons"]:
                keep &= g[pc].to_numpy() >= ic["min_ground_photons"]
            g = g[keep]
            tcol = g.index if "time" not in g.columns else g["time"]
            frames.append(pd.DataFrame({
                "lon": g.geometry.x.to_numpy(), "lat": g.geometry.y.to_numpy(),
                "h_ell": g[hcol].to_numpy(dtype="float64"), "source": "icesat2_atl08",
                "time": pd.Series(tcol).astype(str).to_numpy(),
                "track": (g["rgt"].astype(str) + "_" + g["cycle"].astype(str)).to_numpy()
                         if {"rgt", "cycle"} <= set(g.columns) else "",
            }))

    if "gedi_l2a" in gcfg["datasets"]:
        ge = gcfg["gedi"]
        parms = {"poly": region["poly"], "t0": gcfg["t0"], "t1": gcfg["t1"],
                 "degrade_flag": ge["degrade_flag"], "l2_quality_flag": ge["l2_quality_flag"]}
        log.info("SlideRule: requesting GEDI L2A ...")
        g = gedi.gedi02ap(parms)
        log.info("  GEDI footprints returned: %d", len(g))
        if len(g):
            hcol = _first_col(g, ["elevation_lm", "elev_lowestmode", "elevation_lowestmode"])
            if hcol is None:
                raise RuntimeError(f"no ground elevation column in GEDI result: {list(g.columns)}")
            keep = np.isfinite(g[hcol].to_numpy(dtype="float64"))
            sc = _first_col(g, ["sensitivity"])
            if sc and ge["min_sensitivity"]:
                keep &= g[sc].to_numpy() >= ge["min_sensitivity"]
            g = g[keep]
            tcol = g.index if "time" not in g.columns else g["time"]
            frames.append(pd.DataFrame({
                "lon": g.geometry.x.to_numpy(), "lat": g.geometry.y.to_numpy(),
                "h_ell": g[hcol].to_numpy(dtype="float64"), "source": "gedi_l2a",
                "time": pd.Series(tcol).astype(str).to_numpy(),
                "track": g["beam"].astype(str).to_numpy() if "beam" in g.columns else "",
            }))

    if not frames:
        return _empty()
    df = pd.concat(frames, ignore_index=True)
    return gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.lon, df.lat), crs="EPSG:4326")


def load_local(gcfg: dict):
    import geopandas as gpd
    lc = gcfg["local"]
    path = lc["path"]
    if not path:
        raise ValueError("ground.source = local but ground.local.path is empty")
    if str(path).lower().endswith((".csv", ".txt")):
        df = pd.read_csv(path)
        g = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df[lc["lon_col"]], df[lc["lat_col"]]),
                             crs="EPSG:4326")
    else:
        g = gpd.read_file(path).to_crs(4326)
    out = pd.DataFrame({
        "lon": g.geometry.x, "lat": g.geometry.y, "h_ell": g[lc["h_col"]].astype("float64"),
        "source": g[lc["source_col"]].astype(str) if lc["source_col"] in g else "local",
        "time": g[lc["time_col"]].astype(str) if lc["time_col"] in g else "",
        "track": "",
    })
    gdf = gpd.GeoDataFrame(out, geometry=g.geometry.values, crs="EPSG:4326")
    gdf.attrs["already_orthometric"] = not lc["height_is_ellipsoidal"]
    return gdf


def get_points(aoi_ll, gcfg: dict, out_dir):
    import geopandas as gpd
    cache = out_dir / "ground_points_raw.gpkg"
    if gcfg["source"] == "none":
        return _empty()
    if gcfg["cache"] and cache.exists() and gcfg["source"] == "sliderule":
        log.info("ground points: using cache %s", cache)
        return gpd.read_file(cache)
    if gcfg["source"] == "local":
        g = load_local(gcfg)
    elif gcfg["source"] == "sliderule":
        g = fetch_sliderule(aoi_ll, gcfg)
        if len(g) and gcfg["cache"]:
            g.to_file(cache, driver="GPKG")
    else:
        raise ValueError(f"unknown ground.source: {gcfg['source']}")
    log.info("ground points loaded: %d", len(g))
    return g


# --------------------------------------------------------------------------------------
# 2. Datum and epoch
# --------------------------------------------------------------------------------------
def to_dem_datum(g, gcfg: dict, lc_at_points: np.ndarray | None):
    """Add column H (height in DEM datum, at DEM epoch)."""
    d = gcfg["datum"]
    h = g["h_ell"].to_numpy(dtype="float64")
    if g.attrs.get("already_orthometric") or d["mode"] == "none":
        H = h.copy()
    elif d["mode"] == "constant":
        H = h - float(d["constant_m"])
    elif d["mode"] == "raster":
        import rasterio
        with rasterio.open(d["raster_path"]) as src:
            from rasterio.warp import transform as wtransform
            xs, ys = wtransform("EPSG:4326", src.crs, g["lon"].tolist(), g["lat"].tolist())
            N = np.array([v[0] for v in src.sample(zip(xs, ys))], dtype="float64")
        H = h - N
    elif d["mode"] == "pyproj":
        import pyproj
        pyproj.network.set_network_enabled(True)   # fetch geoid grid from the PROJ CDN if needed
        tr = pyproj.Transformer.from_crs("EPSG:4979", d["target"], always_xy=True)
        _, _, H = tr.transform(g["lon"].to_numpy(), g["lat"].to_numpy(), h)
        H = np.asarray(H, dtype="float64")
        if not np.isfinite(H).any() or np.allclose(H, h):
            raise RuntimeError(
                "pyproj could not apply the geoid model (grid unavailable?). "
                "Use ground.datum.mode: raster with a geoid GeoTIFF, or constant.")
    else:
        raise ValueError(f"unknown ground.datum.mode: {d['mode']}")

    ep = gcfg["epoch"]
    if ep.get("dem_year"):
        years = _year_from_time(g["time"]) if "time" in g else np.full(len(g), np.nan)
        dt = np.nan_to_num(years - float(ep["dem_year"]), nan=0.0)
        rate = np.full(len(g), float(ep.get("subsidence_default_m_per_yr", 0.0)))
        if lc_at_points is not None:
            for k, v in (ep.get("subsidence_m_per_yr_by_class") or {}).items():
                rate[lc_at_points == int(k)] = float(v)
        H = H + rate * dt   # ground was higher at DEM time if the peat has subsided since
    g = g.copy()
    g["H"] = H
    return g


# --------------------------------------------------------------------------------------
# 3-4. Sampling and split
# --------------------------------------------------------------------------------------
def to_grid_xy(g, grid):
    g2 = g.to_crs(grid.crs)
    return g2.geometry.x.to_numpy(), g2.geometry.y.to_numpy()


def sample_arrays(x, y, grid, arrays: dict) -> dict:
    inv = ~grid.transform
    cols, rows = inv * (x, y)
    r = np.floor(rows).astype(int)
    c = np.floor(cols).astype(int)
    inside = (r >= 0) & (r < grid.height) & (c >= 0) & (c < grid.width)
    out = {"inside": inside}
    for name, arr in arrays.items():
        v = np.full(len(x), np.nan)
        if arr is not None:
            v[inside] = arr[r[inside], c[inside]]
        out[name] = v
    return out


def edge_mask(lc: np.ndarray | None, buffer_m: float, res: float) -> np.ndarray | None:
    if lc is None or buffer_m <= 0:
        return None
    from scipy import ndimage
    l = np.nan_to_num(lc, nan=-1)
    change = (ndimage.maximum_filter(l, size=3) != ndimage.minimum_filter(l, size=3))
    it = max(1, int(round(buffer_m / res)))
    return ndimage.binary_dilation(change, iterations=it)


def spatial_split(x, y, cfg_split: dict) -> np.ndarray:
    """True = validation. Whole spatial blocks go to one side."""
    bs = float(cfg_split["block_size_m"])
    bx = np.floor(x / bs).astype(np.int64)
    by = np.floor(y / bs).astype(np.int64)
    keys = bx * 1_000_003 + by
    uniq = np.unique(keys)
    rng = np.random.default_rng(int(cfg_split["seed"]))
    nval = max(1, int(round(len(uniq) * float(cfg_split["validation_fraction"])))) if len(uniq) > 1 else 0
    val_blocks = set(rng.choice(uniq, nval, replace=False).tolist()) if nval else set()
    return np.array([k in val_blocks for k in keys])


# --------------------------------------------------------------------------------------
# 5. Fit
# --------------------------------------------------------------------------------------
def _huber_mean(v, c, iters=10):
    m = np.median(v)
    for _ in range(iters):
        r = v - m
        w = np.where(np.abs(r) <= c, 1.0, c / np.maximum(np.abs(r), 1e-9))
        m = np.sum(w * v) / np.sum(w)
    return float(m)


def _huber_slope_origin(xv, yv, c, iters=15):
    """Robust slope k for y = k * x (through origin)."""
    denom = np.sum(xv * xv)
    if denom <= 0:
        return float("nan"), float("nan")
    k = np.sum(xv * yv) / denom
    for _ in range(iters):
        r = yv - k * xv
        w = np.where(np.abs(r) <= c, 1.0, c / np.maximum(np.abs(r), 1e-9))
        k = np.sum(w * xv * yv) / max(np.sum(w * xv * xv), 1e-12)
    resid = yv - k * xv
    se = float(np.sqrt(np.sum(resid ** 2) / max(len(xv) - 1, 1) / denom))
    return float(k), se


def _huber_line(xv, yv, c, iters=15):
    """Robust y = a + k * x. Returns (a, k, se_k)."""
    A = np.c_[np.ones(len(xv)), xv]
    w = np.ones(len(xv))
    for _ in range(iters):
        sw = np.sqrt(w)
        coef, *_ = np.linalg.lstsq(A * sw[:, None], yv * sw, rcond=None)
        r = yv - A @ coef
        w = np.where(np.abs(r) <= c, 1.0, c / np.maximum(np.abs(r), 1e-9))
    cov = np.linalg.pinv((A * w[:, None]).T @ A) * (np.sum(w * r ** 2) / max(len(xv) - 2, 1))
    return float(coef[0]), float(coef[1]), float(np.sqrt(abs(cov[1, 1])))


def datum_diagnostics(df: pd.DataFrame, dcfg: dict, chm_low_m: float) -> dict:
    """Median DEM - ground per class for low-canopy points. A constant offset shared by all
    classes points to the vertical datum; offsets that differ by class point to the ground
    points, the land cover or the DEM itself."""
    use = df[df["exclude_reason"] == ""]
    r = (use["dem"] - use["H"]).to_numpy()
    low = (use["chm"].to_numpy() < chm_low_m)
    out = {"n_used": int(len(use)), "median_all_m": float(np.median(r)) if len(r) else None,
           "geoid_n_median_m": float(np.median(df["h_ell"] - df["H"])) if "h_ell" in df else None,
           "low_canopy_by_class": {}}
    meds = []
    for c in sorted({int(v) for v in use["lc"].dropna()}):
        m = low & (use["lc"].to_numpy() == c)
        if m.sum() >= 10:
            med = float(np.median(r[m]))
            out["low_canopy_by_class"][str(c)] = {
                "n": int(m.sum()), "median_m": med,
                "nmad_m": float(1.4826 * np.median(np.abs(r[m] - med)))}
            meds.append(med)
    out["class_spread_m"] = float(max(meds) - min(meds)) if len(meds) > 1 else None
    if out["median_all_m"] is not None and abs(out["median_all_m"]) > float(dcfg["warn_median_m"]):
        log.warning("median DEM - ground = %.2f m: vertical datum or ground points need checking",
                    out["median_all_m"])
    if out["class_spread_m"] is not None and out["class_spread_m"] > float(dcfg["warn_class_spread_m"]):
        log.warning("low-canopy offset differs by %.2f m between classes: a single global bias "
                    "will not fit; consider fit.joint_intercept", out["class_spread_m"])
    return out


def fit_k_bias(pts: pd.DataFrame, prior_k: dict, prior_k_default: float, fcfg: dict) -> dict:
    """pts needs columns dem, H, chm, lc. Returns calibration result."""
    r = pts["dem"].to_numpy() - pts["H"].to_numpy()
    chm = np.nan_to_num(pts["chm"].to_numpy(), nan=0.0)
    lc = pts["lc"].to_numpy()
    res = {"n_points": int(len(pts)), "classes": {}}

    low = chm < float(fcfg["bias_from_chm_below_m"])
    if low.sum() >= int(fcfg["min_points_per_class"]):
        b0 = _huber_mean(r[low], float(fcfg["huber_m"]))
        res["bias_global"] = {"value": b0, "n": int(low.sum()), "from": "low-canopy points"}
    else:
        b0 = 0.0
        res["bias_global"] = {"value": 0.0, "n": int(low.sum()),
                              "from": "default (too few low-canopy points)"}

    kmin, kmax = fcfg["k_bounds"]
    classes = sorted({int(c) for c in lc[np.isfinite(lc)]})
    for cls in classes:
        sel = lc == cls
        n = int(sel.sum())
        prior = float(prior_k.get(str(cls), prior_k.get(cls, prior_k_default)))
        entry = {"n": n, "prior_k": prior}
        veg = sel & (chm >= float(fcfg["bias_from_chm_below_m"]))
        if (fcfg.get("joint_intercept") and veg.sum() >= int(fcfg["min_points_per_class"])
                and np.ptp(chm[veg]) > 1.0):
            a, k, se = _huber_line(chm[veg], r[veg], float(fcfg["huber_m"]))
            entry.update({"k": float(np.clip(k, kmin, kmax)), "k_raw": k, "k_se": se,
                          "n_vegetated": int(veg.sum()), "bias": a,
                          "status": "fitted (joint intercept)"})
            if not (kmin <= k <= kmax):
                entry["status"] = "fitted (joint intercept), clipped to bounds"
        elif veg.sum() >= int(fcfg["min_points_per_class"]):
            bias_c = b0
            if fcfg["fit_bias_per_class"] and (sel & low).sum() >= int(fcfg["min_points_per_class"]):
                bias_c = _huber_mean(r[sel & low], float(fcfg["huber_m"]))
            k, se = _huber_slope_origin(chm[veg], r[veg] - bias_c, float(fcfg["huber_m"]))
            entry.update({"k": float(np.clip(k, kmin, kmax)), "k_raw": k, "k_se": se,
                          "n_vegetated": int(veg.sum()), "bias": bias_c, "status": "fitted"})
            if not (kmin <= k <= kmax):
                entry["status"] = "fitted, clipped to bounds"
        else:
            entry.update({"k": prior, "bias": b0, "n_vegetated": int(veg.sum()),
                          "status": "prior (too few vegetated points)"})
        res["classes"][str(cls)] = entry
        log.info("  class %3d: n=%5d  k=%.2f (%s)", cls, n, entry["k"], entry["status"])
    log.info("  global bias = %.2f m (%s)", b0, res["bias_global"]["from"])
    return res


# --------------------------------------------------------------------------------------
# 6. Validation
# --------------------------------------------------------------------------------------
def error_table(H, surfaces: dict, lc, names: dict | None = None) -> pd.DataFrame:
    rows = []
    groups = [("all", np.ones(len(H), bool))]
    for c in sorted({int(v) for v in lc[np.isfinite(lc)]}):
        groups.append((str(c), lc == c))
    for sname, z in surfaces.items():
        e = z - H
        for gname, sel in groups:
            ee = e[sel & np.isfinite(e)]
            if ee.size == 0:
                continue
            med = np.median(ee)
            rows.append({"surface": sname, "class": gname,
                         "class_name": (names or {}).get(gname, gname), "n": int(ee.size),
                         "bias_m": float(ee.mean()), "rmse_m": float(np.sqrt(np.mean(ee ** 2))),
                         "mae_m": float(np.mean(np.abs(ee))),
                         "nmad_m": float(1.4826 * np.median(np.abs(ee - med)))})
    return pd.DataFrame(rows)


WORLDCOVER_NAMES = {"10": "tree cover", "20": "shrubland", "30": "grassland", "40": "cropland",
                    "50": "built-up", "60": "bare/sparse", "70": "snow/ice", "80": "water",
                    "90": "herbaceous wetland", "95": "mangroves", "100": "moss/lichen",
                    "all": "all classes"}
