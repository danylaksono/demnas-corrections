# demnas_dtm: DEMNAS to DTM correction, calibration and smoothing

Turns DEMNAS (a surface model: canopy tops in forest) into a ground estimate, calibrates the correction against ICESat-2 (and optionally GEDI) ground heights, and smooths the result to a dome-scale surface for canal block siting. One command runs everything.

```
DEMNAS tiles ─► mosaic, clip to KHG, reproject to UTM
                     │
land cover + canopy height (cloud COGs by default, or local files)
                     │
ground points (ICESat-2 / GEDI via SlideRule, or local) ─► fit k and bias per class
                     │
DTM = DEM − k[land cover] × canopy height − bias[land cover]
                     │
mask water, built-up, optional vectors (canals, settlements) ─► remove spikes
                     │
smoothing steps (median, gaussian, polynomial trend) ─► fill gaps ─► GeoTIFFs + report
                     │
validation on held-out ground points ─► validation.csv
```

If no ground points are available (no network, no coverage), the run continues with the prior `k` values in the config and says so in the log and report.

## Install

```bash
pip install -r requirements.txt
```

## Run

```bash
# defaults: download corrections, median 300 m then gaussian 150 m
python demnas_dtm.py --dem "demnas/*.tif" --aoi khg_01.gpkg --out out_khg01

# full control through a config file (CLI flags override it)
python demnas_dtm.py --print-default-config > my_config.yaml
python demnas_dtm.py --config my_config.yaml

# offline, with your own correction rasters
python demnas_dtm.py --dem dem.tif --out out --no-download --landcover lc.tif --chm chm.tif

# add GEDI for denser ground points
python demnas_dtm.py --dem "demnas/*.tif" --aoi khg_01.gpkg --out out --datasets icesat2_atl08 gedi_l2a

# your own ground points (CSV with lon, lat, h)
python demnas_dtm.py --dem dem.tif --aoi khg.gpkg --out out --ground-points points.csv

# skip calibration
python demnas_dtm.py --dem dem.tif --out out --ground none
```

## Correction sources (mode: `download`)

| Layer | Default source | Read as |
|---|---|---|
| Land cover | ESA WorldCover 10 m v200 (2021), public COGs on AWS | nearest neighbour onto the DEM grid |
| Canopy height | Meta/WRI global canopy height, public COGs on AWS | average onto the DEM grid, using a COG overview near the target resolution so only small reads are needed |

Only the windows covering the AOI are fetched. Both URL templates are config keys (`corrections.worldcover_url`, `corrections.chm_url`). If a provider moves its files, update the template rather than the code. Any tile that can't be read is logged and skipped:
- Where canopy height is missing, per-class fallback heights apply (`chm_fallback_m`).
- Where land cover is missing, `k_default` applies everywhere.

Land cover classes follow ESA WorldCover codes (10 tree cover, 20 shrubland, 80 water, 90 herbaceous wetland, 95 mangroves, …). A local land cover raster must use the same codes, or `k_by_class` must be rewritten for your classes.

## Ground truthing

1. **Points.** By default the tool asks SlideRule's public service for ICESat-2 ATL08 segments (PhoREAL processing, 100 m segments) inside the AOI. Add `gedi_l2a` for GEDI footprints. Results are cached in `ground_points_raw.gpkg`, so reruns don't fetch again. With `--ground-points` you supply your own CSV (`lon`, `lat`, `h`; optional `time`, `source`) or vector file instead, for example RTK survey points.
2. **Datum.** ICESat-2 and GEDI heights are ellipsoidal. `ground.datum.mode` converts them to the DEM's datum:
   - `pyproj` (default): converts to WGS84 + EGM2008 heights using PROJ, fetching the geoid grid from the PROJ CDN if needed.
   - `raster`: uses your own geoid undulation GeoTIFF, e.g. INAGeoid2020.
   - `constant`: subtracts one value. Only for quick tests.
   - `none`: for points already in the DEM datum.
3. **Epoch (optional).** Set `ground.epoch.dem_year` and subsidence rates per class to shift ground heights back to the DEM acquisition year. Peat that subsided since then was higher when DEMNAS was made.
4. **Filtering.** Points are dropped when they are:
   - in water or built-up land cover;
   - within 50 m of a land-cover change, where footprints mix surfaces;
   - inside the vector masks, such as canal buffers;
   - more than 30 m from the DEM.

   Each dropped point keeps its `exclude_reason`.
5. **Split.** Usable points go to calibration or validation by **2 km spatial blocks**, 30% validation. Neighbouring footprints along a track are correlated, so a random split would overstate accuracy.
6. **Fit.** Model: `DEM − ground = bias + k[class] × canopy height`.
   - **bias:** a robust (Huber) mean over points with canopy below 2 m.
   - **k:** a robust slope through the origin for each class with at least 30 vegetated points, clipped to 0–1. Classes with too few points keep their prior k.
7. **Validation.** Errors at validation points (bias, RMSE, MAE, NMAD), per class and overall, for four surfaces: raw DEM, prior-k DTM, calibrated DTM and smoothed DTM.

## Key parameters

| Key | Default | Meaning |
|---|---|---|
| `target_res_m` | 8 | Output resolution (DEMNAS native ≈ 8 m) |
| `target_crs` | auto-utm | UTM zone from the AOI centre (south or north) |
| `aoi_buffer_m` | 500 | Extra margin so smoothing behaves at the KHG edge |
| `k_by_class` | 0.8 forest/mangrove, 0.5 shrub, 0.4 wetland, 0.2 grass/crop | Share of canopy height present in the DEM (**placeholder**) |
| `bias_by_class` | none | Constant offset per class (m), for later calibration |
| `corrections.chm_min_m` / `chm_max_m` | 2 / 60 | Canopy below 2 m treated as 0; above 60 m clipped |
| `mask.landcover_classes` | [80, 50] | Water and built-up excluded from the smoothing input |
| `mask.vectors` | none | e.g. `[{path: canals.gpkg, buffer_m: 15}]` so canal cuts don't pull the surface down |
| `mask.despike_threshold_m` | 2.5 | Cells further than this from a 400 m local median are removed |
| `smoothing.steps` | median 300 m → gaussian 150 m | Applied in order (see below) |
| `ground.source` | sliderule | `sliderule`, `local` or `none` |
| `ground.datasets` | icesat2_atl08 | add `gedi_l2a` for more points |
| `ground.t0` / `t1` | 2019–2025 | time window for ICESat-2 / GEDI |
| `ground.datum.mode` | pyproj (EGM2008) | must match the DEM's vertical datum |
| `ground.split.block_size_m` | 2000 | spatial block size for the calibration/validation split |
| `ground.fit.min_points_per_class` | 30 | below this, the prior k is kept |
| `ground.fit.joint_intercept` | false | true = fit `DEM − ground = a + k × canopy` per class instead of one global bias. Use when `datum_diagnostics.json` shows offsets that differ by class |
| `ground.diagnostics.warn_median_m` | 1.0 | warn when the median DEM − ground exceeds this (datum or ground point problem) |
| `ground.fail_on_error` | false | true = stop if ground truthing fails |

Smoothing methods, combinable in any order:
- `median` (`window_m`): robust large-window median, fast on large grids.
- `gaussian` (`sigma_m`): ignores masked cells instead of averaging them in as zeros.
- `polynomial` (`degree`, `huber_m`, `sample`): robust trend surface for the overall dome shape. Use degree 3–4 for one dome; follow it with a gentle gaussian if you want some local form back.

## Outputs

| File | Content |
|---|---|
| `dtm_smooth.tif` | Final smoothed ground surface: input for contours, profiles and flow direction |
| `dtm_residual.tif` | Corrected DTM minus smoothed surface: local relief and remaining noise |
| `dtm_corrected.tif` | Vegetation-corrected, unsmoothed DTM |
| `dem_utm.tif` | DEMNAS mosaic on the output grid |
| `correction_m.tif`, `k_factor.tif` | What was subtracted, and the factor used |
| `landcover.tif`, `canopy_height.tif` | Corrections resampled to the DEM grid |
| `excluded_mask.tif` | Cells kept out of smoothing (1 = excluded) |
| `ground_points_raw.gpkg` | Ground points as fetched (cache) |
| `ground_points.gpkg` | All points with datum-corrected height, samples, `exclude_reason`, `split` and errors for each surface |
| `calibration.json` | Fitted k (with standard error and point counts) and bias per class |
| `datum_diagnostics.json` | Median geoid undulation, median DEM − ground, and low-canopy offset per class |
| `validation.csv` | Error statistics on validation points per surface and class |
| `run_report.json` | Config snapshot, sources actually used, k values applied, calibration and validation summary, runtime |

## Things to keep in mind

- **Vertical datum.** Outputs stay in DEMNAS's datum. Ground points are converted to it, so `ground.datum.target` must match the DEMNAS metadata (commonly EGM2008). A median DEM-minus-ground difference above 20 m triggers a warning: almost always a datum mismatch.
- **SlideRule fields.** The ATL08 and GEDI column names (`h_te_median`, `elevation_lm`, `sensitivity`) follow the SlideRule client at the time of writing. If the service changes them, the run stops ground truthing with an error naming the columns it received. You can then export points yourself and use `--ground-points`.
- **Dates differ.** DEMNAS, WorldCover (2021) and the canopy height product come from different years. Clearing or regrowth since DEMNAS was acquired will be over- or under-corrected; the calibration step will absorb part of this.
- **Memory.** A 50 × 50 km KHG at 8 m is about 40 million cells; expect roughly 1–2 GB of RAM. Use `--res 10` or 15 for quick runs.
- **Synthetic test.** On a synthetic dome with 25 m forest canopy, the run started from deliberately wrong priors (k = 0.5 forest, 0.3 shrub). Calibration on simulated ICESat-2 tracks recovered k = 0.80 and 0.49, close to the true 0.8 and 0.5. On held-out points, RMSE went from 12.4 m (raw DEM) and 4.8 m (prior k) to 1.1 m (calibrated) and 0.4 m (smoothed). The SlideRule and PROJ downloads could not be tested from the build environment.
