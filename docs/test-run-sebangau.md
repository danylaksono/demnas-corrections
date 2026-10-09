# Test run on the Sebangau sample: log of runs and decisions

Date: 2026-10-08. Input: `_data_sample/DEMNAS_sebangau.tif` (gitignored).
Purpose: check that `demnas_dtm.py` runs end to end on real DEMNAS data, including ground truthing.

## Input

| Property | Value |
|---|---|
| CRS | EPSG:4326 |
| Size | 17,341 x 9,323 px (about 7.5e-5 deg, roughly 8 m) |
| Extent | lon 113.581 to 114.280, lat -3.469 to -2.169 |
| Nodata | -9999 |
| Valid data | A diagonal strip, not the whole rectangle. Heights up to about 25 m on a coarse read. |

## Environment

- Python 3.14.7 in a throwaway venv (not in the repo). `pip install -r requirements.txt` worked unchanged, with no build problems on 3.14.
- No code in the repo was changed during these runs. The only repo changes are `.gitignore` (ignores `runs/`, `_data_sample/`, `__pycache__/`) and this report.
- Run outputs live in `runs/` (gitignored). Seen from `runs/`:

| Folder | What |
|---|---|
| `sebangau_30m/` | Run 1, full tile, no ground truthing |
| `aoi_window.gpkg` | AOI for run 3 |
| `window_calib/` | Run 3, windowed, with ground truthing |
| `window_calib.log` | Log of run 3 |

## Runs

### Run 1: full tile, 30 m, ground truthing off

Command: `python demnas_dtm.py --dem _data_sample/DEMNAS_sebangau.tif --out runs/sebangau_30m --res 30 --ground none`

**Decision: resolution 30 m and no ground truthing.** The tile is about 160 million cells at native resolution, and a SlideRule request over roughly 78 x 144 km would be very large. The aim was to prove the DEM, corrections and smoothing path first.

Result: completed in 1,123 s. Grid 2607 x 4806 px in EPSG:32749 (UTM 49S). Most of the time went to reading canopy height tiles (about 4-5 min each, no progress logged between tiles). 13.0% of cells excluded in total, 5.4% of them as spikes.

| Surface | Median (m) | Notes |
|---|---|---|
| DEM | 1.9 | |
| Smoothed DTM | 0.6 | min -12.0, max 21.5 |
| Correction | 0 | 95th percentile 13.6 |

k values were the uncalibrated placeholders (forest 0.8 and so on), so this run only shows the pipeline works. Residuals reach -19 m to +25 m. I did not investigate them. Likely cause is mask edges around water and built-up land.

An early attempt piped the output through `tail`, so the log showed nothing until the process ended. It looked hung and was stopped, then rerun with output to a file. It had not been hung.

### Run 2: windowed run, failed (network)

**Decision: use a 15 x 15 km window instead of the full tile**, so the ground-point fetch is a realistic size. The valid data is a diagonal strip, so I picked the centre from a coarse look at the data, 113.905 E, 3.094 S, half-width 0.0675 deg, saved as `runs/aoi_window.gpkg`. A check showed 100% valid DEM cells in the window (heights -2.4 to 10.6 m).

The first attempt at 8 m lost network access (DNS failures for WorldCover, canopy height and SlideRule). The pipeline kept going by design: no land cover, no canopy height, no ground points, `k_default` everywhere. The outputs were therefore meaningless, and I deleted them.

Note for later: with `ground.fail_on_error: false` (the default), a network failure still ends in a normal "done" with a warning-only trail. Check the log for "No land cover available" and "no ground points available".

### Run 3: windowed run with ground truthing

Command: `python demnas_dtm.py --dem _data_sample/DEMNAS_sebangau.tif --aoi runs/aoi_window.gpkg --out runs/window_calib --res 8`
All other settings are the defaults: ESA WorldCover and Meta/WRI canopy height from the cloud, ICESat-2 ATL08 through SlideRule (2019-2025), datum mode `pyproj` (EGM2008), median 300 m then gaussian 150 m smoothing. Completed in 290 s.

**Ground points.** ATL08 returned 7,810 segments, 6,689 were loaded and 1,514 were usable after filtering. Of these 1,073 went to calibration and 441 to validation, split in 2 km spatial blocks.

**Calibration** (`calibration.json`):

| Class | n | k | Status |
|---|---|---|---|
| 10 tree cover | 434 | 0.45 (se 0.019) | fitted, 186 vegetated points |
| 20 shrubland | 236 | 0.50 | prior kept, 18 vegetated points |
| 30 grassland | 403 | 0.20 | prior kept, 3 vegetated points |

Global bias is -3.77 m, from low-canopy points.

**Validation** (441 held-out points):

| Surface | Bias (m) | RMSE (m) | NMAD (m) |
|---|---|---|---|
| Raw DEM | -1.92 | 3.19 | 3.26 |
| Prior k | -3.49 | 4.25 | 1.58 |
| Calibrated | +0.95 | 2.31 | 1.73 |
| Smoothed | +0.95 | 1.97 | 1.69 |

By class (bias in m):

| Class | Raw DEM | Calibrated | Smoothed |
|---|---|---|---|
| 10 tree cover | +0.04 | +2.15 | +2.09 |
| 20 shrubland | -4.13 | -0.39 | -0.34 |
| 30 grassland | -4.13 | -0.40 | -0.36 |

## Findings from run 3 (first reading)

1. **The pipeline works end to end** on the sample: mosaic, reprojection, WorldCover and canopy height fetch, SlideRule fetch, datum conversion, calibration, validation and outputs. The README says the SlideRule and PROJ paths had never been tested against the live services. This run exercised both without errors.
2. **A large DEM - ground offset in the open classes.** Where canopy is near zero (shrubland, grassland) the DEM sits about 4.1 m below the ICESat-2 ground, with a tight spread (NMAD about 0.7 m).
3. **The calibration absorbed that offset into one global bias** (-3.77 m, so +3.77 m on the DTM) and applied it to every class. In forest the raw DEM was already close to ground (+0.04 m), so the calibrated DTM ended up about 2.1 m too high there.
4. **Two of the three classes had too few vegetated points** (18 and 3) to fit k, so they kept the priors.
5. `dtm_corrected.tif` is the unsmoothed vegetation-corrected surface. `dtm_smooth.tif` is the smoothed surface for dome-scale work. Neither should be used as a final product until the offset question below is settled.

## Diagnostics on the cached ICESat-2 points (offline)

These checks reuse `runs/window_calib/ground_points.gpkg` and need no network. They revise the reading above.

| Check | Result | What it rules in or out |
|---|---|---|
| Geoid undulation applied by `pyproj` | N = 43.0 to 43.7 m (median 43.3 m) | The geoid grid is applied. Skipping it (`mode: none`) would be 43 m off, so the 4 m offset is not a missing geoid conversion. |
| Offset by acquisition year (2019-2025, open points) | -4.07, -3.99, -3.51, -3.75, -3.99, -4.47, -3.58 m, no trend | Not subsidence or a drift in time. |
| Offset by DEM height quartile (open points) | -4.68, -4.32, -3.67, -1.58 m from lowest to highest | Not a constant datum shift. The offset correlates 0.86 with DEM height. |
| Open-class ground heights | median ground 5.9 m; in forest 3.9 m. DEM medians are 1.9 m and 4.2 m. | ICESat-2 says the open land is 2 m higher than the forest, the DEM says the reverse. Superseded: see the GEDI section, which finds the same heights from a second sensor. |
| Open points: ground vs DEM | corr -0.32 (ground = -0.25 DEM + 6.2) | ICESat-2 ground in open classes does not track the DEM. |
| Canopy height in "tree cover" points | median 0 m, 75th percentile 6.4 m, 95th 10.6 m; 53% below 2 m | This window is degraded or burnt peat forest, not tall forest. |
| Forest points with canopy >= 2 m (n = 315) | DEM - ground = -0.15 + 0.045 x CHM | DEMNAS is close to ground here (k about 0.05). |
| Forest points with canopy < 2 m | median DEM - ground = -2.5 m | The same class gives -2.5 m at low canopy and -0.15 m at higher canopy. |
| Land-cover edge filter | 5,166 of 6,689 points (77%) dropped | The 50 m edge buffer removes most of the data. |

**Revised reading.** The offset is not a clean constant datum shift. It differs by class (-4 m in open classes, -2.5 m in forest at low canopy, -0.15 m in forest with canopy), it varies with DEM height, and it has no temporal trend. The calibration model assumes one constant bias plus k x canopy height, and that assumption does not hold here. The cause is not resolved. Candidates: a datum difference shared by both lidars, WorldCover class errors, spatial variation in the DEMNAS source data, or a real DEM error. The GEDI section below rules out an ICESat-2-only problem.

## Code changes made after run 3

Both changes are in `ground_truth.py` and `demnas_dtm.py`. Defaults are unchanged, so existing runs reproduce exactly (the `base` experiment below matches run 3 to the second decimal).

- **Joint per-class fit** (`ground.fit.joint_intercept: true`). Fits `DEM - ground = a[class] + k[class] x CHM` with a robust (Huber) line on all vegetated points of a class. Classes without enough vegetated points keep the prior k and the global bias.
- **Datum diagnostics** (`ground.diagnostics`). Every run now writes `datum_diagnostics.json`, and `run_report.json` carries it. It holds the median geoid undulation, the median DEM - ground, and the low-canopy offset per class. It warns when the median offset exceeds 1 m (the old threshold was 20 m, which missed this case) and when low-canopy offsets differ by more than 1 m between classes.
- **Exclusion summary.** The log now lists how many points each exclusion reason removed.
- **SlideRule tiling** (`ground.tile_deg`, default 0.5). Large AOIs are split into tiles, because SlideRule rejects a request that matches more than 300 ATL03 granules. See the Bengkalis notes below.

## Experiments

All experiments reuse the cached ICESat-2 points and the land cover and canopy rasters from run 3 (`corrections.mode: local`), at 8 m on the same window. Configs and logs are in `runs/exp_*.yaml` and `runs/exp_*.log`.

| Run | Edge buffer | Fit | Validation n | Smoothed bias (m) | Smoothed RMSE (m) |
|---|---|---|---|---|---|
| base | 50 m | global bias | 441 | +0.95 | 1.97 |
| joint | 50 m | joint intercept | 441 | -0.09 | 1.24 |
| edge10 | 10 m | joint intercept | 1006 | -0.16 | 1.47 |
| edge0 | none | joint intercept | 1752 | -1.14 | 2.12 |

**Comparability.** `base` and `joint` use the same 441 validation points, so their numbers compare directly. The `edge` runs have different validation sets, so their RMSE cannot be compared with those two.

- **Joint fit:** tree cover k falls from 0.45 to 0.04 with intercept -0.25 m, and the tree cover validation bias goes from +2.1 to +0.2 m. Shrubland and grassland are unchanged (-0.4 m) because they still use the global bias and prior k.
- **k stays near zero in every joint run:** 0.04 to 0.06 for tree cover, and 0 to 0.03 for shrubland and grassland once enough points are available. In this window the canopy height product does not explain the DEM - ground difference. The prior k of 0.8 is far from what the data show.
- **Smaller edge buffers keep more points** but add noisier points near class boundaries, and the offset between classes persists. The edge-free run has a -1.1 m bias.
- **Still unresolved:** the median DEM - ground stays near -3.4 to -3.9 m in all runs, and the warning fires in all of them.

## Second ground source: GEDI L2A

GEDI L2A (lowest-mode elevation, quality flag 1, sensitivity >= 0.95) was fetched through SlideRule for the same window together with ICESat-2 (`runs/exp_gedi.yaml`, joint-intercept fit, 50 m edge buffer, same canopy and land cover rasters). It returned 10,104 footprints, against 6,689 ICESat-2 segments. After filtering, 2,560 GEDI and 1,514 ICESat-2 points were usable, 4,074 in total (24% of 16,793).

**The two sensors agree.** Median DEM - ground by class and canopy, usable points only:

| Source | Class | Canopy | n | Ground H (m) | DEM - ground (m) | NMAD (m) |
|---|---|---|---|---|---|---|
| GEDI | shrubland | < 2 m | 626 | 6.12 | -4.31 | 1.15 |
| ICESat-2 | shrubland | < 2 m | 329 | 6.18 | -4.30 | 0.82 |
| GEDI | grassland | < 2 m | 387 | 5.59 | -4.06 | 0.86 |
| ICESat-2 | grassland | < 2 m | 490 | 5.85 | -4.26 | 0.83 |
| GEDI | tree cover | >= 2 m | 924 | 4.22 | +0.09 | 1.75 |
| ICESat-2 | tree cover | >= 2 m | 315 | 3.78 | +0.49 | 1.02 |
| GEDI | tree cover | < 2 m | 555 | 5.34 | -1.47 | 2.62 |
| ICESat-2 | tree cover | < 2 m | 353 | 5.79 | -2.54 | 1.95 |

In the 247 cells of 100 m that contain both sources, the median GEDI - ICESat-2 ground height is +0.25 m (NMAD 0.86 m). The two lidars see the same ground height to within about a quarter of a metre, and both put it about 4.1 to 4.3 m above DEMNAS in open land.

**What this changes.**
- The earlier suspicion that the open-class ICESat-2 heights were wrong (see the diagnostics table) is now much less likely. Two independent instruments, with different footprints, wavelengths and processing, give the same heights.
- It does not rule out a shared error. Both sources were converted to orthometric height with the same geoid (EGM2008 through PROJ), so a datum difference between DEMNAS and EGM2008 would appear identically in both. The independent LiDAR and GNSS data are what can separate a datum error from a DEM error.
- The pattern that open land has higher lidar ground than forest (about 5.3 to 6.2 m against 3.8 to 4.2 m) also holds in both sensors. The DEM shows the reverse. This is a second thing to explain: either DEMNAS is lower than the ground in open areas, or lidar ground in forest is lower than the DEM's.

**Calibration with both sources** (joint intercept): tree cover k = 0.06, shrubland k = 0.20 (now fitted, 58 and more vegetated points), grassland keeps the prior. Validation on 1,373 held-out points: smoothed bias -0.10 m, RMSE 2.22 m, NMAD 1.32 m (raw DEM: bias -1.91 m, RMSE 3.57 m). The RMSE is higher than the ICESat-2-only joint run (1.24 m) because the validation set is different and larger, with many more tree cover points (769), where the error is largest (smoothed RMSE 2.5 m). Do not compare the two numbers directly.

The median DEM - ground is still -2.9 m, and the datum warning still fires.

## Additional sites: Mendahara and Bengkalis

Two more DEMNAS samples were processed on 2026-10-08 with the same pipeline: `DEMNAS_mendahara.tif` (Jambi, 8,554 x 8,142 px) and `DEMNAS_bengkalis.tif` (Riau, 16,668 x 19,005 px, 1.3 GB). Both are EPSG:4326 at about 8 m with nodata -9999, and about 42-44% of the cells are valid (the rest is sea).

**Method choices.**
- **Memory.** Free RAM was about 1.9 GB, and the pipeline holds several copies of the mosaic. Each tile was first averaged 4x to about 30 m (`runs/prep/DEMNAS_<site>_x4.tif`, nodata-aware average) and then run at 30 m. Sebangau's 30 m run used the native 8 m tile, so the three sites differ slightly in how the DEM was reduced.
- **AOI.** A land-only polygon was built from the valid DEM cells (`runs/aoi_<site>.gpkg`) so the lidar requests skip open sea.
- **Settings.** `ground.fit.joint_intercept: true`, ICESat-2 ATL08 plus GEDI requested, all other settings default. Configs: `runs/mendahara.yaml`, `runs/bengkalis.yaml`, `runs/bengkalis_gt.yaml`.

**What went wrong first, and the fix.**
- **Mendahara GEDI.** The GEDI request hung for about 90 minutes and then failed with a connection error. The pipeline continued with ICESat-2 only. Mendahara therefore has no GEDI.
- **Bengkalis first run.** ICESat-2 returned 0 segments and GEDI returned HTTP 502. The pipeline logged this at INFO level and finished with the unfitted prior k, so `runs/bengkalis_30m/` is uncalibrated. Do not use its smoothed DEM.
- **Bengkalis cause.** On retry SlideRule reported "number of CMR hits <301> exceeded maximum allowed <300> for ATL03": a request over a large AOI matches too many granules and is rejected whole.
- **Fix.** `ground.tile_deg` (default 0.5 degrees) now splits large AOIs into tiles, requests each separately, logs and skips a failed tile, and warns that the result is incomplete. The Bengkalis AOI split into 8 non-empty tiles (of up to 9), none failed. This reran from the land cover and canopy rasters of the first run, so the 30 minute canopy download was not repeated. The calibrated Bengkalis result is `runs/bengkalis_30m_gt/`.

**Results.** Validation points are held out in 2 km spatial blocks (30%). Smoothed DEMs are `dtm_smooth.tif` in each output folder.

| Site | Ground points (usable / loaded) | Sources | Validation n | Raw DEM bias / RMSE (m) | Calibrated bias / RMSE (m) | Smoothed bias / RMSE / NMAD (m) |
|---|---|---|---|---|---|---|
| Sebangau window (15 km) | 1,514 / 6,689 (ICESat-2) | ICESat-2 | 441 | -1.92 / 3.19 | -0.10 / 1.49 | -0.09 / 1.24 / 1.12 |
| Mendahara (full tile) | 20,403 / 68,979 | ICESat-2 | 5,869 | +0.50 / 2.64 | +0.70 / 2.74 | +0.64 / 2.32 / 1.40 |
| Bengkalis (full tile) | 47,530 / 313,544 | ICESat-2 | 13,803 | +2.57 / 4.14 | +0.72 / 3.24 | +0.64 / 2.78 / 1.85 |

The Sebangau row is the joint-intercept ICESat-2 run on the 15 km window. The sites differ in size and in the number of points, so the rows are not a ranking of the sites.

**Datum diagnostics** (median DEM - ground, usable points; geoid undulation N):

| Site | N (m) | Median all (m) | Tree cover low-canopy (m) | Other classes, low canopy (m) |
|---|---|---|---|---|
| Sebangau | 43.3 | -3.4 | -2.5 | shrubland -4.3, grassland -4.3 |
| Mendahara | 10.0 | -0.03 | +0.09 | shrubland -1.88, grassland -0.77, cropland -0.08 |
| Bengkalis | 3.3 | +2.16 | +2.24 | shrubland -0.19, grassland +0.09, cropland -0.34, mangrove +1.06 |

**Findings.**
1. **The offset is site-specific.** The same conversion (PROJ, EGM2008) gives a DEM that is 3-4 m below ground at Sebangau, about right at Mendahara, and 2 m above ground at Bengkalis (where tree cover is dense and the DEM likely sees canopy). A single global datum error would give the same offset everywhere. This points to the DEM itself, or to local conditions, as a cause, not to the geoid conversion alone. It supports the idea that DEMNAS is a mosaic whose behaviour depends on the source data. It does not prove it.
2. **Fitted k is far below the prior.** Tree cover k is 0.10 (Mendahara), 0.21 (Bengkalis) and about 0.05 (Sebangau). The prior of 0.8 is not supported by any site. Using the prior would have made Mendahara much worse (RMSE 5.9 m against 2.6 m raw) and Bengkalis worse too (5.4 m against 4.1 m).
3. **Calibration helps where the DEM sits above ground** (Bengkalis RMSE 4.14 m to 3.24 m, bias +2.57 m to +0.72 m) and does not help where it is already near ground (Mendahara).
4. **Smoothing helps everywhere** that was tested, by 0.4 to 0.5 m of RMSE after calibration (Mendahara 2.74 to 2.32 m, Bengkalis 3.24 to 2.78 m, Sebangau 1.49 to 1.24 m).
5. **None of the sites meets the 0.25 m target** set in the re-peat methodology. The best smoothed RMSE is 1.24 m.
6. **Bengkalis problem classes.** The smoothed DEM has a -5.2 m bias in cropland (88 points) and -2.1 m in shrubland (52 points). Both classes have few points and I did not investigate them. They may come from masks and the 30 m smoothing near the coast or from the class k and intercept fits. Treat Bengkalis cropland and shrubland as unreliable.
7. **Heavy exclusion.** At Bengkalis 171,200 of 313,544 points fell in excluded land cover classes (mostly water and sea), and 75,119 more were near a land-cover edge. Only 15% were usable.

**Limits.** These runs are at 30 m, from DEMs averaged down in the case of Mendahara and Bengkalis. ICESat-2 only. No independent ground reference. Calibration and validation come from the same sites, split by 2 km blocks, which is not an independent test. The CHM and WorldCover dates differ from DEMNAS. Mendahara and Bengkalis include mangroves and coastal peat, where tide and the flat terrain mean 1-2 m errors are large relative to the relief.

## Related literature (Julzarika and others)

Searched on 2026-10-08. Entries marked with a dagger were seen only as a search summary or a title, not read in full. Check the full text before citing them.

- Julzarika, A. and Harintaka (2019). Free global DEM: converting DSM to DTM and its applications. *Int. Arch. Photogramm. Remote Sens. Spatial Inf. Sci.* XLII-4/W16, 319-325. doi 10.5194/isprs-archives-XLII-4-W16-319-2019. Abstract read. DSM from X-SAR, SRTM and TanDEM-X is integrated, referenced to EGM2008 and converted to DTM using canopy tilt angle, vegetation height and a radius around the canopy. Tested on Rote Island, not on peat.
- Julzarika, A., Aditya, T., Subaryono and Harintaka (2020). Comparison of the latest DTM with DEM Pleiades in monitoring the dynamic peatland. *ISRITI 2020*. doi 10.1109/ISRITI51436.2020.9315410. Abstract read. Palangkaraya-Pulang Pisau peatland, the nearest study area to Sebangau. The DTM is from ALOS PALSAR/PALSAR-2 InSAR and Sentinel DInSAR, referenced to EGM2008. Mean height differences: Pleiades vs DTM 0.923 m, Pleiades vs field GNSS 0.557 m, DTM vs field GNSS 0.705 m.
- Julzarika, A. Indonesian DEMNAS: DSM or DTM? ResearchGate (2020) (dagger). Directly on point, but the site blocked access, so the content is unread. This is the first thing to obtain.
- Julzarika's earlier height-model integration work (ALOS PALSAR, X-SAR, SRTM C, ICESat/GLAS) (dagger) reports vertical accuracies of about 1.1 to 2.0 m. Listed from a search summary only.
- Search summaries also say DEMNAS combines TerraSAR-X, IFSAR, ALOS PALSAR and mass points (dagger), that TerraSAR-X X-band stereo DSMs underestimate forest height by about 27% of canopy height because the signal penetrates the canopy (dagger), and give DEMNAS accuracy figures against GPS (8.17 m) and levelling (11.3 m standard deviation) (dagger) from studies I could not identify. None of these is verified.

**Why this matters here.** If DEMNAS is a mosaic of X-band, IFSAR and L-band sources, its canopy penetration and its vertical offset probably differ from tile to tile. That fits our finding that k is near zero and the offset varies, and it suggests testing whether the per-tile source or acquisition date explains the residual. I have not tested it.

## Limits of this evidence

- One 15 x 15 km window and one DEMNAS tile. Two satellite lidar sources (ICESat-2 ATL08 and GEDI L2A) that share the same geoid conversion, and no GNSS or airborne LiDAR points yet.
- A datum error common to both satellite sources cannot be separated from a DEM error without an independent reference (airborne LiDAR or GNSS).
- The canopy height product (Meta/WRI) and WorldCover date from different years than DEMNAS. This window burned in 2015 and 2019.
- The 0.25 m target in the re-peat methodology is not met by any run (best smoothed RMSE 1.24 m, on 441 points).
- The experiments tune and evaluate on the same window and the same ground points. The spatial split makes the comparison fair between runs, but it is not an independent test.

## Open items

- Obtain Julzarika's "Indonesian DEMNAS: DSM or DTM?" and any DEMNAS metadata on source sensors, acquisition years and the stated vertical datum.
- Add an independent ground reference: GNSS spot heights on canal banks and at loggers, or LiDAR, to decide whether the offset is a datum error or a DEM error.
- Run the airborne LiDAR comparison when the data arrive: it is the independent test of the datum question. Use it to check the shared EGM2008 conversion and the open-versus-forest height pattern.
- Test whether the DEM - ground residual maps to DEMNAS source tiles or acquisition dates.
- Report the relative (detrended) error alongside the absolute error, since the hydrology needs relative heights.
- Benchmark FABDEM against the same held-out points.
- Fail loudly (or add a flag) when SlideRule returns 0 points or errors: a failed fetch currently ends as a normal run with the unfitted prior k. The Bengkalis first run shows the cost.
- Add a timeout to GEDI requests: one hung for about 90 minutes at Mendahara.
- Rerun Mendahara with GEDI once SlideRule is stable, and check Bengkalis cropland and shrubland (bias -5.2 m and -2.1 m after smoothing).
- Rerun Mendahara and Bengkalis from the native 8 m DEM, not the 4x averaged one, if memory allows.
- Consider failing loudly on network errors in download mode.
- Cache downloaded canopy windows; the fetch took 4-5 min per tile.

## Mendahara at the native 8 m (2026-10-09)

Rerun of Mendahara from the full-resolution DEM (not the 4x averaged one), on an 8,456 x 8,830 px grid at 8 m, with canopy height and land cover fetched at 8 m. The ICESat-2 points from the 30 m run were reused (`ground_points_raw.gpkg` copied), so GEDI was not requested. Config: `runs/mendahara_8m.yaml`; output: `runs/mendahara_8m/` (637 MB). It completed in about 12 minutes on a machine with little free memory (1.2 GB free of 15 GB at the start).

| Run | Grid | Usable / loaded points | Validation n | Raw DEM RMSE (m) | Smoothed bias / RMSE / NMAD (m) |
|---|---|---|---|---|---|
| mendahara_30m | 30 m (4x averaged DEM) | 20,403 / 68,979 | 5,869 | 2.64 | +0.64 / 2.32 / 1.40 |
| mendahara_8m | 8 m native | 23,426 / 68,979 | 6,660 | 2.83 | +0.72 / 2.32 / 1.54 |

The validation sets differ because more points land on valid cells at 8 m, so the numbers are not strictly comparable. The smoothed RMSE is the same (2.32 m) at both grids. Fitted tree cover k is 0.09 (0.10 at 30 m). The finer grid does not improve the accuracy of `dtm_smooth.tif`. That surface is smoothed over 300 m and 150 m scales by design, so its effective resolution is set by the smoothing, not by the pixel size. The 8 m grid matters for `dtm_corrected.tif`, which keeps cell-level detail.
