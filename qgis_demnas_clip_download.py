# -*- coding: utf-8 -*-
"""
Download & clip DEMNAS (Badan Informasi Geospasial) - QGIS 4 Processing script

The public DEMNAS ImageServer (geoservices.big.go.id) returns a placeholder
value (1.134e38) for pixels to anonymous users, so its exportImage/identify
cannot be used for real elevations. This script instead downloads the actual
DEMNAS GeoTIFF tiles from BIG's Ina-Geoportal (tanahair.indonesia.go.id),
then mosaics and clips them to a polygon.

How it works
  1. You paste ONE fresh DEMNAS download link from the Ina-Geoportal. The
     link carries a JWT token (valid ~1 hour). The token is NOT tied to a
     filename, so the script reuses it for every tile.
  2. The polygon is dissolved and transformed to EPSG:4326.
  3. The ImageServer's footprint catalog (the /query operation - vector
     metadata, which still works) is asked which tiles intersect the area.
     Each footprint gives the map-sheet number (NLP) -> the tile filename.
  4. Each tile is downloaded from tanahair with your token, reused across
     tiles, with a polite delay between requests.
  5. Tiles are mosaicked and clipped to the polygon (optionally reprojected).

Getting the link (needed once per run, token lasts ~1 hour)
  - Log in at https://tanahair.indonesia.go.id/, open the DEMNAS download map,
    click any one tile, and copy the download URL it opens, e.g.
    https://tanahair.indonesia.go.id/api-inageo/unduh/demnas?token=eyJ...&filename=DEMNAS_0917-14_v1.0.tif
  - Paste that whole URL into the "DEMNAS download link" box below. The exact
    tile you picked does not matter; the script works out the right tiles.

Note: you are downloading data you already have free access to, using your own
session token. Be considerate with the request delay on BIG's servers.

Data source: Badan Informasi Geospasial (DEMNAS). Vertical datum EGM2008.
"""

import base64
import json
import os
import re
import tempfile
import time
from urllib.parse import parse_qs, urlencode, urlparse

from osgeo import gdal

from qgis.PyQt.QtCore import QCoreApplication, QEventLoop, QUrl
from qgis.PyQt.QtNetwork import QNetworkRequest
from qgis.core import (
    Qgis,
    QgsBlockingNetworkRequest,
    QgsCoordinateReferenceSystem,
    QgsCoordinateTransform,
    QgsFeature,
    QgsFeatureRequest,
    QgsFileDownloader,
    QgsGeometry,
    QgsPointXY,
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterCrs,
    QgsProcessingParameterEnum,
    QgsProcessingParameterExtent,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterFolderDestination,
    QgsProcessingParameterNumber,
    QgsProcessingParameterRasterDestination,
    QgsProcessingParameterString,
    QgsProcessingUtils,
    QgsRasterFileWriter,
    QgsVectorLayer,
)

gdal.UseExceptions()

SERVICE_URL = (
    "https://geoservices.big.go.id/raster/rest/services/"
    "DEMNAS/DEM_Indonesia/ImageServer"
)
TIFF_MAGIC = (b"II*\x00", b"MM\x00*", b"II+\x00", b"MM\x00+")
NLP_RE = re.compile(r"(\d{4}-\d{2,3})")  # 1:50k = 2 digits, 1:25k = 3 digits


class DemnasClipDownload(QgsProcessingAlgorithm):
    INPUT = "INPUT"
    EXTENT = "EXTENT"
    LINK = "LINK"
    CLIP = "CLIP"
    TARGET_CRS = "TARGET_CRS"
    RESAMPLING = "RESAMPLING"
    VERSION = "VERSION"
    DELAY = "DELAY"
    NODATA = "NODATA"
    KEEP = "KEEP"
    SERVICE = "SERVICE"
    OUTPUT = "OUTPUT"

    RESAMPLING_LABELS = ["Nearest neighbour", "Bilinear", "Cubic"]
    RESAMPLING_GDAL = ["near", "bilinear", "cubic"]

    # ------------------------------------------------------------------ meta
    def tr(self, text):
        return QCoreApplication.translate("DemnasClipDownload", text)

    def createInstance(self):
        return DemnasClipDownload()

    def name(self):
        return "demnas_clip_download"

    def displayName(self):
        return self.tr("Download & clip DEMNAS (BIG / Ina-Geoportal)")

    def group(self):
        return self.tr("BIG / DEMNAS")

    def groupId(self):
        return "big_demnas"

    def shortHelpString(self):
        return self.tr(
            "Downloads the real DEMNAS tiles from BIG's Ina-Geoportal "
            "(tanahair.indonesia.go.id) for the area of a polygon and clips "
            "them to it.\n\n"
            "You must paste ONE fresh DEMNAS download link (it carries a "
            "login token valid for about 1 hour). The token is reused for "
            "every tile, so you only need one.\n\n"
            "Get it: log in at tanahair.indonesia.go.id, open the DEMNAS "
            "download map, click any single tile, and copy the URL that "
            "opens. Any tile works - the script finds the right ones.\n\n"
            "Area of interest: give EITHER a boundary polygon OR an extent "
            "(the extent box lets you draw a rectangle on the canvas, take a "
            "layer's or the map's extent, or type coordinates). If you want a "
            "free-form shape rather than a rectangle, draw a temporary "
            "scratch polygon layer and use it as the boundary polygon.\n\n"
            "Values are heights in metres (vertical datum EGM2008)."
        )

    # ------------------------------------------------------------ parameters
    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT,
                self.tr("Boundary polygon (use this OR the extent below)"),
                [Qgis.ProcessingSourceType.VectorPolygon],
                optional=True,
            )
        )
        self.addParameter(
            QgsProcessingParameterExtent(
                self.EXTENT,
                self.tr("Or area by extent (draw on canvas / layer / typed)"),
                optional=True,
            )
        )
        self.addParameter(
            QgsProcessingParameterString(
                self.LINK,
                self.tr("DEMNAS download link (fresh, from Ina-Geoportal)"),
            )
        )
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.CLIP,
                self.tr("Clip to polygon shape (unchecked = keep whole tiles)"),
                defaultValue=True,
            )
        )
        self.addParameter(
            QgsProcessingParameterCrs(
                self.TARGET_CRS,
                self.tr("Reproject output to (empty = keep native EPSG:4326)"),
                optional=True,
            )
        )

        advanced = [
            QgsProcessingParameterEnum(
                self.RESAMPLING,
                self.tr("Resampling (only used when reprojecting)"),
                options=self.RESAMPLING_LABELS,
                defaultValue=1,
            ),
            QgsProcessingParameterNumber(
                self.DELAY,
                self.tr("Pause between downloads (seconds)"),
                type=Qgis.ProcessingNumberParameterType.Double,
                defaultValue=2.0,
                minValue=0.0,
            ),
            QgsProcessingParameterString(
                self.VERSION,
                self.tr("Version suffix in filenames"),
                defaultValue="v1.0",
            ),
            QgsProcessingParameterNumber(
                self.NODATA,
                self.tr("Output NoData value"),
                type=Qgis.ProcessingNumberParameterType.Double,
                defaultValue=-9999.0,
            ),
            QgsProcessingParameterFolderDestination(
                self.KEEP,
                self.tr("Also keep the downloaded tiles in this folder"),
                optional=True,
                createByDefault=False,
            ),
            QgsProcessingParameterString(
                self.SERVICE,
                self.tr("ImageServer URL (tile index)"),
                defaultValue=SERVICE_URL,
            ),
        ]
        for p in advanced:
            p.setFlags(p.flags() | Qgis.ProcessingParameterFlag.Advanced)
            self.addParameter(p)

        self.addParameter(
            QgsProcessingParameterRasterDestination(
                self.OUTPUT, self.tr("DEMNAS clipped")
            )
        )

    # --------------------------------------------------------------- helpers
    def _get(self, url, feedback, retries=3):
        """GET via QGIS network stack (respects proxy/SSL settings)."""
        last = ""
        for attempt in range(1, retries + 1):
            if feedback.isCanceled():
                raise QgsProcessingException(self.tr("Canceled"))
            wait = self._delay - (time.monotonic() - self._last)
            if wait > 0:
                time.sleep(wait)
            self._last = time.monotonic()
            req = QgsBlockingNetworkRequest()
            err = req.get(QNetworkRequest(QUrl(url)), True, feedback)
            if err == QgsBlockingNetworkRequest.ErrorCode.NoError:
                return req.reply().content().data()
            last = req.errorMessage()
            feedback.pushWarning(
                self.tr(f"Request failed (attempt {attempt}/{retries}): {last}")
            )
            time.sleep(4 * attempt)
        raise QgsProcessingException(self.tr(f"Request failed: {last}"))

    def _parse_link(self, link, feedback):
        """Pull token + download endpoint out of a pasted URL."""
        link = link.strip().strip('"').strip("'")
        u = urlparse(link)
        if not u.scheme or "unduh" not in u.path:
            raise QgsProcessingException(
                self.tr("That does not look like an Ina-Geoportal download link.")
            )
        qs = parse_qs(u.query)
        token = (qs.get("token") or [""])[0]
        if not token:
            raise QgsProcessingException(self.tr("No token found in the link."))
        endpoint = f"{u.scheme}://{u.netloc}{u.path}"

        # Decode the (unverified) JWT payload just to report time left.
        try:
            pl = token.split(".")[1]
            pl += "=" * (-len(pl) % 4)
            exp = json.loads(base64.urlsafe_b64decode(pl)).get("exp")
            if exp:
                left = exp - time.time()
                if left <= 0:
                    raise QgsProcessingException(
                        self.tr("This token has expired - paste a fresh link.")
                    )
                feedback.pushInfo(
                    self.tr(f"Token valid for about {left / 60:.0f} more minute(s).")
                )
        except QgsProcessingException:
            raise
        except Exception:  # noqa: BLE001
            feedback.pushWarning(self.tr("Could not read token expiry; continuing."))

        sample = (qs.get("filename") or [""])[0]
        feedback.pushInfo(self.tr(f"Example filename in link: {sample or '(none)'}"))
        return token, endpoint

    def _query_tiles(self, service, bbox, feedback):
        """Ask the ImageServer footprint catalog for tiles intersecting bbox.

        Returns a list of (nlp, QgsGeometry) for primary rasters.
        """
        params = {
            "where": "Category = 1",  # 1 = Primary (skip overviews)
            "geometry": f"{bbox.xMinimum()},{bbox.yMinimum()},"
                        f"{bbox.xMaximum()},{bbox.yMaximum()}",
            "geometryType": "esriGeometryEnvelope",
            "inSR": 4326,
            "outSR": 4326,
            "spatialRel": "esriSpatialRelIntersects",
            "outFields": "Name",
            "returnGeometry": "true",
            "f": "json",
        }
        url = f"{service}/query?{urlencode(params)}"
        js = json.loads(self._get(url, feedback).decode("utf-8"))
        if "error" in js:
            raise QgsProcessingException(
                self.tr(f"Tile-index query failed: {js['error']}")
            )
        out, seen = [], set()
        for f in js.get("features", []):
            name = str(f.get("attributes", {}).get("Name", ""))
            m = NLP_RE.search(name)
            if not m or m.group(1) in seen:
                continue
            rings = f.get("geometry", {}).get("rings")
            if not rings:
                continue
            poly = [[QgsPointXY(x, y) for x, y in ring] for ring in rings]
            out.append((m.group(1), QgsGeometry.fromPolygonXY(poly)))
            seen.add(m.group(1))
        if not out:
            raw = [str(f.get("attributes", {}).get("Name", "")) for f in js.get("features", [])]
            feedback.pushInfo(self.tr(f"Raw Name values returned: {raw[:20]}"))
            raise QgsProcessingException(
                self.tr(
                    "No map-sheet numbers found in the tile index for this "
                    "area. See the raw Name values just logged."
                )
            )
        feedback.pushInfo(self.tr("Tiles from index: " + ", ".join(sorted(seen))))
        return out

    def _download(self, endpoint, token, filename, dest, feedback, base, span):
        """Stream a tile straight to disk (no full-file RAM buffer)."""
        url = f"{endpoint}?{urlencode({'token': token, 'filename': filename})}"
        wait = self._delay - (time.monotonic() - self._last)
        if wait > 0:
            time.sleep(wait)
        self._last = time.monotonic()

        tmp = dest + ".part"
        if os.path.exists(tmp):
            try:
                os.remove(tmp)
            except OSError:
                pass

        loop = QEventLoop()
        state = {"ok": False, "err": []}
        try:
            dl = QgsFileDownloader(QUrl(url), tmp, delayStart=True)
        except TypeError:  # older constructor signature
            dl = QgsFileDownloader(QUrl(url), tmp, "", True)
        dl.downloadCompleted.connect(lambda *a: state.update(ok=True))
        dl.downloadError.connect(lambda errs: state["err"].extend(list(errs)))
        dl.downloadExited.connect(loop.quit)

        def _p(received, total):
            if feedback.isCanceled():
                dl.cancelDownload()
                return
            if span and total > 0:
                feedback.setProgress(base + span * received / total)

        dl.downloadProgress.connect(_p)
        dl.startDownload()
        loop.exec()

        def _cleanup():
            if os.path.exists(tmp):
                try:
                    os.remove(tmp)
                except OSError:
                    pass

        if feedback.isCanceled():
            _cleanup()
            raise QgsProcessingException(self.tr("Canceled"))
        if not state["ok"] or state["err"]:
            _cleanup()
            raise QgsProcessingException(
                self.tr(
                    f"{filename}: download failed: "
                    f"{'; '.join(state['err']) or 'unknown error'}"
                )
            )
        # Guard against an HTML/JSON error page saved as a .tif
        with open(tmp, "rb") as fh:
            head = fh.read(4)
        if head not in TIFF_MAGIC:
            with open(tmp, "rb") as fh:
                snippet = fh.read(300).decode("utf-8", "replace")
            _cleanup()
            raise QgsProcessingException(
                self.tr(
                    f"{filename}: not a GeoTIFF (token expired or tile "
                    f"missing?). Response: {snippet}"
                )
            )
        os.replace(tmp, dest)

    # ------------------------------------------------------------- algorithm
    def processAlgorithm(self, parameters, context, feedback):
        source = self.parameterAsSource(parameters, self.INPUT, context)
        link = self.parameterAsString(parameters, self.LINK, context)
        clip = self.parameterAsBoolean(parameters, self.CLIP, context)
        target_crs = self.parameterAsCrs(parameters, self.TARGET_CRS, context)
        resampling_idx = self.parameterAsEnum(parameters, self.RESAMPLING, context)
        version = self.parameterAsString(parameters, self.VERSION, context).strip()
        nodata = self.parameterAsDouble(parameters, self.NODATA, context)
        keep = self.parameterAsString(parameters, self.KEEP, context)
        service = self.parameterAsString(parameters, self.SERVICE, context).rstrip("/")
        out_path = self.parameterAsOutputLayer(parameters, self.OUTPUT, context)

        self._delay = self.parameterAsDouble(parameters, self.DELAY, context)
        self._last = 0.0

        token, endpoint = self._parse_link(link, feedback)

        wgs84 = QgsCoordinateReferenceSystem("EPSG:4326")
        reproject = target_crs.isValid() and target_crs != wgs84

        # 1. Build the area of interest in EPSG:4326 ------------------------
        rect = self.parameterAsExtent(parameters, self.EXTENT, context, wgs84)
        if source is not None:
            if not rect.isEmpty():
                feedback.pushInfo(
                    self.tr("Both polygon and extent given; using the polygon.")
                )
            xform = QgsCoordinateTransform(
                source.sourceCrs(), wgs84, context.transformContext()
            )
            geoms = []
            for feat in source.getFeatures(QgsFeatureRequest().setNoAttributes()):
                if feedback.isCanceled():
                    return {}
                g = feat.geometry()
                if g is None or g.isEmpty():
                    continue
                g = QgsGeometry(g)
                g.transform(xform)
                geoms.append(g)
            if not geoms:
                raise QgsProcessingException(self.tr("Input layer has no geometries."))
            mask = QgsGeometry.unaryUnion(geoms)
            if not mask.isGeosValid():
                mask = mask.makeValid()
        elif not rect.isEmpty():
            mask = QgsGeometry.fromRect(rect)
        else:
            raise QgsProcessingException(
                self.tr("Provide either a boundary polygon or an extent.")
            )
        bbox = mask.boundingBox()
        feedback.pushInfo(self.tr(f"Area extent (EPSG:4326): {bbox.toString(6)}"))

        # 2. Find the tiles that actually touch the polygon ------------------
        tiles = self._query_tiles(service, bbox, feedback)
        tiles = [(nlp, g) for nlp, g in tiles if g.intersects(mask)]
        if not tiles:
            raise QgsProcessingException(self.tr("No tiles intersect the polygon."))
        feedback.pushInfo(
            self.tr(f"{len(tiles)} tile(s) to download: "
                    + ", ".join(sorted(n for n, _ in tiles)))
        )

        # 3. Download each tile from tanahair (one token, reused) -----------
        work = tempfile.mkdtemp(prefix="demnas_", dir=QgsProcessingUtils.tempFolder())
        save_dir = keep if keep else work
        if keep:
            os.makedirs(keep, exist_ok=True)
        rasters = []
        for i, (nlp, _) in enumerate(sorted(tiles)):
            if feedback.isCanceled():
                return {}
            filename = f"DEMNAS_{nlp}_{version}.tif"
            dest = os.path.join(save_dir, filename)
            n = len(tiles)
            if keep and os.path.exists(dest) and os.path.getsize(dest) > 0:
                feedback.pushInfo(
                    self.tr(f"Tile {i + 1}/{n}: {filename} (cached)")
                )
            else:
                feedback.pushInfo(self.tr(f"Tile {i + 1}/{n}: {filename}"))
                self._download(
                    endpoint, token, filename, dest, feedback,
                    base=80.0 * i / n, span=80.0 / n,
                )
            rasters.append(dest)

        # 4. Mosaic the tiles, and give it a CRS ----------------------------
        # DEMNAS tiles from tanahair ship WITHOUT an embedded CRS, so the VRT
        # has none either. We assign EPSG:4326 (DEMNAS is geographic), which
        # is what previously broke the clip ("Cannot find source SRS").
        mosaic = os.path.join(work, "mosaic.vrt")
        gdal.BuildVRT(mosaic, rasters)
        ds = gdal.Open(mosaic, gdal.GA_Update)
        if ds is not None and not ds.GetProjectionRef():
            ds.SetProjection(wgs84.toWkt())
            feedback.pushInfo(
                self.tr("Tiles had no CRS; assigned EPSG:4326 to the mosaic.")
            )
        ds = None

        # 5. Clip / reproject with QGIS Processing (GDAL provider) ----------
        # Passing the mask as a real layer lets QGIS handle the cutline CRS,
        # instead of hand-writing GeoTIFF/cutline SRS ourselves.
        import processing

        mask_layer = QgsVectorLayer("Polygon?crs=EPSG:4326", "aoi", "memory")
        mf = QgsFeature()
        mf.setGeometry(mask)
        mask_layer.dataProvider().addFeature(mf)
        mask_layer.updateExtents()

        gtiff_opts = "COMPRESS=DEFLATE|PREDICTOR=3|TILED=YES|BIGTIFF=IF_SAFER"
        stage_in = mosaic

        feedback.pushInfo(self.tr("Mosaicking / clipping / writing output..."))
        if clip:
            dest = out_path if not reproject else os.path.join(work, "clip.tif")
            processing.run(
                "gdal:cliprasterbymasklayer",
                {
                    "INPUT": stage_in,
                    "MASK": mask_layer,
                    "SOURCE_CRS": wgs84,
                    "TARGET_CRS": None,
                    "NODATA": nodata,
                    "CROP_TO_CUTLINE": True,
                    "KEEP_RESOLUTION": True,
                    "MULTITHREADING": True,
                    "DATA_TYPE": 0,  # keep Float32
                    "OPTIONS": gtiff_opts,
                    "OUTPUT": dest,
                },
                context=context, feedback=feedback, is_child_algorithm=True,
            )
            stage_in = dest

        if reproject:
            processing.run(
                "gdal:warpreproject",
                {
                    "INPUT": stage_in,
                    "SOURCE_CRS": wgs84,
                    "TARGET_CRS": target_crs,
                    "RESAMPLING": resampling_idx,
                    "NODATA": nodata,
                    "DATA_TYPE": 0,
                    "MULTITHREADING": True,
                    "OPTIONS": gtiff_opts,
                    "OUTPUT": out_path,
                },
                context=context, feedback=feedback, is_child_algorithm=True,
            )
        elif not clip:
            # neither clip nor reproject: just write the mosaic out
            processing.run(
                "gdal:translate",
                {
                    "INPUT": stage_in,
                    "DATA_TYPE": 0,
                    "OPTIONS": gtiff_opts,
                    "OUTPUT": out_path,
                },
                context=context, feedback=feedback, is_child_algorithm=True,
            )

        feedback.setProgress(100)
        if keep:
            feedback.pushInfo(self.tr(f"Downloaded tiles kept in: {keep}"))
        return {self.OUTPUT: out_path}
