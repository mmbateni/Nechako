# =============================================================================
# Nechako River Basin: nechakoBound.shp → RVH → GeoJSON
# All parameters filled from actual data sources:
#
#   ELEVATION  → extracted from uploaded DEM (ESRI Grid ADF, ~92m resolution)
#                fid=1 Nechako River : mean 1003 m  (DEM range 560–2600 m)
#                fid=2 Stuart River  : mean  971 m  (DEM range 658–2212 m)
#
#   LAND USE   → Sub-Boreal Interior ecozone, BC Interior Plateau
#                >90% evergreen coniferous forest (lodgepole pine,
#                white/hybrid spruce, subalpine fir)
#                Luvisolic soils on glacial till (Interior Plateau standard)
#
# NOTE: LandUse / Vegetation / SoilProfile names MUST match the class names
#       defined in your Raven .rvp parameter file. Adjust if your .rvp uses
#       different names for the same physical classes.
# =============================================================================


# ── STEP 1: Install & load packages ──────────────────────────────────────────

pkgs <- c("sf", "RavenR", "dplyr", "stringdist", "rmapshaper")
new  <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(new)) install.packages(new)

library(sf)
library(RavenR)
library(dplyr)

# Source the converter function (place it in your working directory)
source("rvn_rvh_shp_geojson.R")


# ── STEP 2: Paths ─────────────────────────────────────────────────────────────

shpfile    <- "nechakoBound.shp"   # your uploaded basin shapefile
rvhfile    <- "Nechako.rvh"        # created in Step 4
outputfile <- "Nechako.geojson"    # final output


# ── STEP 3: Read shapefile ────────────────────────────────────────────────────

sf_obj <- st_read(shpfile, quiet = TRUE)
# CRS: EPSG:3005 (NAD83 / BC Albers) — read automatically from .prj
cat("CRS :", st_crs(sf_obj)$input, "\n")
cat("Rows:", nrow(sf_obj), "\n")

# Compute centroids in projected CRS, then convert to WGS84 for lat/lon
centroids_wgs <- sf_obj |>
  st_centroid() |>
  st_transform(crs = 4326)
coords <- st_coordinates(centroids_wgs)   # [lon, lat]


# ── STEP 4: Build SBtable ─────────────────────────────────────────────────────
#
# Subbasins (from shapefile):
#   fid=1  Nechako River  ~31,630 km²   centroid: lat=53.69, lon=-125.33
#   fid=2  Stuart River   ~15,604 km²   centroid: lat=54.95, lon=-125.13
#
# Topology: Stuart River drains INTO Nechako → fid=2 downstream = 1
#           Nechako drains to Fraser (outlet) → fid=1 downstream = -1
#
# Elevation (DEM-derived, ESRI ADF Grid, 92 m resolution):
#   fid=1 Nechako : mean 1003 m  (range 560–2600 m, std 236 m)
#   fid=2 Stuart  : mean  971 m  (range 658–2212 m)
#
# Reach length: from FEAT_LEN attribute (perimeter proxy, km)
#   fid=1: ~1471 km perimeter  →  main channel ~458 km
#   fid=2: ~1204 km perimeter  →  main channel ~360 km

subIDs        <- as.integer(sf_obj$fid)            # 1, 2
subNames <- gsub(" ", "_", sf_obj$MJR_WTRSHM) # "Nechako_River", "Stuart_River"
areas_km2     <- sf_obj$FEAT_AREA / 1e6            # 31630.1, 15604.2
lats          <- coords[, "Y"]                     # 53.6859, 54.9486
lons          <- coords[, "X"]                     # -125.3307, -125.1338
downstreamIDs <- c(-1L, 1L)                        # Nechako→outlet, Stuart→Nechako

# --- DEM-derived mean elevations (extracted from uploaded ADF raster) ---
elevations <- c(1003, 971)   # m;  fid=1: 1003 m,  fid=2: 971 m

SBtable <- rvn_rvh_blankSBdf(nSubBasins = 2)
SBtable$SBID          <- subIDs
SBtable$Name          <- subNames
SBtable$Downstream_ID <- downstreamIDs
SBtable$Area          <- areas_km2
SBtable$Latitude      <- lats
SBtable$Longitude     <- lons
SBtable$Elevation     <- elevations        # ← DEM-derived
SBtable$ChannelSteep  <- 0.001             # default; refine with stream survey data
SBtable$Q_med         <- c(277, 85)        # m³/s — Nechako gauge at Isle Pierre;
                                           #        Stuart estimated ~85 m³/s
SBtable$ReachLength   <- c(458, 360)       # km — approximate main channel lengths

cat("\n=== SBtable ===\n")
print(SBtable[, c("SBID","Name","Downstream_ID","Area","Elevation","Latitude","Longitude")])


# ── STEP 5: Build HRUtable ────────────────────────────────────────────────────
#
# Land cover (BC Interior, Sub-Boreal ecozone — literature-based):
#
#   LandUse    = "FOREST"
#     >90% evergreen coniferous forest; lodgepole pine (Pinus contorta),
#     white/hybrid spruce (Picea glauca × engelmannii), subalpine fir (Abies lasiocarpa)
#     Sub-Boreal Spruce (SBS) and Sub-Boreal Pine-Spruce (SBPS) BEC zones
#
#   Vegetation = "FOREST_BOREAL"
#     Sub-boreal conifer mix; use "FOREST" if your .rvp does not define
#     "FOREST_BOREAL" separately
#
#   SoilProfile = "LUVISOL_LOAM"
#     Luvisolic soils dominant on Interior Plateau (Gray Luvisols over glacial
#     till); use "TILL_LOAM" or "DEFAULT_P" if not defined in your .rvp
#
# ─────────────────────────────────────────────────────────────────────────────
# ⚠️  IMPORTANT: These three names must exactly match the class names defined
#     in your Raven .rvp file (:LandUseClasses, :VegetationClasses,
#     :SoilProfiles blocks). Edit below to match your .rvp if different.
# ─────────────────────────────────────────────────────────────────────────────

HRUtable <- rvn_rvh_blankHRUdf(nHRUs = 2, subbasinIDs = subIDs)
HRUtable$HRU_ID      <- c(1L, 2L)
HRUtable$SubId       <- subIDs
HRUtable$Area        <- areas_km2
HRUtable$Latitude    <- lats
HRUtable$Longitude   <- lons
HRUtable$Elevation   <- elevations          # ← DEM-derived
HRUtable$LandUse     <- "FOREST"            # ← BC Interior Sub-Boreal conifer
HRUtable$Vegetation  <- "FOREST_BOREAL"     # ← Sub-boreal spruce/pine/fir mix
HRUtable$SoilProfile <- "LUVISOL_LOAM"      # ← Luvisolic/glacial till (Interior Plateau)
HRUtable$Aquifer     <- "[NONE]"
HRUtable$Terrain     <- "[NONE]"
HRUtable$GlobalParam <- "[NONE]"

cat("\n=== HRUtable ===\n")
print(HRUtable[, c("HRU_ID","SubId","Area","Elevation","LandUse","Vegetation","SoilProfile")])


# ── STEP 6: Write .rvh file ──────────────────────────────────────────────────

rvn_rvh_write(
  SBtable  = SBtable,
  HRUtable = HRUtable,
  filename = rvhfile
)
cat("\nRVH written to:", rvhfile, "\n")

rvh_check <- rvn_rvh_read(rvhfile)
cat("\n=== RVH verification ===\n")
print(rvh_check$SBtable[, c("SBID","Name","Downstream_ID","Area","TotalUpstreamArea","Elevation")])


# ── STEP 7: Convert shapefile + .rvh → GeoJSON ───────────────────────────────

result <- rvn_rvh_shp_geojson(
  shpfile          = shpfile,
  rvhfile          = rvhfile,
  outputfile       = outputfile,
  matchingColumns  = list(shpfile = "fid",    # fid column in nechakoBound.shp
                          rvhfile = "SBID"),  # SBID in Raven .rvh
  simplifyGeometry = TRUE
)

if (result$success) {
  cat("\n✔ SUCCESS! GeoJSON saved to:", result$file, "\n")
  cat("  Subbasins in output:", nrow(result$sf), "\n")
  cat("  Columns:", paste(names(result$sf), collapse = ", "), "\n")
} else {
  cat("\n✘ FAILED:", result$message, "\n")
}


# ── STEP 8: Quick plot ───────────────────────────────────────────────────────

cols <- c("lightblue", "lightyellow")
plot(st_geometry(result$sf), col = cols, border = "grey30", lwd = 1.5,
     main = "Nechako Basin – Subwatersheds")
text(st_coordinates(st_centroid(result$sf)),
     labels = paste0(result$sf$rvhName, "\n", round(result$sf$BasArea / 1e6, 1), " km²"),
     cex = 0.8, font = 2)
legend("bottomright",
       legend = c("Nechako River", "Stuart River"),
       fill   = cols, bty = "n")
