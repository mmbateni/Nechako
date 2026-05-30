#' Convert a basin shapefile and RVH to GeoJSON (robust, sf-based, improved)
#'
#' This function reads an input basin shapefile and an RVH file (via RavenR::rvn_rvh_read()),
#' matches subbasin IDs between them (with fuzzy matching fallback), copies selected RVH
#' attributes to the shapefile attribute table, computes or applies outlet coordinates, optionally
#' simplifies geometry, and writes a GeoJSON file.
#'
#' @param shpfile path to basin shapefile (any format readable by sf::st_read)
#' @param rvhfile path to RVH file (read by RavenR::rvn_rvh_read)
#' @param outputfile path to output GeoJSON (defaults to ./output.geojson)
#' @param CRSshp optional PROJ4 or EPSG integer to assign if input shapefile has no CRS
#' @param matchingColumns list with elements shpfile and rvhfile giving column names to match on
#' @param outletCoords list with optional elements: ID (vector of SubId to which coords apply), outletLat, outletLng
#' @param simplifyGeometry logical (or numeric keep fraction). If TRUE uses rmapshaper::ms_simplify(keep = 0.05). If a number in (0,1) uses that as keep.
#' @param missing_value value used for missing numeric fields (default -9999). Use NA to keep NA.
#' @param overwrite logical whether to overwrite existing outputfile (default TRUE)
#' @return A list with: $success (logical), $file (path or NULL), $sf (sf object written), $message
#' @examples
#' # rvn_rvh_shp_geojson('basins.shp','model.rvh','basins.geojson')
#' @export
rvn_rvh_shp_geojson <- function(
    shpfile,
    rvhfile,
    outputfile = file.path(getwd(), "output.geojson"),
    CRSshp = NA,
    matchingColumns = list(shpfile = "subid", rvhfile = "SBID"),
    outletCoords = list(ID = NA, outletLat = NA, outletLng = NA),
    simplifyGeometry = TRUE,
    missing_value = -9999,
    overwrite = TRUE
) {
  # ---- Dependencies ----
  required <- c("sf", "RavenR", "stringdist", "rmapshaper", "dplyr")
  missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) {
    stop("Please install required packages: ", paste(missing_pkgs, collapse = ", "))
  }
  # load quietly
  library(sf)
  library(RavenR)
  library(stringdist)
  library(rmapshaper)
  library(dplyr)
  
  # ---- Basic checks ----
  if (missing(shpfile) || missing(rvhfile)) stop("Arguments 'shpfile' and 'rvhfile' are required")
  
  if (!file.exists(shpfile)) stop("Shapefile does not exist: ", shpfile)
  if (!file.exists(rvhfile)) stop("RVH file does not exist: ", rvhfile)
  outdir <- dirname(outputfile)
  if (!dir.exists(outdir)) stop("Output directory does not exist: ", outdir)
  if (file.exists(outputfile) && !overwrite) stop("Output file exists and overwrite = FALSE: ", outputfile)
  
  # ---- Read shapefile (sf) ----
  sf_obj <- tryCatch({
    st_read(shpfile, quiet = TRUE)
  }, error = function(e) stop("Failed to read shapefile: ", e$message))
  
  # assign CRS if missing
  if (is.na(st_crs(sf_obj))) {
    if (!is.na(CRSshp)) {
      st_crs(sf_obj) <- CRSshp
    } else stop("Input shapefile has no CRS and CRSshp not supplied")
  }
  # transform to WGS84
  sf_obj <- st_make_valid(st_transform(sf_obj, crs = 4326))

  # ---- Read RVH ----
  rvh <- tryCatch({
    rvn_rvh_read(rvhfile)
  }, error = function(e) stop("Failed to read RVH file with RavenR::rvn_rvh_read(): ", e$message))
  
  if (!is.data.frame(rvh$SBtable)) stop("Unexpected RVH structure: SBtable missing or not a data.frame")
  
  sb_tbl <- rvh$SBtable
  
  # ---- Identify columns to match (with fuzzy fallback) ----
  shp_match_requested <- matchingColumns$shpfile
  rvh_match_requested <- matchingColumns$rvhfile
  
  find_col <- function(want, choices, side = "shp") {
    if (is.null(want) || is.na(want) || want == "") return(NA_character_)
    if (want %in% choices) return(want)
    # try case-insensitive exact match
    ci <- choices[tolower(choices) == tolower(want)]
    if (length(ci)) return(ci[1])
    # fuzzy match via stringdist
    idx <- amatch(want, choices, maxDist = 3)
    if (!is.na(idx)) return(choices[idx])
    # try common fallbacks
    fallbacks <- list(
      shp = c("subid", "subbasinid", "id", "gid"),
      rvh = c("SBID", "SubId", "subid", "ID")
    )
    for (f in fallbacks[[side]]) {
      ci <- choices[tolower(choices) == tolower(f)]
      if (length(ci)) return(ci[1])
    }
    return(NA_character_)
  }
  
  shp_col <- find_col(shp_match_requested, names(sf_obj), side = "shp")
  rvh_col <- find_col(rvh_match_requested, names(sb_tbl), side = "rvh")
  
  if (is.na(shp_col) || is.na(rvh_col)) {
    stop("Could not determine matching columns. Found shp: ", shp_col, "; rvh: ", rvh_col,
         ". Please set matchingColumns explicitly to existing column names.")
  }
  
  # coerce both to character for matching
  sf_obj <- sf_obj %>% mutate(.match_key = as.character(.data[[shp_col]]))
  sf_obj<-sf_obj[,".match_key"]
  sf_obj$BasArea<-as.vector(st_area(sf_obj))
  sb_tbl[[rvh_col]] <- as.character(sb_tbl[[rvh_col]])
  
  # ---- Join attributes ----
  # select relevant columns from rvh SBtable (keep SubId, Downstream_ID, Name, Area if present)
  rvh_columns <- c(.match_key = rvh_col,
                   DowSubId = ifelse("Downstream_ID" %in% names(sb_tbl), "Downstream_ID",
                                     ifelse("DownstreamID" %in% names(sb_tbl), "DownstreamID", NA)),
                   rvhName = ifelse("Name" %in% names(sb_tbl), "Name", NA))
  # ensure columns exist
  rvh_columns <- rvh_columns[!is.na(rvh_columns)]
  sb_sel <- sb_tbl[, unique(unname(rvh_columns)), drop = FALSE]
  colnames(sb_sel)<-names(rvh_columns)
  sb_sel$.match_key <- as.character(sb_sel$.match_key)
  
  # left_join
  sf_joined <- left_join(sf_obj, sb_sel, by = ".match_key")
  
  # replace NA values for numeric fields with missing_value if requested
  if (!is.na(missing_value)) {
    num_cols <- sapply(sf_joined, is.numeric)
    sf_joined[num_cols] <- lapply(sf_joined[num_cols], function(x) ifelse(is.na(x), missing_value, x))
  }
  
  # ---- Outlet coordinates handling ----
  # Try to build outlet lat/lon per subbasin from rvh HRU table if not provided
  outlet_IDs <- outletCoords$ID
  outletLat <- outletCoords$outletLat
  outletLng <- outletCoords$outletLng
  
  # function to extract HRU lat/lon and aggregate by SubId
  hru_latlon_by_sub <- function(rvh) {
    # try to find Hru table (commonly rvh[[2]] or rvh$HRUtable)
    hru <- NULL
    if (is.data.frame(rvh$HRUtable)) hru <- rvh$HRUtable
    if (is.data.frame(rvh[[2]]) && is.null(hru)) hru <- rvh[[2]]
    if (is.null(hru)) return(NULL)
    # find lat/lon columns
    lat_col <- find_col("Latitude", names(hru), side = "rvh")
    lon_col <- find_col("Longitude", names(hru), side = "rvh")
    id_col <- find_col("SubId", names(hru), side = "rvh")
    if (any(is.na(c(lat_col, lon_col, id_col)))) return(NULL)
    hru[[lat_col]] <- as.numeric(as.character(hru[[lat_col]]))
    hru[[lon_col]] <- as.numeric(as.character(hru[[lon_col]]))
    agg <- hru %>% group_by(.data[[id_col]]) %>% summarize(outletLat = mean(.data[[lat_col]], na.rm = TRUE),
                                                           outletLng = mean(.data[[lon_col]], na.rm = TRUE)) %>%
      ungroup()
    names(agg)[1] <- ".match_key"
    agg$.match_key <- as.character(agg$.match_key)
    agg
  }
  

  # prepare outlet table with .match_key
  outlet_tbl <- NULL
  if (!all(is.na(outlet_IDs))) {
    # user provided IDs; match to .match_key
    outlet_tbl <- data.frame(.match_key = as.character(outlet_IDs), outletLat = outletLat, outletLng = outletLng, stringsAsFactors = FALSE)
  } else if (!all(is.na(outletLat)) && !all(is.na(outletLng)) && length(outletLat) == nrow(sf_joined)) {
    # user provided vector aligned to sf rows
    outlet_tbl <- data.frame(.match_key = sf_joined$.match_key, outletLat = outletLat, outletLng = outletLng, stringsAsFactors = FALSE)
  } else {
    outlet_tbl <- hru_latlon_by_sub(rvh)
  }
  
  if (!is.null(outlet_tbl)) {
    sf_joined <- left_join(sf_joined, outlet_tbl, by = ".match_key")
    # if missing and missing_value specified, replace
    if (!is.na(missing_value)) {
      if (!"outletLat" %in% names(sf_joined)) sf_joined$outletLat <- missing_value
      if (!"outletLng" %in% names(sf_joined)) sf_joined$outletLng <- missing_value
      sf_joined$outletLat[is.na(sf_joined$outletLat)] <- missing_value
      sf_joined$outletLng[is.na(sf_joined$outletLng)] <- missing_value
    }
  } else {
    # no outlet info found -- create columns
    sf_joined$outletLat <- ifelse(is.na(missing_value), NA_real_, missing_value)
    sf_joined$outletLng <- ifelse(is.na(missing_value), NA_real_, missing_value)
    warning("No outlet coordinates could be determined from inputs or RVH HRU table. outletLat/outletLng set to missing_value")
  }
  
  # ---- Optional geometry simplification ----
  geom_to_write <- sf_joined
  if (!identical(simplifyGeometry, FALSE)) {
    if (is.logical(simplifyGeometry) && simplifyGeometry) {
      geom_to_write <- ms_simplify(sf_joined, keep = 0.05, keep_shapes = TRUE)
    } else if (is.numeric(simplifyGeometry) && simplifyGeometry > 0 && simplifyGeometry < 1) {
      geom_to_write <- ms_simplify(sf_joined, keep = simplifyGeometry, keep_shapes = TRUE)
    }
  }
  colnames(geom_to_write)[colnames(geom_to_write)==".match_key"]<-"SubId"
  # ---- Write GeoJSON ----
  # ensure extension
  if (!grepl("\\.geojson$", outputfile, ignore.case = TRUE)) {
    outputfile <- paste0(tools::file_path_sans_ext(outputfile), ".geojson")
  }
  
  # remove existing if overwrite
  if (file.exists(outputfile) && overwrite) file.remove(outputfile)
  
  tryCatch({
    st_write(geom_to_write, outputfile, driver = "GeoJSON", delete_dsn = overwrite, quiet = TRUE)
  }, error = function(e) stop("Failed to write GeoJSON: ", e$message))
  
  success <- file.exists(outputfile)
  message <- if (success) "Successfully written GeoJSON" else "Failed to write GeoJSON"
  
  return(list(success = success, file = if (success) normalizePath(outputfile) else NULL, sf = geom_to_write, message = message))
}
