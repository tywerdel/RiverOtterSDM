######################################################################
# Species Distribution Model (SDM) pipeline — GBIF → Maxnet → Scenarios
#
# AUTHOR / AFFILIATION
#   Ty J. Werdel (ORCID: 0000-0003-4023-9668)
#   Department of Rangeland, Wildlife, and Fisheries Management
#   Texas A&M University, College Station, TX, USA
#
# PURPOSE
#   One script that can run start-to-finish to reproduce the key
#   manuscript-reported results (and to regenerate all intermediate files).
#
# WHAT THIS SCRIPT PRINTS (for manuscript cross-check)
#   - Cleaned independent GBIF records (n)
#   - Counties with ≥1 record (n) + new-county counts by first-detection period
#   - Presence points before/after spatial thinning (n)
#   - K-fold CV tuning summary + selected (fc, rm)
#   - Final omit-10% threshold (cloglog)
#   - Scenario raster paths written (current + noHuman) + threshold reminder
#   - Predictor importance Σ|β| by predictor (table)
#
# NOTES
#   - Step B uses EPSG:5070 + 100 m buffer join (robust county join).
#   - Step H uses a stable terra tempdir to avoid missing temp-file errors.
#   - GIS post-processing (binary maps, areas, connectivity) can remain in ArcGIS Pro.
######################################################################


# --------------------------- 0. HOUSEKEEPING --------------------------
rm(list = ls()); gc()
if (!is.null(dev.list())) dev.off()
closeAllConnections()
cat("\014")  # clear console (RStudio)

# --------------------------- 1. USER PARAMETERS ------------------------
species_name   <- "Lontra canadensis"
region_state   <- "TX"
year_min       <- 1900
year_max       <- 2025

thin_km        <- 10
k_folds        <- 5
seed_value     <- 42
n_background   <- 10000

# Maxnet tuning grid (feature classes must be lowercase for maxnet)
fc_grid <- c("l", "lq", "h", "lqh", "lqhp")
rm_grid <- c(0.5, 1.0, 1.5, 2.0, 3.0)

# overwrite = TRUE ensures full regeneration of all outputs for reproducibility
# overwrite = FALSE = reuse existing outputs; STILL prints results from disk.
overwrite      <- FALSE

# Equal-area CRS used for thinning + county join buffering
crs_area       <- "EPSG:5070"   # NAD83 / Conus Albers

# --------------------------- 2. PATHS / I-O ----------------------------
gbif_file <- "GBIF_Observations_10_1_25_Texas.xlsx"

nhd_dir   <- "00_raw/nhd"
natwat_focal <- file.path(nhd_dir, "natwatfocal.tif")
humwat_focal <- file.path(nhd_dir, "humwatfocal.tif")
allwat_focal <- file.path(nhd_dir, "allwatfocal.tif")

major_rivers_shp <- "00_raw/MajorRivers_TX/MajorRivers_TX.shp"

dir.create("01_intermediate", showWarnings = FALSE, recursive = TRUE)
dir.create("02_occ",          showWarnings = FALSE, recursive = TRUE)
dir.create("03_models",       showWarnings = FALSE, recursive = TRUE)
dir.create("04_pred",         showWarnings = FALSE, recursive = TRUE)
dir.create("05_figs",         showWarnings = FALSE, recursive = TRUE)
dir.create("06_connectivity", showWarnings = FALSE, recursive = TRUE)
dir.create("07_results",      showWarnings = FALSE, recursive = TRUE)

log_file <- file.path("07_results", "pipeline_log.txt")
if (overwrite && file.exists(log_file)) file.remove(log_file)

out_occ_clean     <- file.path("02_occ", "Otter_TX_Cleaned_Geocoded.csv")
out_county_pres   <- file.path("07_results", "Otter_TX_CountyPresence.csv")
out_period_sum    <- file.path("07_results", "Otter_TX_PeriodSummary_records_by_county.csv")
out_newcounties   <- file.path("07_results", "Otter_TX_NewCounties_by_FirstDetectionPeriod.csv")

out_pres_thin     <- file.path("02_occ", sprintf("pres_thin_%dkm.csv", thin_km))
out_pres_env      <- file.path("02_occ", "Otter_TX_Presence_Env.csv")
out_bg_env        <- file.path("02_occ", "Otter_TX_Background_Env.csv")

out_dens_nat      <- file.path("01_intermediate", "dens_nat_km2.tif")
out_dens_hum      <- file.path("01_intermediate", "dens_hum_km2.tif")
out_dens_all      <- file.path("01_intermediate", "dens_all_km2.tif")
out_dist_riv      <- file.path("01_intermediate", "dist_major_rivers_km.tif")

out_tune          <- file.path("03_models", "Otter_Maxnet_kfold_tuning_results.csv")
out_fold_metrics  <- file.path("03_models", "Otter_Maxnet_kfold_fold_metrics.csv")
out_model_rds     <- file.path("03_models", "Otter_Maxnet_NHD_natHum.rds")
out_thr_txt       <- file.path("03_models", "Otter_Maxnet_NHD_omit10.txt")

out_curr_clog     <- file.path("04_pred", "Otter_Maxnet_NHD_current_cloglog.tif")
out_nohum_clog    <- file.path("04_pred", "Otter_Maxnet_NHD_noHuman_cloglog.tif")

out_beta_csv      <- file.path("03_models", "Otter_Maxnet_beta_raw_features.csv")
out_imp_csv       <- file.path("03_models", "Otter_Maxnet_predictor_importance.csv")
out_imp_png       <- file.path("05_figs",   "Otter_Maxnet_predictor_importance.png")
out_beta_png      <- file.path("05_figs",   "Otter_Maxnet_beta_features.png")
out_resp_png      <- file.path("05_figs",   "Otter_Maxnet_response_curves_allWater.png")

out_sessinfo_txt  <- file.path("07_results","sessionInfo_full_pipeline.txt")

# --------------------------- 3. PACKAGES -------------------------------
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(lubridate)
  library(sf)
  library(terra)
  library(tigris)
  library(maxnet)
  library(ggplot2)
  library(tidyr)
})
options(tigris_use_cache = TRUE)
sf::sf_use_s2(FALSE)
terra::terraOptions(progress = 1, memfrac = 0.2, todisk = TRUE)
set.seed(seed_value)

# --------------------------- 4. HELPERS --------------------------
msg <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(line, "\n")
  write(line, file = log_file, append = TRUE)
}

stop_missing <- function(paths) {
  miss <- paths[!file.exists(paths)]
  if (length(miss)) stop("Missing file(s):\n", paste0(" - ", miss, collapse = "\n"))
}

print_outputs <- function(paths, prefix = "Outputs") {
  paths <- unique(paths)
  msg(prefix, ":")
  for (p in paths) {
    msg("  - ", p, ifelse(file.exists(p), "", "  [MISSING]"))
  }
}

# "Run or skip" BUT ALWAYS report from disk if outputs exist.
maybe_run <- function(step_label, out_paths, code_block, report_block = NULL) {
  out_paths <- unique(out_paths)
  
  should_run <- overwrite || any(!file.exists(out_paths))
  
  if (should_run) {
    msg(step_label, " | RUNNING (overwrite=", overwrite, ")")
    force(code_block)
  } else {
    msg(step_label, " | SKIPPING (overwrite=", overwrite, "); outputs already exist.")
  }
  
  # Always print where outputs are supposed to be
  print_outputs(out_paths, prefix = paste0(step_label, " | Output locations"))
  
  # Always report results if possible (read from disk)
  if (!is.null(report_block)) {
    if (all(file.exists(out_paths))) {
      msg(step_label, " | Reporting results (from disk)")
      force(report_block)
    } else {
      msg(step_label, " | Cannot report results (some outputs missing).")
    }
  }
  
  invisible(NULL)
}

# AUC from ranks (no extra packages)
auc_fast <- function(pos_scores, neg_scores) {
  x <- c(pos_scores, neg_scores)
  r <- rank(x, ties.method = "average")
  n1 <- length(pos_scores)
  n0 <- length(neg_scores)
  r1 <- sum(r[seq_len(n1)])
  (r1 - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# 10% training omission evaluated on test presences
omit10_rate <- function(train_pos_scores, test_pos_scores) {
  thr <- as.numeric(stats::quantile(train_pos_scores, probs = 0.10, na.rm = TRUE, names = FALSE))
  mean(test_pos_scores < thr, na.rm = TRUE)
}

# Spatial thinning in km (greedy) using sf distances in crs_area
thin_points_km <- function(pts_vect, thin_km, crs_project = crs_area) {
  if (nrow(pts_vect) <= 1) return(pts_vect)
  
  pts_sf <- sf::st_as_sf(pts_vect)
  pts_sf <- sf::st_transform(pts_sf, crs_project)
  
  keep <- logical(nrow(pts_sf))
  keep[1] <- TRUE
  kept_idx <- 1
  
  for (i in 2:nrow(pts_sf)) {
    d <- as.numeric(sf::st_distance(pts_sf[i, ], pts_sf[kept_idx, ], by_element = FALSE))
    if (min(d, na.rm = TRUE) >= thin_km * 1000) {
      keep[i] <- TRUE
      kept_idx <- which(keep)
    }
  }
  
  pts_keep <- pts_sf[keep, ]
  pts_keep <- sf::st_transform(pts_keep, sf::st_crs(sf::st_as_sf(pts_vect)))
  terra::vect(pts_keep)
}

# Read env tables, enforce shared columns/order, numeric + complete cases
read_env_numeric <- function(pres_csv, bg_csv) {
  pres <- readr::read_csv(pres_csv, show_col_types = FALSE)
  bg   <- readr::read_csv(bg_csv,   show_col_types = FALSE)
  
  pred_cols <- intersect(names(pres), names(bg))
  if (length(pred_cols) < 2) {
    stop("Presence/background tables do not share enough columns.\n",
         "Presence: ", paste(names(pres), collapse = ", "), "\n",
         "Background: ", paste(names(bg), collapse = ", "))
  }
  
  pred_cols <- names(bg)[names(bg) %in% pred_cols]  # canonical order
  pres <- pres[, pred_cols, drop = FALSE]
  bg   <- bg[,   pred_cols, drop = FALSE]
  
  pres_num <- as.data.frame(lapply(pres, function(z) suppressWarnings(as.numeric(z))))
  bg_num   <- as.data.frame(lapply(bg,   function(z) suppressWarnings(as.numeric(z))))
  
  pres_num <- pres_num[stats::complete.cases(pres_num), , drop = FALSE]
  bg_num   <- bg_num[stats::complete.cases(bg_num), , drop = FALSE]
  
  list(pres = pres_num, bg = bg_num, pred_cols = pred_cols)
}

# --------------------------- 5. STEP A: GBIF → CLEAN OCCURRENCES --------
maybe_run(
  step_label = "STEP A",
  out_paths  = c(out_occ_clean),
  code_block = {
    stop_missing(c(gbif_file))
    msg("STEP A | Reading GBIF: ", gbif_file)
    msg("STEP A | Species=", species_name, " years=", year_min, "-", year_max)
    
    if (grepl("\\.xlsx$", gbif_file, ignore.case = TRUE)) {
      suppressPackageStartupMessages(library(readxl))
      gbif <- readxl::read_excel(gbif_file)
    } else if (grepl("\\.csv$", gbif_file, ignore.case = TRUE)) {
      gbif <- readr::read_csv(gbif_file, show_col_types = FALSE)
    } else {
      stop("GBIF file must be .xlsx or .csv: ", gbif_file)
    }
    
    msg("STEP A | Raw GBIF rows: ", nrow(gbif))
    
    nms <- names(gbif); nms_low <- tolower(nms)
    get_col <- function(candidates) {
      idx <- which(nms_low %in% tolower(candidates))
      if (!length(idx)) return(NA_character_)
      nms[idx[1]]
    }
    
    col_sci  <- get_col(c("scientificname", "species", "taxon"))
    col_date <- get_col(c("eventdate", "date", "event_date"))
    col_year <- get_col(c("year", "eventyear", "observationyear"))
    col_lat  <- get_col(c("decimallatitude", "latitude", "lat"))
    col_lon  <- get_col(c("decimallongitude", "longitude", "lon", "long"))
    col_id   <- get_col(c("gbifid", "key", "occurrenceid"))
    
    if (is.na(col_lat) || is.na(col_lon)) stop("Could not find lat/lon columns in GBIF export.")
    
    occ <- gbif %>%
      transmute(
        gbifID           = if (!is.na(col_id)) .data[[col_id]] else NA,
        scientificName   = if (!is.na(col_sci)) .data[[col_sci]] else NA,
        eventDate_raw    = if (!is.na(col_date)) .data[[col_date]] else NA,
        year_raw         = if (!is.na(col_year)) .data[[col_year]] else NA,
        decimalLatitude  = as.numeric(.data[[col_lat]]),
        decimalLongitude = as.numeric(.data[[col_lon]])
      )
    
    if (!all(is.na(occ$scientificName))) {
      occ <- occ %>% filter(scientificName == species_name)
    }
    
    occ <- occ %>%
      mutate(
        year = suppressWarnings(as.numeric(year_raw)),
        eventDate = suppressWarnings(as.Date(eventDate_raw)),
        year = if_else(is.na(year) & !is.na(eventDate), lubridate::year(eventDate), year)
      ) %>%
      select(-eventDate_raw, -year_raw)
    
    occ <- occ %>%
      filter(is.finite(decimalLongitude), is.finite(decimalLatitude)) %>%
      filter(decimalLatitude >= -90, decimalLatitude <= 90) %>%
      filter(decimalLongitude >= -180, decimalLongitude <= 180) %>%
      filter(!is.na(year), year >= year_min, year <= year_max) %>%
      distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE)
    
    write_csv(occ, out_occ_clean)
  },
  report_block = {
    occ <- readr::read_csv(out_occ_clean, show_col_types = FALSE)
    msg("STEP A | Cleaned independent georeferenced records (n): ", nrow(occ),
        "  [manuscript cross-check]")
  }
)

# --------------------------- 6. STEP B: COUNTY + PERIOD SUMS ------------
maybe_run(
  step_label = "STEP B",
  out_paths  = c(out_county_pres, out_period_sum, out_newcounties),
  code_block = {
    occ <- readr::read_csv(out_occ_clean, show_col_types = FALSE)
    
    if (!("year" %in% names(occ))) {
      ycol <- names(occ)[tolower(names(occ)) %in% c("year","eventyear","observationyear")]
      if (length(ycol)) occ$year <- occ[[ycol[1]]]
    }
    occ$year <- suppressWarnings(as.numeric(occ$year))
    
    stopifnot(all(c("decimalLongitude","decimalLatitude") %in% names(occ)))
    occ_sf <- sf::st_as_sf(occ, coords = c("decimalLongitude", "decimalLatitude"),
                           crs = 4326, remove = FALSE)
    
    occ_sf <- occ_sf %>%
      mutate(period = case_when(
        !is.na(year) & year < 1985 ~ "Pre-1985",
        !is.na(year) & year >= 1985 & year < 1995 ~ "1985-1995",
        !is.na(year) & year >= 1995 & year < 2005 ~ "1995-2005",
        !is.na(year) & year >= 2005 & year < 2015 ~ "2005-2015",
        !is.na(year) & year >= 2015 & year <= 2025 ~ "2015-2025",
        TRUE ~ NA_character_
      ))
    
    tx_counties <- tigris::counties(state = region_state, cb = TRUE, year = 2022)
    
    alb <- 5070
    occ_alb <- st_transform(occ_sf, alb)
    txc_alb <- st_transform(tx_counties, alb)
    
    occ_with_cty <- st_join(st_buffer(occ_alb, dist = 100), txc_alb,
                            join = st_intersects, left = TRUE)
    
    nm <- names(occ_with_cty)
    county_col <- dplyr::case_when(
      "NAME"   %in% nm ~ "NAME",
      "NAME.y" %in% nm ~ "NAME.y",
      TRUE ~ NA_character_
    )
    if (is.na(county_col)) stop("County name column not found after join. Columns were:\n", paste(nm, collapse=", "))
    
    occ_tab <- st_drop_geometry(occ_with_cty)
    
    county_summary <- occ_tab %>%
      filter(!is.na(.data[[county_col]]), !is.na(year)) %>%
      group_by(county = .data[[county_col]]) %>%
      summarise(
        n_records    = n(),
        first_record = min(year, na.rm = TRUE),
        last_record  = max(year, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(n_records))
    
    county_presence <- county_summary %>%
      transmute(county, presence = 1, first_record, last_record, n_records)
    
    period_summary <- occ_tab %>%
      filter(!is.na(.data[[county_col]]), !is.na(period)) %>%
      group_by(period, county = .data[[county_col]]) %>%
      summarise(n_records = n(), .groups = "drop") %>%
      arrange(period, desc(n_records))
    
    first_by_county <- county_summary %>%
      transmute(county, first_year = first_record) %>%
      mutate(period = case_when(
        first_year < 1985 ~ "Pre-1985",
        first_year >= 1985 & first_year < 1995 ~ "1985-1995",
        first_year >= 1995 & first_year < 2005 ~ "1995-2005",
        first_year >= 2005 & first_year < 2015 ~ "2005-2015",
        first_year >= 2015 & first_year <= 2025 ~ "2015-2025",
        TRUE ~ NA_character_
      )) %>%
      count(period, name = "n_new_counties") %>%
      mutate(period = factor(period, levels = c("Pre-1985","1985-1995","1995-2005","2005-2015","2015-2025"))) %>%
      arrange(period)
    
    write_csv(county_presence, out_county_pres)
    write_csv(period_summary,  out_period_sum)
    write_csv(first_by_county, out_newcounties)
  },
  report_block = {
    county_presence <- readr::read_csv(out_county_pres, show_col_types = FALSE)
    first_by_county <- readr::read_csv(out_newcounties, show_col_types = FALSE)
    msg("STEP B | Counties with ≥1 record (n): ", nrow(county_presence), "  [manuscript cross-check]")
    msg("STEP B | New counties by first-detection period (table):")
    print(first_by_county)
  }
)

# --------------------------- 7. STEP C: BUILD PREDICTORS ----------------
maybe_run(
  step_label = "STEP C1",
  out_paths  = c(out_dens_nat, out_dens_hum, out_dens_all),
  code_block = {
    stop_missing(c(natwat_focal, humwat_focal, allwat_focal))
    msg("STEP C1 | Reading focal rasters and building density predictors")
    
    r_nat <- rast(natwat_focal)
    r_hum <- rast(humwat_focal)
    r_all <- rast(allwat_focal)
    
    if (!compareGeom(r_nat, r_hum, stopOnError = FALSE)) r_hum <- resample(r_hum, r_nat, method = "bilinear")
    if (!compareGeom(r_nat, r_all, stopOnError = FALSE)) r_all <- resample(r_all, r_nat, method = "bilinear")
    
    radius_km <- 2.93
    window_area_km2 <- pi * radius_km^2
    
    dens_nat_km2 <- r_nat * window_area_km2
    dens_hum_km2 <- r_hum * window_area_km2
    dens_all_km2 <- r_all * window_area_km2
    
    names(dens_nat_km2) <- "dens_nat_km2"
    names(dens_hum_km2) <- "dens_hum_km2"
    names(dens_all_km2) <- "dens_all_km2"
    
    writeRaster(dens_nat_km2, out_dens_nat, overwrite = TRUE)
    writeRaster(dens_hum_km2, out_dens_hum, overwrite = TRUE)
    writeRaster(dens_all_km2, out_dens_all, overwrite = TRUE)
  },
  report_block = {
    msg("STEP C1 | Density rasters ready for extraction/prediction.")
  }
)

maybe_run(
  step_label = "STEP C2",
  out_paths  = c(out_dist_riv),
  code_block = {
    stop_missing(c(out_dens_nat, major_rivers_shp))
    msg("STEP C2 | Computing distance-to-major-river raster (km)")
    
    r_template <- rast(out_dens_nat)
    rivers <- vect(major_rivers_shp)
    if (crs(rivers) != crs(r_template)) rivers <- project(rivers, crs(r_template))
    
    dist_m  <- distance(r_template, rivers)
    dist_km <- dist_m / 1000
    names(dist_km) <- "dist_main_river_km"
    
    writeRaster(dist_km, out_dist_riv, overwrite = TRUE)
  },
  report_block = {
    msg("STEP C2 | Distance raster ready for extraction/prediction.")
  }
)

# --------------------------- 8. STEP D: PRESENCE + THINNING + EXTRACT ----
maybe_run(
  step_label = "STEP D",
  out_paths  = c(out_pres_thin, out_pres_env),
  code_block = {
    occ <- readr::read_csv(out_occ_clean, show_col_types = FALSE)
    
    dens_nat_km2 <- terra::rast(out_dens_nat)
    dens_hum_km2 <- terra::rast(out_dens_hum)
    dist_main    <- terra::rast(out_dist_riv)
    
    stack_curr <- c(dist_main, dens_nat_km2, dens_hum_km2)
    names(stack_curr) <- c("dist_main_river_km", "dens_nat_km2", "dens_hum_km2")
    
    pts <- terra::vect(occ, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
    pts <- terra::project(pts, terra::crs(dens_nat_km2))
    
    mask_nonNA <- !is.na(dens_nat_km2)
    in_mask <- !is.na(terra::extract(mask_nonNA, pts)[, 1])
    pts <- pts[in_mask, ]
    
    n_before <- nrow(pts)
    if (n_before < 20) stop("Too few points after raster footprint filter: ", n_before)
    
    pts_thin <- thin_points_km(pts, thin_km = thin_km, crs_project = crs_area)
    n_after  <- nrow(pts_thin)
    
    pts_thin_ll <- terra::project(pts_thin, "EPSG:4326")
    xy <- as.data.frame(terra::crds(pts_thin_ll))
    names(xy) <- c("decimalLongitude", "decimalLatitude")
    readr::write_csv(xy, out_pres_thin)
    
    pres_vals <- terra::extract(stack_curr, pts_thin) |> as.data.frame()
    if ("ID" %in% names(pres_vals)) pres_vals <- dplyr::select(pres_vals, -ID)
    pres_vals <- pres_vals[, names(stack_curr), drop = FALSE]
    pres_vals <- pres_vals[stats::complete.cases(pres_vals), , drop = FALSE]
    if (nrow(pres_vals) == 0) stop("Presence extraction returned 0 complete cases.")
    
    readr::write_csv(pres_vals, out_pres_env)
    
    msg("STEP D | Presence points before thinning: ", n_before)
    msg("STEP D | Presence points after thinning (", thin_km, " km): ", n_after,
        "  [manuscript cross-check]")
  },
  report_block = {
    pres_xy  <- readr::read_csv(out_pres_thin, show_col_types = FALSE)
    pres_env <- readr::read_csv(out_pres_env,  show_col_types = FALSE)
    msg("STEP D | Thinned presence points (n): ", nrow(pres_xy), "  [manuscript cross-check]")
    msg("STEP D | Presence env table rows (n complete cases): ", nrow(pres_env))
  }
)

# --------------------------- 9. STEP E: BACKGROUND SAMPLE --------------
maybe_run(
  step_label = "STEP E",
  out_paths  = c(out_bg_env),
  code_block = {
    dens_nat_km2 <- terra::rast(out_dens_nat)
    dens_hum_km2 <- terra::rast(out_dens_hum)
    dist_main    <- terra::rast(out_dist_riv)
    
    stack_curr <- c(dist_main, dens_nat_km2, dens_hum_km2)
    names(stack_curr) <- c("dist_main_river_km", "dens_nat_km2", "dens_hum_km2")
    
    mask_nonNA <- !is.na(dens_nat_km2)
    
    set.seed(seed_value)
    bg_pts <- terra::spatSample(
      x         = mask_nonNA,
      size      = n_background,
      method    = "random",
      as.points = TRUE,
      values    = FALSE,
      na.rm     = TRUE
    )
    
    bg_vals <- terra::extract(stack_curr, bg_pts) |> as.data.frame()
    if ("ID" %in% names(bg_vals)) bg_vals <- dplyr::select(bg_vals, -ID)
    
    bg_vals <- bg_vals[, names(stack_curr), drop = FALSE]
    bg_vals <- bg_vals[stats::complete.cases(bg_vals), , drop = FALSE]
    if (nrow(bg_vals) == 0) stop("Background extraction returned 0 complete cases.")
    
    readr::write_csv(bg_vals, out_bg_env)
  },
  report_block = {
    bg_env <- readr::read_csv(out_bg_env, show_col_types = FALSE)
    msg("STEP E | Background samples extracted (n complete cases): ", nrow(bg_env))
  }
)

# --------------------------- 10. STEP F: TUNING + K-FOLD CV -------------
maybe_run(
  step_label = "STEP F",
  out_paths  = c(out_tune, out_fold_metrics),
  code_block = {
    env <- read_env_numeric(out_pres_env, out_bg_env)
    pres_vals_num <- env$pres
    bg_vals_num   <- env$bg
    
    if (nrow(pres_vals_num) < k_folds) stop("Too few presence rows after NA drop.")
    if (nrow(bg_vals_num) < 100) stop("Too few background rows after NA drop.")
    
    set.seed(seed_value)
    fold_id <- sample(rep(1:k_folds, length.out = nrow(pres_vals_num)))
    
    tune_grid <- tidyr::crossing(fc = fc_grid, rm = rm_grid) %>%
      mutate(model_id = sprintf("fc_%s_rm_%s", fc, rm))
    
    fold_out <- list()
    tune_out <- list()
    
    for (i in seq_len(nrow(tune_grid))) {
      fc  <- tolower(tune_grid$fc[i])
      rm  <- as.numeric(tune_grid$rm[i])
      mid <- tune_grid$model_id[i]
      
      msg("STEP F | CV tuning: ", mid)
      
      fold_list <- vector("list", k_folds)
      
      for (fold in 1:k_folds) {
        test_idx  <- which(fold_id == fold)
        train_idx <- which(fold_id != fold)
        
        x_train <- rbind(pres_vals_num[train_idx, , drop = FALSE], bg_vals_num)
        p_train <- c(rep(1, length(train_idx)), rep(0, nrow(bg_vals_num)))
        
        x_test  <- rbind(pres_vals_num[test_idx, , drop = FALSE], bg_vals_num)
        
        f <- maxnet::maxnet.formula(p_train, x_train, classes = fc)
        m <- maxnet::maxnet(p = p_train, data = x_train, f = f, regmult = rm)
        
        s_train <- as.numeric(predict(m, x_train, type = "cloglog"))
        s_test  <- as.numeric(predict(m, x_test,  type = "cloglog"))
        
        s_test_pos  <- s_test[seq_len(length(test_idx))]
        s_test_bg   <- s_test[(length(test_idx) + 1):length(s_test)]
        s_train_pos <- s_train[seq_len(length(train_idx))]
        
        fold_list[[fold]] <- data.frame(
          model_id = mid, fc = fc, rm = rm, fold = fold,
          n_train_pres = length(train_idx),
          n_test_pres  = length(test_idx),
          auc_test     = auc_fast(s_test_pos, s_test_bg),
          omit10_test  = omit10_rate(s_train_pos, s_test_pos)
        )
      }
      
      fold_df <- bind_rows(fold_list)
      fold_out[[mid]] <- fold_df
      
      tune_out[[mid]] <- fold_df %>%
        summarise(
          model_id = first(model_id),
          fc = first(fc),
          rm = first(rm),
          mean_auc_test = mean(auc_test, na.rm = TRUE),
          sd_auc_test   = sd(auc_test, na.rm = TRUE),
          mean_omit10   = mean(omit10_test, na.rm = TRUE),
          sd_omit10     = sd(omit10_test, na.rm = TRUE),
          .groups = "drop"
        )
    }
    
    fold_metrics <- bind_rows(fold_out)
    tune_metrics <- bind_rows(tune_out) %>%
      arrange(desc(mean_auc_test), mean_omit10, rm)
    
    write_csv(fold_metrics, out_fold_metrics)
    write_csv(tune_metrics, out_tune)
  },
  report_block = {
    tune_metrics <- readr::read_csv(out_tune, show_col_types = FALSE)
    best <- tune_metrics[1, , drop = FALSE]
    msg("STEP F | Best CV model (from disk): fc=", best$fc[1], " rm=", best$rm[1],
        " | mean AUC=", round(best$mean_auc_test[1], 4), " ± ", round(best$sd_auc_test[1], 4),
        " | mean omit10=", round(best$mean_omit10[1], 4), " ± ", round(best$sd_omit10[1], 4),
        "  [manuscript cross-check]")
  }
)

# --------------------------- 11. STEP G: FINAL MODEL + THRESHOLD --------
maybe_run(
  step_label = "STEP G",
  out_paths  = c(out_model_rds, out_thr_txt),
  code_block = {
    env <- read_env_numeric(out_pres_env, out_bg_env)
    pres_vals_num <- env$pres
    bg_vals_num   <- env$bg
    
    tune_metrics <- read_csv(out_tune, show_col_types = FALSE)
    best_row <- tune_metrics[1, , drop = FALSE]
    
    best_fc <- tolower(as.character(best_row$fc[1]))
    best_rm <- as.numeric(best_row$rm[1])
    
    msg("STEP G | Selected final model: fc=", best_fc, " rm=", best_rm)
    
    x_all <- rbind(pres_vals_num, bg_vals_num)
    p_all <- c(rep(1L, nrow(pres_vals_num)), rep(0L, nrow(bg_vals_num)))
    
    f_all <- maxnet::maxnet.formula(p_all, x_all, classes = best_fc)
    mod_maxnet <- maxnet::maxnet(p = p_all, data = x_all, f = f_all, regmult = best_rm)
    
    saveRDS(mod_maxnet, out_model_rds)
    
    pred_pres <- as.numeric(predict(mod_maxnet, pres_vals_num, type = "cloglog"))
    thr_omit10 <- as.numeric(quantile(pred_pres, probs = 0.10, na.rm = TRUE, names = FALSE))
    
    writeLines(sprintf("omit10_cloglog_threshold: %.6f", thr_omit10), con = out_thr_txt)
  },
  report_block = {
    thr_line <- readLines(out_thr_txt)[1]
    thr <- as.numeric(sub(".*: *", "", thr_line))
    msg("STEP G | Omit-10% threshold (cloglog): ", sprintf("%.6f", thr),
        "  [manuscript cross-check]")
  }
)

# --------------------------- 12. STEP H: EXPORT SCENARIO PREDICTIONS -----
maybe_run(
  step_label = "STEP H",
  out_paths  = c(out_curr_clog, out_nohum_clog),
  code_block = {
    # Stable terra tempdir to prevent missing-temp-file errors on Windows
    terra_tmp <- file.path(getwd(), "_terra_tmp")
    dir.create(terra_tmp, showWarnings = FALSE, recursive = TRUE)
    terra::terraOptions(tempdir = terra_tmp, todisk = TRUE, progress = 1)
    
    mod_maxnet <- readRDS(out_model_rds)
    
    thr_line   <- readLines(out_thr_txt)[1]
    thr_omit10 <- as.numeric(sub(".*: *", "", thr_line))
    stopifnot(is.finite(thr_omit10))
    
    dens_nat_km2 <- terra::rast(out_dens_nat)
    dens_hum_km2 <- terra::rast(out_dens_hum)
    dist_main    <- terra::rast(out_dist_riv)
    
    names(dens_nat_km2) <- "dens_nat_km2"
    names(dens_hum_km2) <- "dens_hum_km2"
    names(dist_main)    <- "dist_main_river_km"
    
    stack_curr <- c(dist_main, dens_nat_km2, dens_hum_km2)
    pred_names <- names(stack_curr)
    
    dens_hum_zero <- dens_hum_km2 * 0
    names(dens_hum_zero) <- "dens_hum_km2"
    
    stack_nohum <- c(dist_main, dens_nat_km2, dens_hum_zero)
    names(stack_nohum) <- pred_names
    
    # Avoid overwrite locks: if target exists, write a timestamped copy
    safe_outfile <- function(path) {
      if (!file.exists(path)) return(path)
      stem <- tools::file_path_sans_ext(path)
      ext  <- tools::file_ext(path)
      ts   <- format(Sys.time(), "%Y%m%d_%H%M%S")
      paste0(stem, "_", ts, ".", ext)
    }
    out_curr_use  <- safe_outfile(out_curr_clog)
    out_nohum_use <- safe_outfile(out_nohum_clog)
    
    msg("STEP H | Exporting scenario rasters (cloglog). Threshold for GIS: ", sprintf("%.6f", thr_omit10))
    
    terra::predict(
      object    = stack_curr,
      model     = mod_maxnet,
      type      = "cloglog",
      na.rm     = TRUE,
      filename  = out_curr_use,
      overwrite = TRUE
    )
    
    terra::predict(
      object    = stack_nohum,
      model     = mod_maxnet,
      type      = "cloglog",
      na.rm     = TRUE,
      filename  = out_nohum_use,
      overwrite = TRUE
    )
    
    # If we wrote to timestamped outputs, update the canonical path *pointer* file for humans
    writeLines(out_curr_use,  file.path("04_pred", "LATEST_current_cloglog_path.txt"))
    writeLines(out_nohum_use, file.path("04_pred", "LATEST_noHuman_cloglog_path.txt"))
    
    msg("STEP H | Wrote current: ", out_curr_use)
    msg("STEP H | Wrote noHuman:  ", out_nohum_use)
    msg("STEP H | Also wrote pointers: 04_pred/LATEST_*_path.txt")
  },
  report_block = {
    thr_line <- readLines(out_thr_txt)[1]
    thr <- as.numeric(sub(".*: *", "", thr_line))
    msg("STEP H | Reminder: threshold for GIS binary maps (omit-10% cloglog) = ", sprintf("%.6f", thr))
    
    # Prefer the pointer files (because Step H may timestamp outputs)
    p1 <- file.path("04_pred", "LATEST_current_cloglog_path.txt")
    p2 <- file.path("04_pred", "LATEST_noHuman_cloglog_path.txt")
    if (file.exists(p1)) msg("STEP H | Current cloglog raster: ", readLines(p1)[1])
    if (file.exists(p2)) msg("STEP H | NoHuman cloglog raster: ",  readLines(p2)[1])
  }
)

# --------------------------- 13. STEP I: MODEL INTERPRETATION OUTPUTS ----
maybe_run(
  step_label = "STEP I",
  out_paths  = c(out_beta_csv, out_imp_csv, out_imp_png, out_beta_png, out_resp_png),
  code_block = {
    mod_maxnet <- readRDS(out_model_rds)
    
    beta_vec <- mod_maxnet$betas
    beta_df <- data.frame(
      feature = names(beta_vec),
      beta    = as.numeric(beta_vec),
      stringsAsFactors = FALSE
    ) %>%
      filter(feature != "(Intercept)") %>%
      mutate(abs_beta = abs(beta)) %>%
      arrange(desc(abs_beta))
    
    write_csv(beta_df, out_beta_csv)
    
    beta_grouped <- beta_df %>%
      mutate(
        predictor = case_when(
          grepl("dist_main_river_km", feature) ~ "Distance to major river (km)",
          grepl("dens_nat_km2",       feature) ~ "Natural water density (km²)",
          grepl("dens_hum_km2",       feature) ~ "Anthropogenic water density (km²)",
          TRUE ~ "Other"
        )
      )
    
    var_imp <- beta_grouped %>%
      group_by(predictor) %>%
      summarise(sum_abs_beta = sum(abs_beta), .groups = "drop") %>%
      arrange(desc(sum_abs_beta))
    
    write_csv(var_imp, out_imp_csv)
    
    p_imp <- ggplot(var_imp, aes(x = reorder(predictor, sum_abs_beta), y = sum_abs_beta)) +
      geom_col() +
      coord_flip() +
      labs(x = NULL, y = "Σ |β| across features", title = "Predictor importance – otter Maxnet model") +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(hjust = 0.5), panel.grid.major.y = element_blank())
    
    ggsave(out_imp_png, p_imp, width = 6, height = 4, dpi = 600)
    
    p_feat <- ggplot(beta_df, aes(x = reorder(feature, abs_beta), y = beta)) +
      geom_col() +
      coord_flip() +
      labs(x = NULL, y = "Coefficient (β)", title = "Maxnet feature coefficients – otter SDM") +
      theme_minimal(base_size = 10) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.y = element_text(size = 6),
            panel.grid.major.y = element_blank())
    
    ggsave(out_beta_png, p_feat, width = 7, height = 5, dpi = 600)
    
    dens_nat_km2 <- rast(out_dens_nat)
    dens_hum_km2 <- rast(out_dens_hum)
    dist_main    <- rast(out_dist_riv)
    stack_curr <- c(dist_main, dens_nat_km2, dens_hum_km2)
    names(stack_curr) <- c("dist_main_river_km", "dens_nat_km2", "dens_hum_km2")
    
    preds <- c("dist_main_river_km", "dens_nat_km2", "dens_hum_km2")
    
    set.seed(seed_value)
    samp_vec <- terra::spatSample(
      stack_curr,
      size      = 50000,
      method    = "random",
      as.points = TRUE,
      values    = TRUE,
      na.rm     = TRUE
    )
    samp <- as.data.frame(samp_vec)
    meds <- apply(samp[, preds, drop = FALSE], 2, stats::median, na.rm = TRUE)
    
    make_response_curve <- function(pred_name, model, raster_stack, medians, n_points = 200) {
      mm   <- terra::minmax(raster_stack[[pred_name]])
      vmin <- mm[1, 1]
      vmax <- mm[2, 1]
      seq_vals <- seq(vmin, vmax, length.out = n_points)
      
      base_row <- data.frame(
        dist_main_river_km = as.numeric(medians["dist_main_river_km"]),
        dens_nat_km2       = as.numeric(medians["dens_nat_km2"]),
        dens_hum_km2       = as.numeric(medians["dens_hum_km2"])
      )
      
      newdat <- base_row[rep(1, n_points), ]
      newdat[[pred_name]] <- seq_vals
      newdat$pred <- as.numeric(predict(model, newdat, type = "cloglog"))
      
      data.frame(predictor_value = seq_vals, pred = newdat$pred, predictor = pred_name)
    }
    
    resp_df <- bind_rows(lapply(
      preds,
      make_response_curve,
      model        = mod_maxnet,
      raster_stack = stack_curr,
      medians      = meds
    ))
    
    samp$dens_all_km2 <- samp$dens_nat_km2 + samp$dens_hum_km2
    rng_all <- range(samp$dens_all_km2, na.rm = TRUE)
    seq_all <- seq(rng_all[1], rng_all[2], length.out = 200)
    
    ratio_nat <- as.numeric(mean(samp$dens_nat_km2, na.rm = TRUE) / mean(samp$dens_all_km2, na.rm = TRUE))
    ratio_hum <- 1 - ratio_nat
    
    newdat_all <- data.frame(
      dist_main_river_km = as.numeric(meds["dist_main_river_km"]),
      dens_nat_km2       = as.numeric(seq_all * ratio_nat),
      dens_hum_km2       = as.numeric(seq_all * ratio_hum)
    )
    rownames(newdat_all) <- NULL
    
    resp_all <- data.frame(
      predictor_value = seq_all,
      pred            = as.numeric(predict(mod_maxnet, newdat_all, type = "cloglog")),
      predictor       = "dens_all_km2"
    )
    
    resp_df <- bind_rows(resp_df, resp_all)
    
    resp_df$predictor <- dplyr::recode(resp_df$predictor,
                                       "dist_main_river_km" = "Distance to Major River (km)",
                                       "dens_nat_km2"       = "Natural Water Density (km²)",
                                       "dens_hum_km2"       = "Anthropogenic Water Density (km²)",
                                       "dens_all_km2"       = "Water Density (km²)")
    resp_df$predictor <- factor(resp_df$predictor,
                                levels = c("Water Density (km²)",
                                           "Anthropogenic Water Density (km²)",
                                           "Natural Water Density (km²)",
                                           "Distance to Major River (km)"))
    
    p_resp <- ggplot(resp_df, aes(x = predictor_value, y = pred)) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~ predictor, scales = "free_x", ncol = 1) +
      labs(x = "Predictor Value", y = "Predicted Suitability (cloglog)",
           title = "Maxent response curves – river otter SDM") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5),
            strip.text = element_text(size = 11, face = "bold"))
    
    ggsave(out_resp_png, p_resp, width = 7, height = 10, dpi = 600)
  },
  report_block = {
    var_imp <- readr::read_csv(out_imp_csv, show_col_types = FALSE)
    msg("STEP I | Predictor importance Σ|β| (from disk; manuscript cross-check):")
    print(var_imp)
    msg("STEP I | Figures written: ", out_imp_png, " ; ", out_beta_png, " ; ", out_resp_png)
  }
)
 
# --------------------------- 14. STEP K: COST / ACCESSIBILITY SUMMARIES ----
# Expected CSV columns:
#   NatDistAcc, CurrDistAcc, DeltaCost   (DeltaCost can be computed if absent)
cost_csv <- file.path("07_results", "CostAccessibility.csv")  # <-- rename if needed

out_cost_stats_csv  <- file.path("07_results", "CostAccessibility_summary_stats.csv")
out_cost_corr_csv   <- file.path("07_results", "CostAccessibility_correlations.csv")
out_cost_tests_txt  <- file.path("07_results", "CostAccessibility_tests.txt")

out_hist_delta_png  <- file.path("05_figs",   "DeltaCost_histogram.png")
out_dens_delta_png  <- file.path("05_figs",   "DeltaCost_density.png")
out_box_png         <- file.path("05_figs",   "Nat_vs_Curr_boxplot.png")
out_scatter_png     <- file.path("05_figs",   "Nat_vs_Curr_scatter.png")
out_ba_png          <- file.path("05_figs",   "BlandAltman_Curr_minus_Nat.png")

maybe_run(
  step_label = "STEP K",
  out_paths  = c(out_cost_stats_csv, out_cost_corr_csv, out_cost_tests_txt,
                 out_hist_delta_png, out_dens_delta_png, out_box_png, out_scatter_png, out_ba_png),
  code_block = {
    stop_missing(c(cost_csv))
    dat <- readr::read_csv(cost_csv, show_col_types = FALSE)
    
    # Enforce expected names; compute DeltaCost if needed
    req <- c("NatDistAcc", "CurrDistAcc")
    if (!all(req %in% names(dat))) {
      stop("STEP K | Missing required columns. Need: ", paste(req, collapse = ", "),
           "\nFound: ", paste(names(dat), collapse = ", "))
    }
    if (!("DeltaCost" %in% names(dat))) {
      dat <- dat %>% mutate(DeltaCost = CurrDistAcc - NatDistAcc)
      msg("STEP K | DeltaCost column not found; computed as CurrDistAcc - NatDistAcc.")
    }
    
    # Coerce numeric + complete cases for the three columns
    dat_use <- dat %>%
      transmute(
        NatDistAcc  = suppressWarnings(as.numeric(NatDistAcc)),
        CurrDistAcc = suppressWarnings(as.numeric(CurrDistAcc)),
        DeltaCost   = suppressWarnings(as.numeric(DeltaCost))
      ) %>%
      filter(complete.cases(.))
    
    if (nrow(dat_use) < 5) stop("STEP K | Too few complete rows after cleaning: ", nrow(dat_use))
    msg("STEP K | Rows used (complete cases): ", nrow(dat_use))
    
    # --------- Summary stats (mean/SD + median/IQR + quantiles) ----------
    qfun <- function(x) stats::quantile(x, probs = c(0, .25, .5, .75, 1), na.rm = TRUE, names = FALSE)
    summ <- dat_use %>%
      tidyr::pivot_longer(cols = everything(), names_to = "metric", values_to = "value") %>%
      group_by(metric) %>%
      summarise(
        n      = sum(is.finite(value)),
        mean   = mean(value, na.rm = TRUE),
        sd     = stats::sd(value, na.rm = TRUE),
        min    = qfun(value)[1],
        q1     = qfun(value)[2],
        median = qfun(value)[3],
        q3     = qfun(value)[4],
        max    = qfun(value)[5],
        iqr    = q3 - q1,
        .groups = "drop"
      )
    
    readr::write_csv(summ, out_cost_stats_csv)
    
    # --------- Correlations (Pearson + Spearman) ----------
    cor_df <- tidyr::crossing(
      x = c("NatDistAcc","CurrDistAcc","DeltaCost"),
      y = c("NatDistAcc","CurrDistAcc","DeltaCost")
    ) %>%
      filter(x < y) %>%
      rowwise() %>%
      mutate(
        pearson  = stats::cor(dat_use[[x]], dat_use[[y]], method = "pearson",  use = "complete.obs"),
        spearman = stats::cor(dat_use[[x]], dat_use[[y]], method = "spearman", use = "complete.obs")
      ) %>%
      ungroup()
    
    readr::write_csv(cor_df, out_cost_corr_csv)
    
    # --------- Paired tests (Curr vs Nat) + Delta summary text ----------
    tt <- stats::t.test(dat_use$CurrDistAcc, dat_use$NatDistAcc, paired = TRUE)
    wt <- stats::wilcox.test(dat_use$CurrDistAcc, dat_use$NatDistAcc, paired = TRUE, exact = FALSE)
    
    # Effect size (paired Cohen's d on differences)
    diffs <- dat_use$CurrDistAcc - dat_use$NatDistAcc
    d_paired <- mean(diffs, na.rm = TRUE) / stats::sd(diffs, na.rm = TRUE)
    
    # Write a compact text file for manuscript cross-check
    con <- file(out_cost_tests_txt, open = "wt")
    writeLines(c(
      paste0("Rows used (complete cases): ", nrow(dat_use)),
      "",
      "Paired comparison: CurrDistAcc vs NatDistAcc",
      paste0("Mean difference (Curr - Nat): ", signif(mean(diffs), 6)),
      paste0("SD difference: ", signif(stats::sd(diffs), 6)),
      paste0("Cohen's d (paired, on differences): ", signif(d_paired, 6)),
      "",
      "Paired t-test:",
      paste0("  t = ", signif(unname(tt$statistic), 6),
             ", df = ", unname(tt$parameter),
             ", p = ", signif(tt$p.value, 6)),
      paste0("  95% CI of mean diff: [", signif(tt$conf.int[1], 6), ", ", signif(tt$conf.int[2], 6), "]"),
      "",
      "Wilcoxon signed-rank test:",
      paste0("  V = ", signif(unname(wt$statistic), 6),
             ", p = ", signif(wt$p.value, 6))
    ), con)
    close(con)
    
    # --------- Plots ----------
    # 1) Histogram of DeltaCost
    p_hist <- ggplot(dat_use, aes(x = DeltaCost)) +
      geom_histogram(bins = 40) +
      labs(x = "DeltaCost", y = "Count", title = "Histogram of DeltaCost") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(out_hist_delta_png, p_hist, width = 6.5, height = 4.5, dpi = 600)
    
    # 2) Density of DeltaCost
    p_dens <- ggplot(dat_use, aes(x = DeltaCost)) +
      geom_density(linewidth = 0.8) +
      labs(x = "DeltaCost", y = "Density", title = "Density of DeltaCost") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(out_dens_delta_png, p_dens, width = 6.5, height = 4.5, dpi = 600)
    
    # 3) Boxplots: Nat vs Curr (long format)
    long <- dat_use %>%
      select(NatDistAcc, CurrDistAcc) %>%
      tidyr::pivot_longer(everything(), names_to = "Scenario", values_to = "DistAcc")
    p_box <- ggplot(long, aes(x = Scenario, y = DistAcc)) +
      geom_boxplot(outlier.alpha = 0.5) +
      labs(x = NULL, y = "Distance/Accessibility metric", title = "NatDistAcc vs CurrDistAcc") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(out_box_png, p_box, width = 5.5, height = 4.5, dpi = 600)
    
    # 4) Scatter: Nat vs Curr with 1:1 line
    p_scatter <- ggplot(dat_use, aes(x = NatDistAcc, y = CurrDistAcc)) +
      geom_point(alpha = 0.6) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      labs(x = "NatDistAcc", y = "CurrDistAcc", title = "NatDistAcc vs CurrDistAcc (1:1 line)") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(out_scatter_png, p_scatter, width = 5.5, height = 5.5, dpi = 600)
    
    # 5) Bland–Altman style: diff vs mean (Curr - Nat vs mean)
    ba <- dat_use %>%
      mutate(mean_acc = (CurrDistAcc + NatDistAcc)/2,
             diff_acc = (CurrDistAcc - NatDistAcc))
    md <- mean(ba$diff_acc, na.rm = TRUE)
    sd_d <- stats::sd(ba$diff_acc, na.rm = TRUE)
    loa_u <- md + 1.96*sd_d
    loa_l <- md - 1.96*sd_d
    
    p_ba <- ggplot(ba, aes(x = mean_acc, y = diff_acc)) +
      geom_point(alpha = 0.6) +
      geom_hline(yintercept = md) +
      geom_hline(yintercept = loa_u, linetype = "dashed") +
      geom_hline(yintercept = loa_l, linetype = "dashed") +
      labs(x = "Mean of (Curr, Nat)", y = "Difference (Curr - Nat)",
           title = "Bland–Altman: Curr − Nat vs Mean") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(out_ba_png, p_ba, width = 6.5, height = 4.8, dpi = 600)
    
    msg("STEP K | Wrote summary stats: ", out_cost_stats_csv)
    msg("STEP K | Wrote correlations:  ", out_cost_corr_csv)
    msg("STEP K | Wrote tests:         ", out_cost_tests_txt)
  },
  report_block = {
    summ <- readr::read_csv(out_cost_stats_csv, show_col_types = FALSE)
    msg("STEP K | Summary stats (from disk):")
    print(summ)
    msg("STEP K | Plots written to 05_figs/")
  }
)

# --------------------------- 15. STEP J: SESSION INFO -------------------
maybe_run(
  step_label = "STEP J",
  out_paths  = c(out_sessinfo_txt),
  code_block = {
    sink(out_sessinfo_txt)
    cat("Run time:", as.character(Sys.time()), "\n\n")
    cat("Species:", species_name, "\n")
    cat("Years:", year_min, "-", year_max, "\n")
    cat("Thinning (km):", thin_km, "\n")
    cat("K-folds:", k_folds, "\n")
    cat("Background n:", n_background, "\n\n")
    print(sessionInfo())
    sink()
  },
  report_block = {
    msg("STEP J | Session info written: ", out_sessinfo_txt)
  }
)

msg("PIPELINE COMPLETE.")
msg("Log written to: ", log_file)
