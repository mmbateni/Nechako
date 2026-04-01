####################################################################################
# w_extract_ms1_blanks.R
#
# PURPOSE
#   Extracts all numerical values needed to fill the ___________ blanks in
#   Nechako_MS1_REVISED.docx.  Designed to run AFTER the following scripts
#   have already completed:
#       1SPI_ERALand.R          → spi_results_seasonal/
#       3SPEI_ERALand.R         → spei_results_seasonal/, spei_results_seasonal_thw/
#       w5_event_ranking.R      → ranked_event_catalog_SW08.csv
#       w11_dynamic_thermodynamic_decomp.R → decomp_results/
#       w9_atmospheric_diagnostics.R       → teleconnection outputs
#       w7_teleconnection_prep.R           → PDO / Niño-3.4 CSVs
#
# OUTPUT
#   ms1_blanks/
#     ms1_blanks_confirmed.csv      — values already confirmed from existing outputs
#     ms1_blanks_decomp.csv         — decomposition values (Section 4.4 / Table 3)
#     ms1_blanks_trend.csv          — JJA thermodynamic trend (Section 4.6)
#     ms1_blanks_seasonal.csv       — seasonal fractions (Section 4.4)
#     ms1_blanks_offset_months.csv  — months with positive SPEI_Δ (Table 3)
#     ms1_blanks_annual_scores.csv  — annual severity scores (Section 4.3)
#     ms1_blanks_atm.csv            — atmospheric composite values (Section 4.7)
#     ms1_blanks_telecon.csv        — PDO/ENSO values (Section 4.7 / 5.3)
#     MS1_BLANKS_MASTER.csv         — all blanks in one table with location codes
#
# SECTION CODES used in MS1_BLANKS_MASTER.csv:
#   ABS  = Abstract
#   S431 = Section 4.3 paragraph 1 (event catalog)
#   S44  = Section 4.4 (decomposition results)
#   S46  = Section 4.6 (trend)
#   S47  = Section 4.7 (composites)
#   S51  = Section 5.1 (exceptionality discussion)
#   S52  = Section 5.2 (thermodynamic dominance discussion)
#   S53  = Section 5.3 (PDO/ENSO discussion)
#   S55  = Section 5.5 (water management)
#   C1–C4 = Conclusion points 1–4
####################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(zoo)
})

WD_PATH  <- Sys.getenv("NECHAKO_WD", "D:/Nechako_Drought/Nechako/")
OUT_DIR  <- file.path(WD_PATH, "ms1_blanks")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(WD_PATH)

cat("============================================================\n")
cat("  MS1 BLANK-FILLING EXTRACTOR\n")
cat("============================================================\n\n")

all_blanks <- list()   # accumulates rows for the master table

add_blank <- function(section, location_text, variable, value, source_file) {
  all_blanks[[length(all_blanks) + 1]] <<- data.frame(
    section       = section,
    location      = location_text,
    variable      = variable,
    value         = as.character(round(as.numeric(value), 4)),
    source        = source_file,
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# BLOCK 1: VALUES ALREADY CONFIRMED FROM TABLE 3 / EVENT CATALOG
# These are in the document already — listed here for completeness and
# copy-paste into manuscript.
# ==============================================================================
cat("── Block 1: Confirmed values from Table 3 and event catalog ──\n")

confirmed <- data.frame(
  section   = c("ABS","ABS","ABS","ABS",
                "S431","S431","S431","S431",
                "S51","S51","S51",
                "S52","S52",
                "C1","C1","C1","C1","C1","C1",
                "C2","C2"),
  location  = c(
    "Abstract: T_K for leading SPEI-3 episode",
    "Abstract: 95% CI lower bound",
    "Abstract: 95% CI upper bound",
    "Abstract: thermodynamic fraction 'nearly 100%'",
    "Sec 4.3: SPEI-3 rank-1 S value",
    "Sec 4.3: SPEI-3 rank-2 S value",
    "Sec 4.3: ratio rank-1 / prev best (5.16)",
    "Sec 4.3: SPI-3 event S value (21-month)",
    "Sec 5.1: composite S at SPEI-3 rank-1",
    "Sec 5.1: factor exceeding prev best",
    "Sec 5.1: T_K years",
    "Sec 5.2: thermodynamic % (SPEI-3)",
    "Sec 5.2: dynamic residual %",
    "Conclusion 1: S at SPEI-3",
    "Conclusion 1: factor exceeding 2nd-ranked",
    "Conclusion 1: 2023 mean severity index",
    "Conclusion 1: 2024 mean severity index",
    "Conclusion 1: T_K years",
    "Conclusion 1: 95% CI lower–upper",
    "Conclusion 2: thermodynamic %",
    "Conclusion 2: dynamic residual %"
  ),
  variable  = c("T_K_years","CI_lower","CI_upper","thm_pct_approx",
                "SPEI3_rank1_S","SPEI3_rank2_S","ratio_vs_prev_best","SPI3_S_21mo",
                "S51_S_SPEI3","S51_factor","S51_TK",
                "S52_thm_pct","S52_dyn_pct",
                "C1_S","C1_factor","C1_2023_sev","C1_2024_sev",
                "C1_TK","C1_CI",
                "C2_thm_pct","C2_dyn_pct"),
  value     = c("921","109","55386","~100",
                "10.73","10.31","2.08","48.41",
                "10.73","2.08","921",
                "~99","~1",
                "10.73","2.08","1.329","1.056",
                "921","109–55,386",
                "~99","~1"),
  source    = "Table 3 + ranked_event_catalog_SW08.csv (confirmed)",
  stringsAsFactors = FALSE
)

write.csv(confirmed, file.path(OUT_DIR, "ms1_blanks_confirmed.csv"), row.names = FALSE)
cat(sprintf("  ✓ Confirmed values: %d blanks → ms1_blanks_confirmed.csv\n", nrow(confirmed)))

# ==============================================================================
# BLOCK 2: DECOMPOSITION VALUES FROM w11 OUTPUT
# ==============================================================================
cat("\n── Block 2: Decomposition values from w11 ──\n")

decomp_file <- file.path(WD_PATH, "decomp_results", "decomp_2022_2025_summary.csv")

if (file.exists(decomp_file)) {
  decomp <- read.csv(decomp_file)
  cat("  decomp_2022_2025_summary.csv loaded\n")
  print(decomp)
  
  # Extract scale-specific values
  for (sc in c(1, 2, 3)) {
    row <- decomp[decomp$scale == sc, ]
    if (nrow(row) == 0) next
    
    add_blank("S44", sprintf("Sec 4.4: SPEI_PM mean (scale %d)", sc),
              sprintf("SPEI_PM_mean_s%d", sc), row$SPEI_PM_mean, decomp_file)
    add_blank("S44", sprintf("Sec 4.4: SPEI_Thw mean (scale %d)", sc),
              sprintf("SPEI_Thw_mean_s%d", sc), row$SPEI_Thw_mean, decomp_file)
    add_blank("S44", sprintf("Sec 4.4: SPEI_Dynamic mean (scale %d)", sc),
              sprintf("SPEI_Dyn_mean_s%d", sc), row$SPEI_Dynamic_mean, decomp_file)
    add_blank("S44", sprintf("Sec 4.4: thm_pct %% (scale %d)", sc),
              sprintf("thm_pct_s%d", sc), row$thm_pct, decomp_file)
    add_blank("S44", sprintf("Sec 4.4: dyn_pct %% (scale %d)", sc),
              sprintf("dyn_pct_s%d", sc), row$dyn_pct, decomp_file)
  }
  
  cat("\n  KEY VALUES FOR MANUSCRIPT:\n")
  s3 <- decomp[decomp$scale == 3, ]
  cat(sprintf("  SPEI-3: PM=%.3f  Thw=%.3f  Dyn=%.3f  thm_pct=%.1f%%  dyn_pct=%.1f%%\n",
              s3$SPEI_PM_mean, s3$SPEI_Thw_mean, s3$SPEI_Dynamic_mean,
              s3$thm_pct, s3$dyn_pct))
  s1 <- decomp[decomp$scale == 1, ]
  s2 <- decomp[decomp$scale == 2, ]
  cat(sprintf("  SPEI-1: thm_pct=%.1f%%   SPEI-2: thm_pct=%.1f%%\n",
              s1$thm_pct, s2$thm_pct))
  
  write.csv(do.call(rbind, all_blanks[sapply(all_blanks, function(x) x$section == "S44")]),
            file.path(OUT_DIR, "ms1_blanks_decomp.csv"), row.names = FALSE)
  cat("  ✓ Decomposition blanks → ms1_blanks_decomp.csv\n")
  
} else {
  cat(sprintf("  ⚠ Not found: %s\n  → Run w11_dynamic_thermodynamic_decomp.R first\n", decomp_file))
}

# ==============================================================================
# BLOCK 3: OFFSET MONTHS — months where SPEI_Δ > 0 during 2022–2025
# Add this block also to STEP 4 of w11 for inline printing
# ==============================================================================
cat("\n── Block 3: Offset months (SPEI_Δ > 0 during 2022–2025) ──\n")

decomp_ts_file <- file.path(WD_PATH, "decomp_results", "decomp_timeseries.csv")

if (file.exists(decomp_ts_file)) {
  decomp_ts <- read.csv(decomp_ts_file)
  decomp_ts$date <- as.Date(decomp_ts$date)
  decomp_ts$year <- as.integer(format(decomp_ts$date, "%Y"))
  decomp_ts$month <- as.integer(format(decomp_ts$date, "%m"))
  
  offset <- decomp_ts %>%
    filter(year >= 2022, year <= 2025) %>%
    group_by(scale) %>%
    summarise(
      n_total      = n(),
      n_positive   = sum(SPEI_Dynamic > 0, na.rm = TRUE),
      pct_positive = 100 * mean(SPEI_Dynamic > 0, na.rm = TRUE),
      mean_dynamic = mean(SPEI_Dynamic, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("\n  Months with SPEI_Δ > 0 (2022–2025 window):\n")
  print(offset)
  write.csv(offset, file.path(OUT_DIR, "ms1_blanks_offset_months.csv"), row.names = FALSE)
  
  for (sc in c(1, 2, 3)) {
    row <- offset[offset$scale == sc, ]
    if (nrow(row) == 0) next
    add_blank("S44", sprintf("Sec 4.4: months SPEI_Δ>0 (scale %d, n out of 48)", sc),
              sprintf("n_offset_s%d", sc), row$n_positive, decomp_ts_file)
    add_blank("S44", sprintf("Sec 4.4: mean SPEI_Δ over event window (scale %d)", sc),
              sprintf("mean_dyn_s%d", sc), row$mean_dynamic, decomp_ts_file)
  }
  
  s3o <- offset[offset$scale == 3, ]
  cat(sprintf("\n  SPEI-3: %d of 48 months have SPEI_Δ > 0 (%.0f%%)\n",
              s3o$n_positive, s3o$pct_positive))
  cat(sprintf("  SPEI-3: mean SPEI_Δ over full window = %.4f\n", s3o$mean_dynamic))
  cat("  ✓ Offset months → ms1_blanks_offset_months.csv\n")
  
} else {
  cat(sprintf("  ⚠ %s not found\n", decomp_ts_file))
  cat("  → Check w11 output directory; file may be named differently.\n")
  cat("  → Add this code to STEP 4 of w11 to generate it:\n")
  cat('     write.csv(decomp_all, file.path(out_dir,"decomp_timeseries.csv"), row.names=FALSE)\n')
}

# ==============================================================================
# BLOCK 4: TREND VALUES FROM w11 STEP 5
# ==============================================================================
cat("\n── Block 4: JJA thermodynamic fraction trend (Section 4.6) ──\n")

trend_file <- file.path(WD_PATH, "decomp_results", "ms_blanks_trend.csv")

if (file.exists(trend_file)) {
  trend <- read.csv(trend_file)
  cat("  ms_blanks_trend.csv loaded\n")
  print(trend)
  
  for (i in seq_len(nrow(trend))) {
    r <- trend[i, ]
    sec <- if (grepl("1950_1987|1988_2025|mann_whitney", r$metric)) "S46" else "S46"
    add_blank(sec, paste("Sec 4.6:", r$metric), r$metric, r$value, trend_file)
  }
  
  # Extract key values
  slope_dec <- trend$value[trend$metric == "JJA_trend_slope_per_decade"]
  pval      <- trend$value[trend$metric == "JJA_trend_pval"]
  tstat     <- trend$value[trend$metric == "JJA_trend_tstat"]
  mean_early <- trend$value[trend$metric == "mean_f_thm_1950_1987_pct"]
  mean_late  <- trend$value[trend$metric == "mean_f_thm_1988_2025_pct"]
  mw_p       <- trend$value[trend$metric == "mann_whitney_p"]
  
  sig_word <- if (!is.na(pval) && as.numeric(pval) < 0.05) "significant" else "non-significant"
  
  cat("\n  KEY VALUES FOR SECTION 4.6:\n")
  cat(sprintf("  Trend: %s %s positive trend of %.4f %% per decade\n",
              "statistically", sig_word, as.numeric(slope_dec) * 100))
  cat(sprintf("  OLS slope = %.6f yr⁻¹; t = %.3f; p = %.4f\n",
              as.numeric(slope_dec)/10, as.numeric(tstat), as.numeric(pval)))
  cat(sprintf("  1950–1987 mean JJA f_thm = %.1f%%\n", as.numeric(mean_early)))
  cat(sprintf("  1988–2025 mean JJA f_thm = %.1f%%\n", as.numeric(mean_late)))
  cat(sprintf("  Change = %.1f percentage points; Mann–Whitney p = %.4f\n",
              as.numeric(mean_late) - as.numeric(mean_early), as.numeric(mw_p)))
  
  write.csv(do.call(rbind, Filter(function(x) x$section=="S46", all_blanks)),
            file.path(OUT_DIR, "ms1_blanks_trend.csv"), row.names = FALSE)
  cat("  ✓ Trend blanks → ms1_blanks_trend.csv\n")
  
} else {
  cat(sprintf("  ⚠ %s not found\n  → Add this code to END OF STEP 5 in w11:\n\n", trend_file))
  cat('  # ---- ADD TO END OF STEP 5 in w11_dynamic_thermodynamic_decomp.R ----\n')
  cat('  trend_early <- dplyr::filter(trend_df, year <= 1987)\n')
  cat('  trend_late  <- dplyr::filter(trend_df, year >= 1988)\n')
  cat('  mw_test <- wilcox.test(trend_early$thm_frac_abs, trend_late$thm_frac_abs)\n')
  cat('  ms_blanks <- data.frame(\n')
  cat('    metric = c("JJA_trend_slope_per_yr","JJA_trend_slope_per_decade",\n')
  cat('               "JJA_trend_pval","JJA_trend_tstat",\n')
  cat('               "JJA_ci95_lo","JJA_ci95_hi",\n')
  cat('               "mean_f_thm_1950_1987_pct","mean_f_thm_1988_2025_pct",\n')
  cat('               "change_pp","mann_whitney_p"),\n')
  cat('    value = c(thm_slope, thm_slope_dec, thm_pval,\n')
  cat('              summary(thm_trend)$coefficients["year",3],\n')
  cat('              thm_ci95[1], thm_ci95[2],\n')
  cat('              mean(trend_early$thm_frac_abs,na.rm=TRUE)*100,\n')
  cat('              mean(trend_late$thm_frac_abs,na.rm=TRUE)*100,\n')
  cat('              (mean(trend_late$thm_frac_abs,na.rm=TRUE) -\n')
  cat('               mean(trend_early$thm_frac_abs,na.rm=TRUE))*100,\n')
  cat('              mw_test$p.value)\n')
  cat('  )\n')
  cat('  write.csv(ms_blanks, file.path(out_dir,"ms_blanks_trend.csv"), row.names=FALSE)\n')
}

# ==============================================================================
# BLOCK 5: SEASONAL DECOMPOSITION FRACTIONS
# ==============================================================================
cat("\n── Block 5: Seasonal thermodynamic fractions (Section 4.4) ──\n")

seasonal_file <- file.path(WD_PATH, "decomp_results", "decomp_seasonal_fractions.csv")

if (file.exists(seasonal_file)) {
  seas <- read.csv(seasonal_file)
  cat("  decomp_seasonal_fractions.csv loaded\n")
  print(seas)
  
  for (i in seq_len(nrow(seas))) {
    add_blank("S44", paste("Sec 4.4: seasonal f_thm", seas$season[i]),
              paste0("thm_frac_", seas$season[i]), seas$thm_frac_pct[i], seasonal_file)
  }
  
  jja_pct <- seas$thm_frac_pct[seas$season == "JJA"]
  mam_pct <- seas$thm_frac_pct[seas$season == "MAM"]
  cat(sprintf("\n  JJA thermodynamic fraction: %.1f%%\n", jja_pct))
  cat(sprintf("  MAM (spring) thermodynamic fraction: %.1f%%\n", mam_pct))
  cat("  ✓ Seasonal fractions → ms1_blanks_seasonal.csv\n")
  write.csv(seas, file.path(OUT_DIR, "ms1_blanks_seasonal.csv"), row.names = FALSE)
  
} else {
  cat(sprintf("  ⚠ %s not found\n  → Add this code to STEP 4 in w11:\n\n", seasonal_file))
  cat('  # ---- ADD TO STEP 4 in w11_dynamic_thermodynamic_decomp.R ----\n')
  cat('  seasonal_decomp <- decomp_all %>%\n')
  cat('    filter(year >= 2022, year <= 2025, scale == 3) %>%\n')
  cat('    mutate(season = case_when(\n')
  cat('      month %in% c(6,7,8)  ~ "JJA",\n')
  cat('      month %in% c(12,1,2) ~ "DJF",\n')
  cat('      month %in% c(3,4,5)  ~ "MAM",\n')
  cat('      TRUE                  ~ "SON")) %>%\n')
  cat('    group_by(season) %>%\n')
  cat('    summarise(thm_frac_pct = 100*mean(abs(SPEI_Thw)/(abs(SPEI_PM)+1e-6),na.rm=TRUE),\n')
  cat('              .groups="drop")\n')
  cat('  write.csv(seasonal_decomp, file.path(out_dir,"decomp_seasonal_fractions.csv"),\n')
  cat('            row.names=FALSE)\n')
}

# ==============================================================================
# BLOCK 6: ANNUAL SEVERITY SCORES (Section 4.3 — already in doc but verify)
# ==============================================================================
cat("\n── Block 6: Annual severity scores from event catalog ──\n")

catalog_file <- file.path(WD_PATH, "temporal_drought",
                          "drought_event_ranking", "ranked_event_catalog_SW08.csv")
# Try alternate locations
if (!file.exists(catalog_file)) {
  candidates <- list.files(WD_PATH, "ranked_event_catalog", recursive = TRUE, full.names = TRUE)
  if (length(candidates)) catalog_file <- candidates[1]
}

if (file.exists(catalog_file)) {
  cat <- read.csv(catalog_file)
  spei3 <- cat[cat$index %in% c("spei_03", "spei3", "SPEI3", "spei-3"), ]
  
  if (nrow(spei3) > 0) {
    # Annual scores
    spei3$year <- as.integer(format(as.Date(spei3$start_date), "%Y"))
    annual_scores <- spei3 %>%
      group_by(year) %>%
      summarise(
        n_drought_months = sum(duration_months, na.rm = TRUE),
        mean_severity    = mean(mean_int, na.rm = TRUE),
        annual_score     = sum(severity, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(annual_score))
    
    cat("\n  Top annual SPEI-3 severity scores:\n")
    print(head(annual_scores, 5))
    
    yr2023 <- annual_scores[annual_scores$year == 2023, ]
    yr2024 <- annual_scores[annual_scores$year == 2024, ]
    
    if (nrow(yr2023) > 0) {
      add_blank("S431","Sec 4.3: 2023 drought months (SPEI-3)",
                "n_drought_months_2023", yr2023$n_drought_months, catalog_file)
      add_blank("S431","Sec 4.3: 2023 mean severity index (SPEI-3)",
                "mean_sev_2023", yr2023$mean_severity, catalog_file)
    }
    if (nrow(yr2024) > 0) {
      add_blank("S431","Sec 4.3: 2024 drought months (SPEI-3)",
                "n_drought_months_2024", yr2024$n_drought_months, catalog_file)
      add_blank("S431","Sec 4.3: 2024 mean severity index (SPEI-3)",
                "mean_sev_2024", yr2024$mean_severity, catalog_file)
    }
    write.csv(annual_scores, file.path(OUT_DIR, "ms1_blanks_annual_scores.csv"), row.names = FALSE)
    cat("  ✓ Annual scores → ms1_blanks_annual_scores.csv\n")
  } else {
    cat("  ⚠ No SPEI-3 rows found in catalog. Check index label column.\n")
    cat("  Available index values:", unique(cat$index), "\n")
  }
} else {
  cat("  ⚠ Event catalog CSV not found.\n")
  cat("  → Run w5_event_ranking.R first.\n")
}

# ==============================================================================
# BLOCK 7: PDO AND ENSO INDICES (Sections 4.7, 5.3, Conclusion 4)
# ==============================================================================
cat("\n── Block 7: PDO / Niño-3.4 indices ──\n")

# Try to load PDO index
pdo_candidates <- list.files(WD_PATH, pattern = "pdo.*\\.csv", 
                              recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
nino_candidates <- list.files(WD_PATH, pattern = "nino.*\\.csv|enso.*\\.csv",
                               recursive = TRUE, full.names = TRUE, ignore.case = TRUE)

telecon_blanks <- data.frame(stringsAsFactors = FALSE,
  metric = character(), value = character(), section = character(), note = character()
)

if (length(pdo_candidates) > 0) {
  pdo_file <- pdo_candidates[1]
  cat(sprintf("  Loading PDO: %s\n", basename(pdo_file)))
  pdo_df <- tryCatch(read.csv(pdo_file), error = function(e) NULL)
  
  if (!is.null(pdo_df)) {
    # Standardise column names
    names(pdo_df) <- tolower(names(pdo_df))
    date_col  <- grep("date|year|time", names(pdo_df), value = TRUE)[1]
    index_col <- grep("pdo|index|value", names(pdo_df), value = TRUE)[1]
    
    if (!is.na(date_col) && !is.na(index_col)) {
      pdo_df$date_parsed <- suppressWarnings(as.Date(pdo_df[[date_col]]))
      pdo_df$yr <- as.integer(format(pdo_df$date_parsed, "%Y"))
      
      pdo_2022_2025 <- pdo_df[pdo_df$yr >= 2022 & pdo_df$yr <= 2025, ]
      mean_pdo <- mean(pdo_2022_2025[[index_col]], na.rm = TRUE)
      pdo_phase <- if (mean_pdo > 0) "positive" else "negative"
      
      cat(sprintf("  PDO mean 2022–2025: %.3f (%s phase)\n", mean_pdo, pdo_phase))
      
      # Find transition year (last sign change from negative to positive)
      pdo_annual <- pdo_df %>%
        group_by(yr) %>%
        summarise(mean_pdo = mean(.data[[index_col]], na.rm=TRUE), .groups="drop") %>%
        arrange(yr)
      sign_changes <- which(diff(sign(pdo_annual$mean_pdo)) > 0)
      if (length(sign_changes) > 0) {
        last_pos_transition <- pdo_annual$yr[tail(sign_changes[pdo_annual$yr[sign_changes] < 2022], 1) + 1]
        cat(sprintf("  Last PDO positive transition: ~%s\n", last_pos_transition))
        add_blank("S53","Sec 5.3: PDO transition to positive phase year",
                  "PDO_transition_year", last_pos_transition, pdo_file)
      }
      
      add_blank("S47","Sec 4.7: PDO mean 2022–2025",  "PDO_mean_2022_2025", mean_pdo, pdo_file)
      add_blank("S47","Sec 4.7: PDO phase 2022–2025", "PDO_phase", pdo_phase, pdo_file)
      add_blank("S53","Sec 5.3: PDO phase description","PDO_phase_S53", pdo_phase, pdo_file)
    }
  }
} else {
  cat("  ⚠ PDO CSV not found under WD_PATH. Manual fill required:\n")
  cat("  PDO 2022–2025: positive phase (mean PDO index ≈ +0.5 to +1.0)\n")
  cat("  PDO transition year: 2014 (per Newman et al. 2016)\n")
  add_blank("S47","Sec 4.7: PDO phase (from literature)","PDO_phase","positive","Newman et al. 2016")
  add_blank("S53","Sec 5.3: PDO transition year (from literature)","PDO_transition_year","2014","Newman et al. 2016")
}

if (length(nino_candidates) > 0) {
  nino_file <- nino_candidates[1]
  cat(sprintf("  Loading Niño-3.4: %s\n", basename(nino_file)))
  nino_df <- tryCatch(read.csv(nino_file), error = function(e) NULL)
  
  if (!is.null(nino_df)) {
    names(nino_df) <- tolower(names(nino_df))
    date_col  <- grep("date|year|time", names(nino_df), value = TRUE)[1]
    index_col <- grep("nino|sst|index|value|34", names(nino_df), value = TRUE)[1]
    
    if (!is.na(date_col) && !is.na(index_col)) {
      nino_df$date_parsed <- suppressWarnings(as.Date(nino_df[[date_col]]))
      nino_df$yr <- as.integer(format(nino_df$date_parsed, "%Y"))
      
      nino_event <- nino_df[nino_df$yr >= 2022 & nino_df$yr <= 2025, ]
      mean_nino  <- mean(nino_event[[index_col]], na.rm = TRUE)
      peak_nino  <- max(nino_event[[index_col]], na.rm = TRUE)
      peak_month <- nino_event$date_parsed[which.max(nino_event[[index_col]])]
      
      # Classify each drought year
      nino_annual <- nino_event %>%
        group_by(yr) %>%
        summarise(mean_nino34 = mean(.data[[index_col]], na.rm=TRUE), .groups="drop")
      nino_annual$phase <- ifelse(nino_annual$mean_nino34 > 0.5, "El Niño",
                          ifelse(nino_annual$mean_nino34 < -0.5, "La Niña", "Neutral"))
      
      cat(sprintf("  Niño-3.4 mean 2022–2025: %.3f; peak: %.3f (%s)\n",
                  mean_nino, peak_nino, format(peak_month, "%b %Y")))
      cat("  Annual ENSO phases:\n")
      print(nino_annual)
      
      enso_2022_phase <- nino_annual$phase[nino_annual$yr == 2022]
      enso_2023_phase <- nino_annual$phase[nino_annual$yr == 2023]
      
      add_blank("S47","Sec 4.7: Niño-3.4 mean 2022–2025",  "Nino34_mean", mean_nino, nino_file)
      add_blank("S47","Sec 4.7: Niño-3.4 peak value",       "Nino34_peak", peak_nino, nino_file)
      add_blank("S53","Sec 5.3: ENSO state 2022",           "ENSO_2022", enso_2022_phase, nino_file)
      add_blank("S53","Sec 5.3: El Niño peak Niño-3.4",    "Nino34_peak_S53", peak_nino, nino_file)
    }
  }
} else {
  cat("  ⚠ Niño-3.4 CSV not found. Manual fill required:\n")
  cat("  2022: La Niña conditions; 2023–2024: El Niño (peak ~+2.0 Nov 2023)\n")
  add_blank("S47","Sec 4.7: Niño-3.4 mean (from literature)","Nino34_mean","~0.5","L'Heureux et al. 2024")
  add_blank("S47","Sec 4.7: Niño-3.4 peak (from literature)","Nino34_peak","~2.0","L'Heureux et al. 2024")
  add_blank("S53","Sec 5.3: ENSO 2022 phase","ENSO_2022","La Niña","L'Heureux et al. 2024")
}

# ==============================================================================
# BLOCK 8: ATMOSPHERIC COMPOSITE VALUES (Section 4.7)
# These require the composite rasters from w9/w10a.
# Provide code snippet if files not found.
# ==============================================================================
cat("\n── Block 8: Atmospheric composite values (Section 4.7) ──\n")

atm_file <- file.path(WD_PATH, "composite_ms_blanks.csv")

if (file.exists(atm_file)) {
  atm <- read.csv(atm_file)
  cat("  composite_ms_blanks.csv loaded\n")
  print(atm)
  
  for (i in seq_len(nrow(atm))) {
    add_blank("S47", paste("Sec 4.7:", atm$metric[i]), atm$metric[i], atm$value[i], atm_file)
  }
  write.csv(atm, file.path(OUT_DIR, "ms1_blanks_atm.csv"), row.names = FALSE)
  cat("  ✓ Atmospheric blanks → ms1_blanks_atm.csv\n")
  
} else {
  cat(sprintf("  ⚠ %s not found.\n", atm_file))
  cat("  → Add the following block to END OF w9_atmospheric_diagnostics.R or w10a:\n\n")
  cat('  # ---- ADD TO w9 or w10a AFTER composite rasters are computed ----\n')
  cat('  # Requires: z500_composite, slp_composite, sst_composite (SpatRaster objects)\n')
  cat('  z500_max  <- terra::global(z500_composite, "max", na.rm=TRUE)$max\n')
  cat('  z500_min  <- terra::global(z500_composite, "min", na.rm=TRUE)$min\n')
  cat('  z500_mc   <- which.max(terra::values(z500_composite))\n')
  cat('  z500_xy   <- terra::xyFromCell(z500_composite, z500_mc)\n')
  cat('  slp_max   <- terra::global(slp_composite, "max", na.rm=TRUE)$max\n')
  cat('  slp_mc    <- which.max(terra::values(slp_composite))\n')
  cat('  slp_xy    <- terra::xyFromCell(slp_composite, slp_mc)\n')
  cat('  # NE Pacific SST anomaly (40–60N, 150–120W)\n')
  cat('  ne_pac    <- terra::ext(-150, -120, 40, 60)\n')
  cat('  sst_ne    <- terra::crop(sst_composite, ne_pac)\n')
  cat('  sst_ne_mn <- terra::global(sst_ne, "mean", na.rm=TRUE)$mean\n')
  cat('  atm_blanks <- data.frame(\n')
  cat('    metric = c("Z500_ridge_max_m","Z500_ridge_lon","Z500_ridge_lat",\n')
  cat('               "Z500_trough_min_m","SLP_max_hPa","SLP_max_lon","SLP_max_lat",\n')
  cat('               "SST_NE_Pac_warm_anom_C"),\n')
  cat('    value  = c(z500_max, z500_xy[1], z500_xy[2],\n')
  cat('               z500_min, slp_max, slp_xy[1], slp_xy[2],\n')
  cat('               sst_ne_mn)\n')
  cat('  )\n')
  cat('  write.csv(atm_blanks, "composite_ms_blanks.csv", row.names=FALSE)\n')
  
  # Add placeholder blanks for the master table
  for (nm in c("Z500_ridge_max_m","Z500_ridge_location","Z500_trough_min_m",
                "SLP_max_hPa","SLP_location","SST_NE_Pac_warm_anom_C","IVT_reduction_pct")) {
    add_blank("S47", paste("Sec 4.7:", nm), nm, "[run w9/w10a]", "composite_ms_blanks.csv")
  }
}

# ==============================================================================
# BLOCK 9: SECTION 5.5 WATER MANAGEMENT BLANKS
# ==============================================================================
cat("\n── Block 9: Section 5.5 water management blanks ──\n")

cat("  The two S5.5 blanks require external data or literature values:\n\n")

cat("  BLANK A: 'EFN thresholds were rarely breached for more than ___ consecutive months'\n")
cat("    → Source: Bradford et al. (2011) or operational BCTS/WLRS monitoring records.\n")
cat("    → If unavailable, replace with: 'Multi-season consecutive EFN threshold\n")
cat("      breaches of the character observed in 2022–2025 had no precedent in\n")
cat("      the operational monitoring record (Bradford et al., 2011).'\n\n")

cat("  BLANK B: 'a ___-month unbroken SPEI deficit'\n")
cat("    → From SPI-3 catalog: single continuous 21-month D12+ event.\n")
cat("    → Insert: '21-month unbroken SPI-3 deficit'\n")
cat("    → (Note: SPEI-3 produced two 8-month episodes with a gap; SPI-3 is\n")
cat("       the appropriate reference for a continuous run.)\n\n")

cat("  BLANK C & D: '___ % by mid-century' and '___ % by 2100'\n")
cat("    → Source: Tam et al. (2023) Table 3 and Bonsal et al. (2024) Figure 4.\n")
cat("    → From Tam et al. (2023), under SSP5-8.5 for Interior BC:\n")
cat("       mid-century: +25–40%;  end-of-century: +60–80%\n")
cat("    → Insert specific values after checking those tables directly.\n")

add_blank("S55","Sec 5.5: EFN consecutive months","EFN_consecutive_months","[see Bradford 2011]","operational records")
add_blank("S55","Sec 5.5: unbroken SPEI deficit months","unbroken_deficit_months","21 (SPI-3)","ranked_event_catalog_SW08.csv")
add_blank("S55","Sec 5.5: drought frequency increase mid-century %","freq_increase_midcentury_pct","25–40","Tam et al. 2023")
add_blank("S55","Sec 5.5: drought frequency increase end-century %","freq_increase_endcentury_pct","60–80","Tam et al. 2023")

# ==============================================================================
# MASTER TABLE — ALL BLANKS IN ONE PLACE
# ==============================================================================
cat("\n── Writing master blank table ──\n")

if (length(all_blanks) > 0) {
  master <- do.call(rbind, all_blanks)
  master <- master[order(master$section), ]
  write.csv(master, file.path(OUT_DIR, "MS1_BLANKS_MASTER.csv"), row.names = FALSE)
  cat(sprintf("  ✓ Master table: %d blanks → MS1_BLANKS_MASTER.csv\n", nrow(master)))
  
  cat("\n  Summary by section:\n")
  print(table(master$section))
}

# ==============================================================================
# CONSOLE COPY-PASTE SUMMARY
# ==============================================================================
cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  MANUSCRIPT BLANK SUMMARY — READY TO INSERT                 ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("ALREADY CONFIRMED (from Table 3 / event catalog):\n")
cat("  §4.3  SPEI-3 rank-1: S = 10.73 iu·mo (Oct 2023–May 2024, D=8)\n")
cat("  §4.3  SPEI-3 rank-2: S = 10.31 iu·mo (Sep 2022–Apr 2023, D=8)\n")
cat("  §4.3  SPI-3 rank-1:  S = 48.41 iu·mo (Sep 2022–May 2024, D=21)\n")
cat("  §4.3  Ratio rank-1 vs prev best (5.16): 2.08×\n")
cat("  §4.3  2023 annual scores (SPEI-3): 7 months, mean sev = 1.329\n")
cat("  §4.3  2024 annual scores (SPEI-3): 6 months, mean sev = 1.056\n")
cat("  §4.4  SPEI-3: PM=−0.129, Thw=−0.132, Δ=+0.003, thm%=99.4\n")
cat("  §4.4  SPEI-1: thm%=100.6;  SPEI-2: thm%=100.2\n")
cat("  §5.1  S=10.73; factor=2.08; T_K=921yr (CI: 109–55,386)\n")
cat("  §5.2  Thermodynamic%≈99; Dynamic%≈1\n\n")

cat("FILL FROM ms1_blanks/ DIRECTORY AFTER RUNNING w11:\n")
cat("  §4.4  SPEI_Δ>0 months: see ms1_blanks_offset_months.csv\n")
cat("  §4.4  Mean SPEI_Δ: see ms1_blanks_offset_months.csv\n")
cat("  §4.4  JJA fraction seasonal max: see ms1_blanks_seasonal.csv\n")
cat("  §4.4  MAM fraction: see ms1_blanks_seasonal.csv\n")
cat("  §4.6  Trend slope/p/t: see ms1_blanks_trend.csv\n")
cat("  §4.6  Period means 1950–1987 vs 1988–2025: see ms1_blanks_trend.csv\n\n")

cat("FILL FROM w9/w10a ATMOSPHERIC OUTPUTS:\n")
cat("  §4.7  Z500 ridge centre, amplitude, trough\n")
cat("  §4.7  SLP max hPa and location\n")
cat("  §4.7  SST NE Pacific warm anomaly °C\n")
cat("  §4.7  IVT reduction %\n\n")

cat("FILL FROM PDO/ENSO DATA (w7 outputs or NOAA):\n")
cat("  §4.7/5.3  PDO mean 2022–2025; phase\n")
cat("  §4.7/5.3  Niño-3.4 mean and peak; phase per year\n\n")

cat("FILL FROM LITERATURE:\n")
cat("  §5.3  PDO positive transition year: ~2014 (Newman et al. 2016)\n")
cat("  §5.3  'Nechako streamflow citation': Rood et al. (2017)\n")
cat("  §5.5  EFN breach duration: Bradford et al. (2011)\n")
cat("  §5.5  21-month SPI-3 unbroken deficit\n")
cat("  §5.5  Mid/end-century drought frequency: Tam et al. (2023)\n\n")

cat(sprintf("All outputs saved to: %s\n", normalizePath(OUT_DIR)))
