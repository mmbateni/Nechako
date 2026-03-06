####################################################################################
# w6_teleconnection_prep.R  ·  TELECONNECTION DATA PREPARATION
# Run ONCE before w7_moderated_mediation.R
####################################################################################
# Downloads and parses three NOAA teleconnection indices:
#   • ONI  (El Niño–Southern Oscillation, 3-month centred running mean of
#            ERSST.v5 Niño-3.4 anomalies)  — NOAA CPC
#   • PDO  (Pacific Decadal Oscillation, leading PC of monthly North Pacific
#            SST anomalies poleward of 20°N)  — NOAA PSL / JISAO
#   • PNA  (Pacific-North American pattern, normalised monthly index)  — NOAA CPC
#
# OUTPUTS (all to {WD_PATH}/teleconnections/):
#   oni_monthly.csv          — date, value
#   pdo_monthly.csv          — date, value
#   pna_monthly.csv          — date, value
#   teleconnections_merged.csv — wide format: date, ONI, PDO, PNA
#   teleconnections_QC_plots.pdf — diagnostic time series + cross-plots
#
# Run BEFORE: w7_moderated_mediation.R
# Depends on: DROUGHT_ANALYSIS_utils.R, utils_teleconnection_addon.R
####################################################################################

rm(list = ls())
gc()

source("DROUGHT_ANALYSIS_utils.R")
source("utils_teleconnection_addon.R")

utils_load_packages(c(
  "ggplot2", "patchwork", "dplyr", "lubridate",
  "openxlsx", "corrplot", "GGally"
))

setwd(WD_PATH)
dir.create(TELE_DIR, showWarnings = FALSE, recursive = TRUE)

cat("\n╔═══════════════════════════════════════════════╗\n")
cat("║  w6  TELECONNECTION DATA PREPARATION          ║\n")
cat("╚═══════════════════════════════════════════════╝\n\n")

# ── Analysis period (must match your drought index record) ───────────────────
START_YEAR <- 1950L
END_YEAR   <- 2024L

####################################################################################
# STEP 1: DOWNLOAD / LOAD EACH INDEX
####################################################################################
cat("── STEP 1: Downloading / loading indices ──\n\n")

# force_download = FALSE  uses cached files if they exist.
# Set to TRUE to refresh from NOAA.
oni_df <- load_teleconnection("oni", start_year = START_YEAR,
                               end_year = END_YEAR, force_download = FALSE)
pdo_df <- load_teleconnection("pdo", start_year = START_YEAR,
                               end_year = END_YEAR, force_download = FALSE)
pna_df <- load_teleconnection("pna", start_year = START_YEAR,
                               end_year = END_YEAR, force_download = FALSE)

# Rename value columns for clarity
colnames(oni_df) <- c("date", "ONI")
colnames(pdo_df) <- c("date", "PDO")
colnames(pna_df) <- c("date", "PNA")

####################################################################################
# STEP 2: MERGE INTO WIDE FORMAT
####################################################################################
cat("\n── STEP 2: Merging into wide format ──\n")

tc_wide <- Reduce(function(a, b) merge(a, b, by = "date", all = FALSE),
                  list(oni_df, pdo_df, pna_df))

tc_wide <- tc_wide[order(tc_wide$date), ]
tc_wide <- tc_wide[stats::complete.cases(tc_wide), ]

cat(sprintf("  ✓ Merged: %d months (%s – %s)\n",
            nrow(tc_wide), min(tc_wide$date), max(tc_wide$date)))
cat(sprintf("  Columns: %s\n", paste(names(tc_wide), collapse = ", ")))

# Save wide CSV — this is the primary input to w7
write.csv(tc_wide,
          file.path(TELE_DIR, "teleconnections_merged.csv"),
          row.names = FALSE)
cat("  ✓ Saved: teleconnections_merged.csv\n")

####################################################################################
# STEP 3: DESCRIPTIVE STATISTICS AND CORRELATION
####################################################################################
cat("\n── STEP 3: Descriptive statistics ──\n")

desc_stats <- data.frame(
  Index   = c("ONI", "PDO", "PNA"),
  N       = sapply(c("ONI","PDO","PNA"), function(v) sum(!is.na(tc_wide[[v]]))),
  Mean    = round(sapply(c("ONI","PDO","PNA"), function(v) mean(tc_wide[[v]], na.rm=TRUE)), 4),
  SD      = round(sapply(c("ONI","PDO","PNA"), function(v) sd(tc_wide[[v]], na.rm=TRUE)), 4),
  Min     = round(sapply(c("ONI","PDO","PNA"), function(v) min(tc_wide[[v]], na.rm=TRUE)), 4),
  Max     = round(sapply(c("ONI","PDO","PNA"), function(v) max(tc_wide[[v]], na.rm=TRUE)), 4),
  Skewness = round(sapply(c("ONI","PDO","PNA"), function(v) {
    x <- tc_wide[[v]][!is.na(tc_wide[[v]])]
    n <- length(x); m <- mean(x); s <- sd(x)
    (sum((x - m)^3) / n) / s^3
  }), 4),
  stringsAsFactors = FALSE
)

print(desc_stats)

# Pairwise Pearson and Spearman correlations
cat("\n  Pairwise correlations:\n")
cor_pearson  <- round(cor(tc_wide[, c("ONI","PDO","PNA")],
                          method = "pearson",  use = "complete.obs"), 3)
cor_spearman <- round(cor(tc_wide[, c("ONI","PDO","PNA")],
                          method = "spearman", use = "complete.obs"), 3)
cat("  Pearson:\n");  print(cor_pearson)
cat("  Spearman:\n"); print(cor_spearman)

# Save stats
stats_wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(stats_wb, "Descriptive")
openxlsx::writeData(stats_wb, "Descriptive", desc_stats)
openxlsx::addWorksheet(stats_wb, "Pearson_Cor")
openxlsx::writeData(stats_wb, "Pearson_Cor", as.data.frame(cor_pearson))
openxlsx::addWorksheet(stats_wb, "Spearman_Cor")
openxlsx::writeData(stats_wb, "Spearman_Cor", as.data.frame(cor_spearman))
openxlsx::saveWorkbook(stats_wb,
                       file.path(TELE_DIR, "teleconnections_stats.xlsx"),
                       overwrite = TRUE)
cat("  ✓ Saved: teleconnections_stats.xlsx\n")

####################################################################################
# STEP 4: QC DIAGNOSTIC PLOTS
####################################################################################
cat("\n── STEP 4: Generating QC plots ──\n")

# ── Long format for ggplot ────────────────────────────────────────────────────
tc_long <- tidyr::pivot_longer(tc_wide, cols = c("ONI","PDO","PNA"),
                                names_to = "Index", values_to = "Value")
# (tidyr fallback if not loaded)
if (!exists("tc_long") || !is.data.frame(tc_long)) {
  tc_long <- do.call(rbind, lapply(c("ONI","PDO","PNA"), function(idx) {
    data.frame(date  = tc_wide$date,
               Index = idx,
               Value = tc_wide[[idx]],
               stringsAsFactors = FALSE)
  }))
}

# Colour scheme consistent with index_colours in utils.R
tele_colours <- c(ONI = "#E74C3C", PDO = "#2980B9", PNA = "#27AE60")

# ── Fig 1: Full time series ────────────────────────────────────────────────
p_ts <- ggplot2::ggplot(tc_long, ggplot2::aes(x = date, y = Value, colour = Index)) +
  ggplot2::geom_line(linewidth = 0.6, alpha = 0.85) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  ggplot2::facet_wrap(~ Index, ncol = 1, scales = "free_y") +
  ggplot2::scale_colour_manual(values = tele_colours) +
  ggplot2::scale_x_date(date_breaks = "10 years", date_labels = "%Y") +
  shared_ts_theme(11) +
  ggplot2::labs(
    title    = "Teleconnection Indices — Nechako Moderated Mediation Analysis",
    subtitle = sprintf("Monthly values %d–%d | Source: NOAA CPC / PSL",
                       START_YEAR, END_YEAR),
    x = NULL, y = "Index Value"
  ) +
  ggplot2::theme(legend.position = "none",
                 strip.text      = ggplot2::element_text(face = "bold"))

# ── Fig 2: Scatter matrix with correlation ────────────────────────────────
p_scatter <- GGally::ggpairs(
  tc_wide[, c("ONI", "PDO", "PNA")],
  upper = list(continuous = GGally::wrap("cor", method = "spearman", size = 4.5)),
  lower = list(continuous = GGally::wrap("points", alpha = 0.25, size = 0.8)),
  diag  = list(continuous = GGally::wrap("densityDiag", fill = "steelblue", alpha = 0.4)),
  title = "Pairwise relationships among teleconnection predictors"
)

# ── Fig 3: ONI × PDO joint density — key moderation visualisation ─────────
p_oniXpdo <- ggplot2::ggplot(tc_wide, ggplot2::aes(x = ONI, y = PDO)) +
  ggplot2::geom_bin2d(bins = 35) +
  ggplot2::scale_fill_viridis_c(option = "plasma", name = "Count") +
  ggplot2::geom_hline(yintercept =  0.5, linetype = "dashed") +
  ggplot2::geom_hline(yintercept = -0.5, linetype = "dashed") +
  ggplot2::geom_vline(xintercept =  0.5, linetype = "dashed") +
  ggplot2::geom_vline(xintercept = -0.5, linetype = "dashed") +
  ggplot2::annotate("text", x = 1.8, y = 1.2,
                    label = "El Niño\nWarm PDO\n(max drought risk)",
                    size = 3.2, colour = "#C0392B", fontface = "bold") +
  ggplot2::annotate("text", x = -1.8, y = -1.2,
                    label = "La Niña\nCool PDO\n(wet conditions)",
                    size = 3.2, colour = "#1A5276", fontface = "bold") +
  shared_ts_theme(11) +
  ggplot2::labs(
    title    = "Joint distribution of ONI and PDO",
    subtitle = "Quadrant labels indicate compound drought / wet forcing",
    x = "ONI (El Niño–La Niña axis)", y = "PDO Index"
  )

# ── Fig 4: Monthly climatology of each index ─────────────────────────────
tc_clim <- tc_long %>%
  dplyr::mutate(Month = as.integer(format(date, "%m"))) %>%
  dplyr::group_by(Index, Month) %>%
  dplyr::summarise(Mean = mean(Value, na.rm = TRUE),
                   SD   = sd(Value,   na.rm = TRUE),
                   .groups = "drop")

p_clim <- ggplot2::ggplot(tc_clim,
                           ggplot2::aes(x = Month, y = Mean,
                                        colour = Index, fill = Index)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = Mean - SD, ymax = Mean + SD),
                       alpha = 0.15, colour = NA) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  ggplot2::scale_colour_manual(values = tele_colours) +
  ggplot2::scale_fill_manual(values = tele_colours) +
  ggplot2::scale_x_continuous(breaks = 1:12,
                               labels = MONTH_NAMES) +
  shared_ts_theme(11) +
  ggplot2::labs(
    title    = "Seasonal climatology of teleconnection indices",
    subtitle = "Mean ± 1 SD by calendar month",
    x = "Month", y = "Mean Index Value"
  )

# ── Save PDF ──────────────────────────────────────────────────────────────
pdf_file <- file.path(TELE_DIR, "teleconnections_QC_plots.pdf")
if (safe_pdf(pdf_file, width = 13, height = 9)) {
  print(p_ts)
  print(p_scatter)
  print(p_oniXpdo)
  print(p_clim)
  grDevices::dev.off()
  cat(sprintf("  ✓ Saved: teleconnections_QC_plots.pdf (%d pages)\n", 4))
}

####################################################################################
# STEP 5: PDO PHASE CLASSIFICATION TABLE
####################################################################################
cat("\n── STEP 5: PDO phase classification ──\n")

tc_wide$pdo_phase <- cut(tc_wide$PDO,
                          breaks = c(-Inf, -0.5, 0.5, Inf),
                          labels = c("Cool", "Neutral", "Warm"),
                          right  = FALSE)

phase_summary <- tc_wide %>%
  dplyr::group_by(pdo_phase) %>%
  dplyr::summarise(
    N_months    = dplyr::n(),
    Pct_total   = round(100 * dplyr::n() / nrow(tc_wide), 1),
    Mean_ONI    = round(mean(ONI, na.rm = TRUE), 3),
    Mean_PNA    = round(mean(PNA, na.rm = TRUE), 3),
    Mean_PDO    = round(mean(PDO, na.rm = TRUE), 3),
    .groups = "drop"
  )

print(phase_summary)
write.csv(phase_summary,
          file.path(TELE_DIR, "pdo_phase_summary.csv"),
          row.names = FALSE)
cat("  ✓ Saved: pdo_phase_summary.csv\n")

####################################################################################
# STEP 6: SAVE FINAL MERGED WIDE CSV (rewrite with phase column)
####################################################################################
write.csv(tc_wide,
          file.path(TELE_DIR, "teleconnections_merged.csv"),
          row.names = FALSE)

cat("\n╔═══════════════════════════════════════════════╗\n")
cat("║  w6 COMPLETE                                  ║\n")
cat("╚═══════════════════════════════════════════════╝\n")
cat("  Next: run w7_moderated_mediation.R\n\n")
