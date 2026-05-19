##############################################################################
# FigS_Seasonality_P_PET.R
# ─────────────────────────────────────────────────────────────────────────────
# PURPOSE  Two-panel bar chart showing the seasonal cycle of the Nechako Basin:
#   Panel (a): mean monthly precipitation (P) vs Penman-Monteith PET (mm month⁻¹)
#              — exactly 12 bars per variable, averaged over all 1950–2025 years.
#   Panel (b): climatic water balance P − PET (mm month⁻¹).
#
# DATA SOURCES
#   Precipitation : temporal_drought/basin_averaged_timeseries/spi_01_basin_average.csv
#                   (written by w2_basin_timeseries.R; values in mm month⁻¹)
#   PM-PET        : monthly_data_direct/ERA5Land_Nechako_PET_monthly_summary.csv
#                   (written by 2preq_PET_ERALand.R; mean_pet column in mm day⁻¹;
#                    converted here to mm month⁻¹ using average days per month)
#
# PREREQ  2preq_PET_ERALand.R and w2_basin_timeseries.R (Part 2) must have run.
# OUTPUT  FigS_Seasonality_P_PET.pdf  +  .png  in BASIN_PLOT_DIR
##############################################################################
source("DROUGHT_ANALYSIS_utils.R")
utils_load_packages(c("ggplot2", "dplyr", "patchwork", "data.table"))
if (!dir.exists(WD_PATH)) stop("Working directory not found: ", WD_PATH)
setwd(WD_PATH)

BASIN_TS_DIR   <- file.path(WD_PATH, "temporal_drought", "basin_averaged_timeseries")
BASIN_PLOT_DIR <- file.path(WD_PATH, "temporal_drought", "basin_averaged_plots")
PET_CSV        <- file.path(WD_PATH, "monthly_data_direct",
                            "ERA5Land_Nechako_PET_monthly_summary.csv")
PRECIP_CSV     <- file.path(BASIN_TS_DIR, "spi_01_basin_average.csv")
dir.create(BASIN_PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

MON_LABELS   <- month.abb                                # Jan … Dec
# Mean days per calendar month (non-leap average for mm day⁻¹ → mm month⁻¹)
DAYS_PER_MON <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# ─── 1. PRECIPITATION: 12 monthly means (mm month⁻¹) ─────────────────────────
if (!file.exists(PRECIP_CSV))
  stop("Precip CSV not found: ", PRECIP_CSV,
       "\n  Run w2_basin_timeseries.R (Part 2) first.")

prec_raw       <- as.data.frame(data.table::fread(PRECIP_CSV))
prec_raw$date  <- as.Date(prec_raw$date)
prec_raw$month <- as.integer(format(prec_raw$date, "%m"))

cm_prec <- prec_raw |>
  dplyr::group_by(month) |>
  dplyr::summarise(mean_val = mean(value, na.rm = TRUE),
                   sd_val   = sd(value,   na.rm = TRUE),
                   n_yrs    = dplyr::n(),
                   .groups  = "drop") |>
  dplyr::mutate(month_lab = factor(MON_LABELS[month], levels = MON_LABELS),
                variable  = "Precipitation (P)")

cat(sprintf("  Precip: %d yrs, mean annual = %.0f mm yr\u207b\u00b9\n",
            cm_prec$n_yrs[1], sum(cm_prec$mean_val)))

# ─── 2. PM-PET: 12 monthly means (mm month⁻¹) ────────────────────────────────
# ERA5Land_Nechako_PET_monthly_summary.csv
#   columns: date (YYYY-MM), year, month, mean_pet [mm day⁻¹], ...
if (!file.exists(PET_CSV))
  stop("PET CSV not found: ", PET_CSV,
       "\n  Run 2preq_PET_ERALand.R first.")

pet_raw       <- as.data.frame(data.table::fread(PET_CSV))
# date can be "YYYY-MM" or "YYYY-MM-DD"
if (!inherits(pet_raw$date, "Date"))
  pet_raw$date <- as.Date(paste0(substr(as.character(pet_raw$date), 1, 7), "-01"))
pet_raw$month  <- as.integer(format(pet_raw$date, "%m"))
pet_raw$pet_mm <- pet_raw$mean_pet * DAYS_PER_MON[pet_raw$month]  # mm month⁻¹

cm_pet <- pet_raw |>
  dplyr::group_by(month) |>
  dplyr::summarise(mean_val = mean(pet_mm, na.rm = TRUE),
                   sd_val   = sd(pet_mm,   na.rm = TRUE),
                   n_yrs    = dplyr::n(),
                   .groups  = "drop") |>
  dplyr::mutate(month_lab = factor(MON_LABELS[month], levels = MON_LABELS),
                variable  = "PM-PET")

cat(sprintf("  PM-PET: %d yrs, mean annual = %.0f mm yr\u207b\u00b9\n",
            cm_pet$n_yrs[1], sum(cm_pet$mean_val)))
cat(sprintf("  Aridity index (PET/P) = %.2f\n",
            sum(cm_pet$mean_val) / sum(cm_prec$mean_val)))

# ─── 3. Panel (a): grouped bar chart P vs PM-PET ─────────────────────────────
COL_P   <- "#4393c3"
COL_PET <- "#d7301f"

cm_all          <- rbind(cm_prec, cm_pet)
cm_all$variable <- factor(cm_all$variable,
                          levels = c("Precipitation (P)", "PM-PET"))

nyrs_label <- cm_prec$n_yrs[1]

p1 <- ggplot2::ggplot(
  cm_all,
  ggplot2::aes(x = month_lab, y = mean_val,
               fill = variable, group = variable)) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(0.72),
    width = 0.65, colour = "white", linewidth = 0.15) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
    position = ggplot2::position_dodge(0.72),
    width = 0.22, colour = "grey30", linewidth = 0.40) +
  ggplot2::scale_fill_manual(
    values = c("Precipitation (P)" = COL_P, "PM-PET" = COL_PET),
    name   = NULL) +
  ggplot2::scale_y_continuous(
    name   = expression(mm~month^{-1}),
    expand = ggplot2::expansion(mult = c(0, 0.09))) +
  ggplot2::labs(
    title    = "(a)  Mean monthly precipitation and Penman\u2013Monteith PET",
    subtitle = paste0(
      nyrs_label, "-year mean (1950\u20132025)  |  Error bars = \u00b11 SD  |  ",
      "Annual totals:  P = ", round(sum(cm_prec$mean_val)), " mm;  ",
      "PM-PET = ", round(sum(cm_pet$mean_val)), " mm"),
    x = NULL) +
  ggplot2::theme_classic(base_size = 9.5) +
  ggplot2::theme(
    legend.position    = c(0.86, 0.88),
    legend.background  = ggplot2::element_rect(fill = "white", colour = NA),
    legend.text        = ggplot2::element_text(size = 8.5),
    legend.key.size    = ggplot2::unit(0.40, "cm"),
    plot.title         = ggplot2::element_text(size = 9.5, face = "bold"),
    plot.subtitle      = ggplot2::element_text(size = 7.5, colour = "grey40"),
    axis.text          = ggplot2::element_text(size = 8.5),
    axis.title.y       = ggplot2::element_text(size = 8.5),
    panel.grid.major.y = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
    plot.margin        = ggplot2::margin(4, 6, 2, 4))

# ─── 4. Panel (b): water balance P − PET ──────────────────────────────────────
cm_wb <- dplyr::inner_join(
  cm_prec[, c("month","month_lab","mean_val")],
  cm_pet[,  c("month","mean_val")],
  by = "month", suffix = c("_p","_pet")) |>
  dplyr::mutate(wb       = mean_val_p - mean_val_pet,
                fill_col = ifelse(wb >= 0, COL_P, COL_PET))

p2 <- ggplot2::ggplot(
  cm_wb,
  ggplot2::aes(x = month_lab, y = wb, fill = fill_col)) +
  ggplot2::geom_col(colour = "white", linewidth = 0.15, width = 0.65) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.55) +
  ggplot2::scale_fill_identity() +
  ggplot2::scale_y_continuous(
    name   = expression(P - PET~~(mm~month^{-1})),
    expand = ggplot2::expansion(mult = c(0.12, 0.10))) +
  # Surplus / deficit labels placed at midpoints of each group of bars
  ggplot2::annotate("text",
                    x     = mean(which(cm_wb$wb > 0)),
                    y     = max(cm_wb$wb[cm_wb$wb > 0], na.rm = TRUE) * 0.82,
                    label = "P surplus", colour = COL_P,
                    size = 3.0, fontface = "italic") +
  ggplot2::annotate("text",
                    x     = mean(which(cm_wb$wb < 0)),
                    y     = min(cm_wb$wb[cm_wb$wb < 0], na.rm = TRUE) * 0.78,
                    label = "PET deficit", colour = COL_PET,
                    size = 3.0, fontface = "italic") +
  ggplot2::labs(
    title    = "(b)  Climatic water balance  P \u2212 PM-PET",
    subtitle = "Blue = monthly surplus (P > PET);  red = monthly deficit (PET > P)",
    x        = "Calendar month") +
  ggplot2::theme_classic(base_size = 9.5) +
  ggplot2::theme(
    plot.title         = ggplot2::element_text(size = 9.5, face = "bold"),
    plot.subtitle      = ggplot2::element_text(size = 7.5, colour = "grey40"),
    axis.text          = ggplot2::element_text(size = 8.5),
    axis.title.y       = ggplot2::element_text(size = 8.5),
    panel.grid.major.y = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
    plot.margin        = ggplot2::margin(2, 6, 4, 4))

# ─── 5. Assemble and save ─────────────────────────────────────────────────────
fig_seas <- p1 / p2 +
  patchwork::plot_annotation(
    title    = paste0(
      "Nechako River Basin \u2014 Seasonal cycle of precipitation ",
      "and Penman\u2013Monteith PET"),
    subtitle = paste0(
      "Each bar = mean of ", nyrs_label, " years (1950\u20132025).  ",
      "Annual totals: P = ", round(sum(cm_prec$mean_val)), " mm;  ",
      "PM-PET = ", round(sum(cm_pet$mean_val)), " mm;  ",
      "aridity index = ", round(sum(cm_pet$mean_val) / sum(cm_prec$mean_val), 2), "."),
    theme = ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 11, face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 8.5, colour = "grey35",
                                            hjust = 0,
                                            margin = ggplot2::margin(b = 4))))

cat("\n── Saving FigS_Seasonality_P_PET...\n")
for (ext in c("pdf", "png")) {
  out_f <- file.path(BASIN_PLOT_DIR, paste0("FigS_Seasonality_P_PET.", ext))
  tryCatch(
    ggplot2::ggsave(out_f, fig_seas,
                    width  = 7.2,
                    height = 7.0,
                    units  = "in",
                    dpi    = if (ext == "png") 300 else "print",
                    device = ext),
    error = function(e) cat(sprintf("  \u26a0 %s: %s\n", ext, e$message)))
  cat(sprintf("  \u2713 %s\n", basename(out_f)))
}
cat("\n\u2713 FigS_Seasonality_P_PET.R complete.\n")
cat(sprintf("  Output: %s\n", normalizePath(BASIN_PLOT_DIR)))