# ####################################################################################
# 4pr_pet_trends_visualization.r  ·  Trend Maps + Temporal Change + Dual-PET Plots
# ─────────────────────────────────────────────────────────────────────────────────
# Reads results produced by 4pr_pet_trends.r (4-variable version) and generates:
#   BASIN-WIDE MAPS  (Pr, PET_PM, PET_Thw, Temperature)
# • Publication maps:  tau, Sen's slope, direction  (VC & TFPW) — all 4 vars
# • Comparison map:    Pr / PET_PM / PET_Thw / Tair side-by-side
# • Regime shift maps: decade of first changepoint — all 4 vars
# • Spectral analysis maps: dominant period — all 4 vars
# • Spatial variability maps: CV of Sen's slope — all 4 vars
# • Method comparison:  VC vs TFPW scatter + agreement bar charts
# [NEW] CALENDAR-MONTHLY TEMPORAL CHANGE PLOTS  (Functions 13a–13d)
# For each variable, a 12-panel figure where each panel shows the inter-annual
# time series of all values for that calendar month (e.g. all Januaries 1950–
#                                                    2025), with OLS trend line + significance annotation.
# 13a  Precipitation — 12-panel calendar-month temporal series
# 13b  PET Penman-Monteith — same layout
# 13c  PET Thornthwaite — same layout
# 13d  Temperature — same layout
# 13e  Four-variable comparison:  Sen's slope profile (all months × all vars)
# [NEW] PET_PM vs PET_Thw COMPARISON  (Function 14)
# 14a  Time series overlay (basin averages, monthly + 12-mo rolling mean)
# 14b  Per-month OLS slope comparison bar chart (PM vs Thw side by side)
# 14c  Basin trend maps side-by-side for every calendar month (Sen's slope)
# 14d  Bias (PM − Thw) barplot by calendar month
# 14e  Per-month scatter: PM value vs Thw value coloured by year
# SPECIFIC-POINT PLOTS
# • 4-panel time series per point: Pr / PET_PM / PET_Thw / Temperature
# • Per-month trend bars: all 4 variables
# • Trend statistics summary CSV + console print
# BASIN TIME SERIES  (Function 8 — extended for PET_Thw)
# • Pr + PET_PM + PET_Thw (triple-panel) + Temperature (fourth panel)
# TEMPERATURE DEDICATED ANALYSIS  (Function 12, unchanged)
####################################################################################
####################################################################################
library(terra)
library(data.table)
library(sf)
library(zoo)
library(ggplot2)
library(patchwork)
library(scales)
setwd("D:/Nechako_Drought/Nechako/")
out_dir <- "trend_analysis_pr_pet"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)
####################################################################################
# ── LOAD DATA
####################################################################################
cat("📦 Loading results, basin boundary, and metadata...\n")
all_results    <- readRDS(file.path(out_dir, "all_results.rds"))
basin_boundary <- readRDS(file.path(out_dir, "basin_boundary.rds"))
metadata       <- readRDS(file.path(out_dir, "analysis_metadata.rds"))
names(all_results)                <- trimws(names(all_results))
names(metadata$basin_avg_monthly) <- trimws(names(metadata$basin_avg_monthly))
names(metadata$basin_avg_annual)  <- trimws(names(metadata$basin_avg_annual))
basin_avg_monthly <- metadata$basin_avg_monthly
basin_avg_annual  <- metadata$basin_avg_annual
SPECIFIC_PTS      <- metadata$specific_pts
pt_trend_file <- file.path(out_dir, "point_trend_stats.csv")
pt_ts_file    <- file.path(out_dir, "point_monthly_timeseries.csv")
if (!file.exists(pt_trend_file))
  warning("point_trend_stats.csv not found — skipping specific-point plots")
if (!file.exists(pt_ts_file))
  warning("point_monthly_timeseries.csv not found — skipping specific-point TS plots")
point_trend_stats <- if (file.exists(pt_trend_file)) fread(pt_trend_file) else NULL
point_monthly_ts  <- if (file.exists(pt_ts_file))    fread(pt_ts_file)    else NULL
clipped_template_file <- file.path(out_dir, "clipped_template.rds")
if (!file.exists(clipped_template_file))
  stop("clipped_template.rds not found — re-run 4pr_pet_trends.r")
clipped_template <- readRDS(clipped_template_file)
template_bc      <- clipped_template
cat(sprintf("✓ Template: %d x %d, res=%.0f m\n",
            nrow(template_bc), ncol(template_bc), res(template_bc)[1]))
####################################################################################
# ── SANITY CHECKS
####################################################################################
.vars_present <- unique(trimws(all_results$variable))
cat("\n📊 Variables in all_results:", paste(.vars_present, collapse=", "), "\n")
has_pet_thw <- "PET_Thw" %in% .vars_present
if (!has_pet_thw)
  cat("\n⚠️  PET_Thw NOT in all_results — re-run 4pr_pet_trends.r with Thw enabled.\n",
      "   Dual-PET comparison plots (Function 14) will be skipped.\n\n")
n_coords_ok <- sum(!all_results$is_basin_average & !is.na(all_results$x) & !is.na(all_results$y))
cat(sprintf("📍 Pixel rows with valid x/y: %d\n", n_coords_ok))
if (n_coords_ok == 0) stop("All pixel coordinates are NA — re-run the processing script.")
required_cols   <- c("x", "y", "variable", "period", "tau_vc", "p_value_vc", "sl_vc",
                     "tau_tfpw", "p_value_tfpw", "filtered_vc", "is_basin_average")
missing_cols    <- required_cols[!required_cols %in% names(all_results)]
if (length(missing_cols))
  stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse=", ")))
####################################################################################
# ── HELPER FUNCTIONS
####################################################################################
create_raster_from_table <- function(results_dt, template, value_col) {
  value_col <- trimws(value_col)
  if (!value_col %in% names(results_dt))
    stop(sprintf("Column '%s' not found.", value_col))
  results_dt <- results_dt[!is.na(x) & !is.na(y)]
  if (!nrow(results_dt)) { r <- template; values(r) <- NA; return(r) }
  pts <- vect(as.data.frame(results_dt), geom=c("x","y"), crs=crs(template))
  r   <- rasterize(pts, template, field=value_col, touches=TRUE)
  mask(r, vect(basin_boundary))
}
smooth_raster <- function(r) {
  if (all(is.na(values(r)))) return(r)
  focal(r, w=matrix(1,3,3), fun=mean, na.rm=TRUE)
}
plot_raster_panel <- function(r, main, zlim=NULL, col, breaks=NULL,
                              legend=TRUE, categorical=FALSE, legend_title="") {
  valid_pct <- 0
  if (!is.null(r)) {
    vals <- values(r)
    valid_pct <- sum(!is.na(vals)) / length(vals) * 100
  }
  if (!categorical) r <- smooth_raster(r)
  if (is.null(r) || all(is.na(values(r)))) {
    plot.new()
    text(0.5, 0.5, "No valid trend data\n(Pixels filtered by quality criteria)",
         cex=1.1, col="gray40", font=3)
    return(invisible(NULL))
  }
  plg_args <- list(cex=0.9)
  if (legend_title != "") plg_args$title <- legend_title
  if (!is.null(breaks))
    plot(r, main=main, cex.main=1.0, col=col, breaks=breaks,
         axes=FALSE, box=FALSE, legend=legend, plg=plg_args, colNA=NA)
  else
    plot(r, main=main, cex.main=1.0, col=col, zlim=zlim,
         axes=FALSE, box=FALSE, legend=legend, plg=plg_args, colNA=NA)
  plot(st_geometry(basin_boundary), add=TRUE, col=NA, border="black", lwd=2.5)
  if (valid_pct < 10) {
    mtext("Note: Most pixels filtered due to data quality criteria\n(see methodology: ties, min-values, autocorrelation)",
          side=1, line=0.5, cex=0.7, col="gray40")
  }
}
.pdf_safe <- function(pdf_path, width, height, family='Helvetica', expr_fn) {
  opened <- FALSE
  tryCatch({
    pdf(pdf_path, width=width, height=height, family=family)
    opened <- TRUE
    expr_fn()
    dev.off()
    opened <- FALSE
    invisible(pdf_path)
  }, error = function(e) {
    if (opened && dev.cur() > 1) dev.off()
    cat(sprintf('  ⚠ PDF failed (%s): %s\n', basename(pdf_path), conditionMessage(e)))
    invisible(NULL)
  })
}
prep_annual   <- function(var_name) {
  dt   <- all_results[variable==var_name & period=="annual" & !filtered_vc & !is_basin_average]
  if (!nrow(dt)) return(NULL)
  dt[, combined_vc    := fifelse(p_value_vc < 0.05, abs(tau_vc), 0)]
  dt[, combined_tfpw  := fifelse(p_value_tfpw < 0.05, abs(tau_tfpw), 0)]
  dt[, direction_vc   := fifelse(p_value_vc < 0.05 & tau_vc < 0, -1L,
                                 fifelse(p_value_vc < 0.05 & tau_vc > 0, 1L, 0L))]
  dt[, direction_tfpw := fifelse(p_value_tfpw < 0.05 & tau_tfpw < 0, -1L,
                                 fifelse(p_value_tfpw < 0.05 & tau_tfpw > 0, 1L, 0L))]
  dt
}
cat("⚙️  Preparing trend metrics...\n")
precip_annual  <- prep_annual("Precipitation")
pet_annual     <- prep_annual("PET")
pet_thw_annual <- if (has_pet_thw) prep_annual("PET_Thw") else NULL
tair_annual    <- prep_annual("Temperature")
make_rasters   <- function(dt, prefix) {
  if (is.null(dt)) return(NULL)
  list(
    comb_vc   = create_raster_from_table(dt, template_bc, "combined_vc"),
    comb_tfpw = create_raster_from_table(dt, template_bc, "combined_tfpw"),
    dir_vc    = create_raster_from_table(dt, template_bc, "direction_vc"),
    dir_tfpw  = create_raster_from_table(dt, template_bc, "direction_tfpw"),
    mag_vc    = create_raster_from_table(dt, template_bc, "sl_vc"),
    mag_tfpw  = create_raster_from_table(dt, template_bc, "sl_tfpw"),
    tau_vc    = create_raster_from_table(dt, template_bc, "tau_vc"),
    tau_tfpw  = create_raster_from_table(dt, template_bc, "tau_tfpw")
  )
}
cat("🖼️  Generating rasters (Pr, PET_PM, PET_Thw, Temperature)...\n")
r_precip   <- make_rasters(precip_annual, "precip")
r_pet      <- make_rasters(pet_annual, "pet")
r_pet_thw  <- if (has_pet_thw) make_rasters(pet_thw_annual, "pet_thw") else NULL
r_tair     <- make_rasters(tair_annual, "tair")
####################################################################################
# FUNCTION 1: Publication map per variable
####################################################################################
create_publication_figure <- function(var_name, method="vc") {
  rlist <- switch(var_name,
                  PET=r_pet, PET_Thw=r_pet_thw,
                  Precipitation=r_precip, Temperature=r_tair, NULL)
  if (is.null(rlist)) { cat("⚠️ No rasters for", var_name, "— skipping\n"); return(NULL) }
  pdf_path   <- file.path(out_dir,
                          sprintf("%s_%s_publication_map.pdf",
                                  gsub("_", "",tolower(var_name)), method))
  r_comb   <- if (method=="vc") rlist$comb_vc else rlist$comb_tfpw
  r_dir    <- if (method=="vc") rlist$dir_vc else rlist$dir_tfpw
  r_mag    <- if (method=="vc") rlist$mag_vc else rlist$mag_tfpw
  r_tau    <- if (method=="vc") rlist$tau_vc else rlist$tau_tfpw
  result  <- .pdf_safe(pdf_path, width=10, height=6.67, expr_fn=function() {
    par(mfrow=c(1,3), mar=c(2,1.5,2.5,1), oma=c(1,1,3,1))
    if (var_name == "Temperature") {
      plot_raster_panel(r_tau, "Kendall's Tau", c(-0.3, 0.3),
                        hcl.colors(101, "RdBu",rev=TRUE), legend_title="τ")
      plot_raster_panel(r_mag, "Sen's Slope (°C/yr)", range(values(r_mag), na.rm=TRUE),
                        hcl.colors(101, "RdBu",rev=TRUE), legend_title="°C/yr")
      plot_raster_panel(r_dir, "Direction (p < 0.05)", NULL,
                        c("#4575b4", "#f0f0f0", "#d73027"),
                        breaks=c(-1.5,-0.5,0.5,1.5), legend=FALSE, categorical=TRUE)
      legend("bottomright", legend=c("Cooling", "Non-sig.", "Warming"),
             fill=c("#4575b4", "#f0f0f0", "#d73027"), bty="n", cex=1.1)
    } else if (grepl("^PET", var_name)) {
      plot_raster_panel(r_comb, "Combined Trend Strength", c(0,0.25),
                        hcl.colors(101, "YlOrRd"), legend_title="|τ|")
      mtext("Shows |tau| only for significant trends (p < 0.05); others zero",
            side=3, line=0.1, cex=0.7, col="gray50")
      plot_raster_panel(r_mag, "Sen's Slope (mm/yr)", c(0,1.5),
                        hcl.colors(101, "YlOrRd"), legend_title="mm/yr")
      mtext("Magnitude shown for all valid pixels (significance in other panels)",
            side=1, line=0.1, cex=0.7, col="gray40")
      plot_raster_panel(r_dir, "Direction (p < 0.05)", NULL,
                        c("#4575b4", "#f0f0f0", "#d73027"),
                        breaks=c(-1.5,-0.5,0.5,1.5), legend=FALSE, categorical=TRUE)
      legend("bottomright", legend=c("Decreasing", "Non-sig.", "Increasing"),
             fill=c("#4575b4", "#f0f0f0", "#d73027"), bty="n", cex=1.1)
    } else {
      plot_raster_panel(r_comb, "Combined Trend Strength", c(-0.15,0.15),
                        hcl.colors(101, "RdBu",rev=TRUE), legend_title="τ (sig.)")
      mtext("Shows |tau| only for significant trends (p < 0.05); others zero",
            side=3, line=0.1, cex=0.7, col="gray50")
      plot_raster_panel(r_mag, "Sen's Slope (mm/yr)", c(-2,2),
                        hcl.colors(101, "RdBu",rev=TRUE), legend_title="mm/yr")
      plot_raster_panel(r_tau, "Kendall's Tau", c(-0.3,0.3),
                        hcl.colors(101, "RdBu",rev=TRUE), legend_title="τ")
    }
    mtext(sprintf("%s Annual Trends (1950-2025) | %s Method",
                  var_name, ifelse(method=="vc", "Variance-Corrected MK", "TFPW-Corrected MK")),
          outer=TRUE, cex=1.1, font=2, line=0.5)
  })
  if (!is.null(result)) cat(sprintf("✅ %s %s map: %s\n", var_name, method, basename(pdf_path)))
  invisible(result)
}
####################################################################################
# FUNCTION 2: 4-variable comparison map
####################################################################################
create_comparison_figure <- function() {
  vars_avail <- Filter(function(v) {
    r_name <- paste0("r_", v)
    if (!exists(r_name)) return(FALSE)
    r_obj <- get(r_name)
    !is.null(r_obj) && !is.null(r_obj$tau_vc) &&
      sum(!is.na(values(r_obj$tau_vc)), na.rm=TRUE) > 0
  }, c(precip="precip", pet="pet", pet_thw="pet_thw", tair="tair"))
  r_list <- lapply(vars_avail, function(v) get(paste0("r_",v))$tau_vc)
  var_labels <- c(precip="Precipitation τ", pet="PET_PM τ",
                  pet_thw="PET_Thw τ", tair="Temperature τ")
  missing_vars <- setdiff(c("precip", "pet", "pet_thw", "tair"), names(vars_avail))
  n <- length(r_list)
  if (n == 0) { cat("⚠️ No variables available for comparison map.\n"); return(NULL) }
  pdf_name <- if (n < 4)
    file.path(out_dir, sprintf("%s_variable_trend_comparison.pdf", n))
  else
    file.path(out_dir, "four_variable_trend_comparison.pdf")
  result <- .pdf_safe(pdf_name, width=5*n, height=6, expr_fn=function() {
    par(mfrow=c(1,n), mar=c(2,1.5,2.5,1), oma=c(1,1,3,1))
    for (v in names(r_list))
      plot_raster_panel(r_list[[v]], var_labels[v], c(-0.3,0.3),
                        hcl.colors(101,"RdBu",rev=TRUE), legend_title="τ")
    main_title <- "Nechako Basin Annual Trends (1950-2025) | VC-MK"
    if (n < 4) {
      main_title <- paste0(main_title, "\n",
                           paste("Missing:", paste(toupper(missing_vars), collapse=", "), "(filtered/no data)"))
    }
    mtext(main_title, outer=TRUE, cex=1.2, font=2, line=0.5)
  })
  if (!is.null(result)) cat(sprintf("✅ %s-variable comparison map: %s\n", n, basename(pdf_name)))
  invisible(result)
}
####################################################################################
# FUNCTION 3: Regime shift maps (all 4 variables)
####################################################################################
create_regime_shift_maps <- function() {
  vars_here <- intersect(c("Precipitation", "PET", "PET_Thw", "Temperature"),
                         .vars_present)
  pdf_path <- file.path(out_dir, "regime_shifts_changepoints.pdf")
  result <- .pdf_safe(pdf_path, width=5*length(vars_here), height=5, expr_fn=function() {
    par(mfrow=c(1,length(vars_here)), mar=c(2,1.5,2.5,1), oma=c(1,1,3,1))
    for (var in vars_here) {
      var_data <- all_results[variable==var & period=="annual" &
                                !is_basin_average & changepoint_detected==TRUE]
      if (!nrow(var_data)) {
        plot.new(); text(0.5,0.5,paste("No regime shifts\n",var),cex=1.2); next
      }
      var_data[, shift_decade_num := as.numeric(cut(
        first_changepoint_year, breaks=seq(1950,2030,by=10), include.lowest=TRUE))]
      r_decade <- create_raster_from_table(var_data, template_bc, "shift_decade_num")
      plot_raster_panel(r_decade, paste0(var, ": Regime Shift Decade"),
                        c(1,8), hcl.colors(8, "Spectral",rev=TRUE))
    }
    mtext("Regime Shift Detection (PELT) — 1950-2025",
          outer=TRUE, cex=1.1, font=2, line=0.5)
  })
  if (!is.null(result)) cat(sprintf("✅ Regime shift map: %s\n", basename(pdf_path)))
  invisible(result)
}
####################################################################################
# FUNCTION 4: Spectral maps (all 4 variables)
####################################################################################
create_spectral_analysis_maps <- function() {
  vars_here <- intersect(c("Precipitation", "PET", "PET_Thw", "Temperature"),
                         .vars_present)
  pdf_path <- file.path(out_dir, "spectral_analysis.pdf")
  result <- .pdf_safe(pdf_path, width=5*length(vars_here), height=5, expr_fn=function() {
    par(mfrow=c(1,length(vars_here)), mar=c(2,1.5,2.5,1), oma=c(1,1,3,1))
    for (var in vars_here) {
      var_data <- all_results[variable==var & period=="annual" &
                                !is_basin_average & n_spectral_peaks > 0]
      if (!nrow(var_data)) {
        plot.new(); text(0.5,0.5,paste("No spectral peaks\n",var),cex=1.2); next
      }
      r_p <- create_raster_from_table(var_data, template_bc, "dominant_period")
      plot_raster_panel(r_p, paste0(var, ": Dominant Cycle (yr)"),
                        range(var_data$dominant_period, na.rm=TRUE),
                        hcl.colors(101, "plasma"), legend_title="Years")
    }
    mtext("Spectral Analysis: Dominant Periodicities",
          outer=TRUE, cex=1.1, font=2, line=0.5)
  })
  if (!is.null(result)) cat(sprintf("✅ Spectral map: %s\n", basename(pdf_path)))
  invisible(result)
}
####################################################################################
# FUNCTION 5: Spatial pattern analysis
####################################################################################
create_spatial_pattern_analysis <- function() {
  vars_here <- intersect(c("Precipitation", "PET", "PET_Thw", "Temperature"),
                         .vars_present)
  pdf_path <- file.path(out_dir, "spatial_patterns.pdf")
  result <- .pdf_safe(pdf_path, width=5*length(vars_here), height=5, expr_fn=function() {
    par(mfrow=c(1,length(vars_here)), mar=c(2,1.5,2.5,1), oma=c(1,1,3,1))
    for (var in vars_here) {
      var_data <- all_results[variable==var & period=="annual" & !is_basin_average & !filtered_vc]
      if (!nrow(var_data)) {
        plot.new(); text(0.5,0.5,paste("No data\n",var),cex=1.2); next
      }
      var_data[, spatial_cv := abs(sl_vc) / (abs(mean(sl_vc, na.rm=TRUE)) + 0.001)]
      r_cv <- create_raster_from_table(var_data, template_bc, "spatial_cv")
      plot_raster_panel(r_cv, paste0(var, ": Spatial Variability (CV)"),
                        c(0,5), hcl.colors(101, "YlOrRd"), legend_title="CV")
    }
    mtext("Spatial Variability of Annual Trends (CV)",
          outer=TRUE, cex=1.1, font=2, line=0.5)
  })
  if (!is.null(result)) cat(sprintf("✅ Spatial pattern map: %s\n", basename(pdf_path)))
  invisible(result)
}
####################################################################################
# FUNCTION 6: Method comparison (VC vs TFPW)
####################################################################################
create_method_comparison <- function() {
  comp_data <- all_results[period=="annual" & !is_basin_average &
                             !filtered_vc & !filtered_tfpw &
                             variable %in% c("Precipitation", "PET", "PET_Thw", "Temperature")]
  if (!nrow(comp_data)) { cat("⚠️ No valid pixels for method comparison.\n"); return(NULL) }
  comp_data[, vc_sig := p_value_vc < 0.05]
  comp_data[, tfpw_sig := p_value_tfpw < 0.05]
  tau_range <- range(c(comp_data$tau_vc, comp_data$tau_tfpw), na.rm=TRUE)
  var_cols <- c(Precipitation="#4575b4", PET="#d73027",
                PET_Thw="#e08214", Temperature="#1a9641")
  p1 <- ggplot(comp_data, aes(x=tau_vc, y=tau_tfpw, color=variable)) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="gray40") +
    geom_point(alpha=0.35, size=1.0) +
    facet_wrap(~variable, labeller=labeller(variable=c(
      Precipitation="Precipitation", PET="PET (PM)",
      PET_Thw="PET (Thw)", Temperature="Temperature"))) +
    scale_color_manual(values=var_cols) +
    coord_equal(xlim=tau_range, ylim=tau_range) +
    theme_bw(base_size=10) + theme(legend.position="none") +
    labs(title="A) Kendall's Tau: VC vs TFPW", x="VC Method (τ)", y="TFPW Method (τ)")
  agreement <- comp_data[, .(
    Both_Sig = sum(vc_sig & tfpw_sig),
    Only_VC = sum(vc_sig & !tfpw_sig),
    Only_TFPW = sum(!vc_sig & tfpw_sig),
    Neither = sum(!vc_sig & !tfpw_sig)
  ), by=variable]
  ag_l <- melt(agreement, id.vars="variable")
  p2 <- ggplot(ag_l, aes(x=variable, y=value, fill=variable)) +
    geom_col() + facet_wrap(~variable, scales="free_x") +
    geom_text(aes(label=sprintf("%.0f%%", value/sum(value)*100), group=variable),
              position=position_stack(vjust=0.5), size=2.8) +
    scale_fill_manual(values=var_cols) +
    theme_bw(base_size=10) + theme(legend.position="none",
                                   axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    labs(title="B) Method Agreement", x=NULL, y="N Pixels")
  combined <- p1 / p2 +
    plot_annotation(title="Mann-Kendall Method Comparison (VC vs TFPW)",
                    theme=theme(plot.title=element_text(size=12,face="bold",hjust=0.5)))
  ggsave(file.path(out_dir,"method_comparison_vc_vs_tfpw.png"),
         combined, width=14, height=10, dpi=300)
  cat("✅ Method comparison plot saved\n")
  invisible(combined)
}
####################################################################################
####################################################################################
# FUNCTION 7: Statistics table
####################################################################################
create_statistics_table <- function() {
  has_rho_vc <- "rho1_vc" %in% names(all_results)
  has_rho_tfpw <- "rho1_tfpw" %in% names(all_results)
  stats_summary <- all_results[period=="annual" & !is_basin_average, .(
    N_Pixels = .N,
    N_Valid_VC = sum(!filtered_vc),
    N_Valid_TFPW = sum(!filtered_tfpw),
    Pct_Sig_VC = {
      n <- sum(!filtered_vc)
      if (n > 0) sprintf("%.1f%%", sum(p_value_vc < 0.05 & !filtered_vc, na.rm=TRUE) / n * 100)
      else "N/A"
    },
    Pct_Sig_TFPW = {
      n <- sum(!filtered_tfpw)
      if (n > 0) sprintf("%.1f%%", sum(p_value_tfpw < 0.05 & !filtered_tfpw, na.rm=TRUE) / n * 100)
      else "N/A"
    },
    Median_Tau_VC = { v <- tau_vc[!filtered_vc]; if (length(v) > 0) sprintf("%.3f", median(v, na.rm=TRUE)) else "N/A" },
    Median_Slope_VC = { v <- sl_vc[!filtered_vc]; if (length(v) > 0) sprintf("%.4f", median(v, na.rm=TRUE)) else "N/A" },
    Mean_Autocorr = if (has_rho_vc) { v <- rho1_vc[!filtered_vc]; if (length(v) > 0) sprintf("%.3f", mean(v, na.rm=TRUE)) else "N/A" } else NA_character_
  ), by=.(variable)]
  fwrite(stats_summary, file.path(out_dir, "summary_statistics.csv"))
  cat("✅ Statistics table saved\n"); print(as.data.frame(stats_summary))
  invisible(stats_summary)
}
####################################################################################
# FUNCTION 8: Basin time series (Pr + PET_PM + PET_Thw + Temperature)
####################################################################################
create_basin_timeseries_plot <- function() {
  has_tair <- "tair_degC_month" %in% names(basin_avg_monthly)
  has_pet_thw <- "pet_thw_mm_month" %in% names(basin_avg_monthly)
  .ols_line <- function(dates, vals) {
    yr <- as.numeric(format(dates, "%Y")) +
      (as.numeric(format(dates, "%m")) - 0.5) / 12
    ok <- !is.na(yr) & !is.na(vals)
    if (sum(ok) < 10) return(NULL)
    lm(vals[ok] ~ yr[ok])
  }
  n_panels <- 2L + has_pet_thw + has_tair
  png_path <- file.path(out_dir, "basin_average_timeseries.png")
  png(png_path, width=12, height=4*n_panels, units="in", res=300)
  layout(matrix(seq_len(n_panels), nrow=n_panels), heights=rep(1, n_panels))
  pr_roll <- rollmean(basin_avg_monthly$precip_mm_month, k=12, fill=NA, align="center")
  par(mar=c(3,5,3.5,2))
  plot(basin_avg_monthly$date, basin_avg_monthly$precip_mm_month,
       type="l", col=rgb(0.2,0.4,0.8,0.3), lwd=0.8,
       xlab="", ylab="Precipitation (mm/month)",
       main="Nechako Basin Climate Variables (1950–2025)",
       ylim=c(0, max(basin_avg_monthly$precip_mm_month, na.rm=TRUE)*1.1),
       cex.lab=1.1, cex.main=1.2)
  lines(basin_avg_monthly$date, pr_roll, col="blue", lwd=2.5)
  legend("topright", legend=c("Pr (monthly)", "Pr 12-mo mean"),
         col=c(rgb(0.2,0.4,0.8,0.3), "blue"), lwd=c(0.8,2.5), bty="n", cex=0.9)
  pet_pm <- basin_avg_monthly$pet_mm_month
  ylim_p <- c(0, max(c(pet_pm, if (has_pet_thw) basin_avg_monthly$pet_thw_mm_month else 0),
                     na.rm=TRUE) * 1.15)
  plot(basin_avg_monthly$date, pet_pm,
       type="l", col=rgb(0.8,0.2,0.2,0.25), lwd=0.6,
       xlab="", ylab="PET (mm/month)", ylim=ylim_p, cex.lab=1.1,
       main="Potential Evapotranspiration — Penman-Monteith vs Thornthwaite")
  lines(basin_avg_monthly$date, rollmean(pet_pm, k=12, fill=NA, align="center"),
        col="red3", lwd=2.5)
  if (has_pet_thw) {
    pet_thw <- basin_avg_monthly$pet_thw_mm_month
    lines(basin_avg_monthly$date, pet_thw, col=rgb(0.9,0.5,0,0.25), lwd=0.6)
    lines(basin_avg_monthly$date, rollmean(pet_thw, k=12, fill=NA, align="center"),
          col="#e08214", lwd=2.5)
    legend("topright",
           legend=c("PET_PM (monthly)", "PET_PM 12-mo", "PET_Thw (monthly)", "PET_Thw 12-mo"),
           col=c(rgb(0.8,0.2,0.2,0.3), "red3",rgb(0.9,0.5,0,0.3), "#e08214"),
           lwd=c(0.6,2.5,0.6,2.5), bty="n", cex=0.85)
  } else {
    legend("topright", legend=c("PET_PM (monthly)", "PET_PM 12-mo"),
           col=c(rgb(0.8,0.2,0.2,0.3), "red3"), lwd=c(0.6,2.5), bty="n", cex=0.9)
  }
  if (has_tair) {
    tair <- basin_avg_monthly$tair_degC_month
    tair_roll <- rollmean(tair, k=12, fill=NA, align="center")
    par(mar=c(if (n_panels >3) 3 else 4, 5, 2, 2))
    plot(basin_avg_monthly$date, tair,
         type="l", col=rgb(0.6,0.2,0.6,0.25), lwd=0.6,
         xlab="", ylab="Temperature (°C)", cex.lab=1.1,
         ylim=range(tair, na.rm=TRUE) + c(-0.5,0.5),
         main="Basin-Mean 2-m Air Temperature")
    lines(basin_avg_monthly$date, tair_roll, col="purple", lwd=2.5)
    fit_t <- tryCatch({
      yr_frac <- as.numeric(format(basin_avg_monthly$date, "%Y")) +
        (as.numeric(format(basin_avg_monthly$date, "%m"))-0.5)/12
      ok <- !is.na(yr_frac) & !is.na(tair)
      if (sum(ok) < 10) NULL else lm(tair[ok] ~ yr_frac[ok])
    }, error=function(e) NULL)
    if (!is.null(fit_t)) {
      yr_seq <- seq(min(as.numeric(format(basin_avg_monthly$date, "%Y"))),
                    max(as.numeric(format(basin_avg_monthly$date, "%Y"))) + 1, length.out=200)
      pred_y <- coef(fit_t)[1] + coef(fit_t)[2] * (yr_seq + 0.5/12)
      lines(as.Date(paste0(floor(yr_seq), "-07-01")), pred_y, col="#d73027", lwd=2, lty=2)
      mtext(sprintf("OLS: %+.4f °C/yr", coef(fit_t)[2]),
            side=3, line=0, adj=1, cex=0.85, col="#d73027", font=2)
    }
    legend("topleft", legend=c("Tair (monthly)", "12-mo mean", "OLS trend"),
           col=c(rgb(0.6,0.2,0.6,0.25), "purple", "#d73027"),
           lwd=c(0.6,2.5,2), lty=c(1,1,2), bty="n", cex=0.9)
  }
  dev.off()
  cat(sprintf("✅ Basin time series: %s\n", basename(png_path)))
  invisible(png_path)
}
####################################################################################
# FUNCTION 13: CALENDAR-MONTHLY TEMPORAL CHANGE PLOTS  [NEW]
####################################################################################
.build_calmonth_df <- function(dates, monthly_vec) {
  data.frame(
    year = as.integer(format(dates, "%Y")),
    month = as.integer(format(dates, "%m")),
    month_abb = factor(month.abb[as.integer(format(dates, "%m"))],
                       levels=month.abb),
    value = as.numeric(monthly_vec),
    stringsAsFactors = FALSE
  )
}
create_calmonth_temporal_changes <- function(var_label, monthly_vec, dates,
                                             y_unit, color_line, color_ols,
                                             pdf_stem) {
  df <- .build_calmonth_df(dates, monthly_vec)
  df <- df[!is.na(df$value), ]
  if (!nrow(df)) {
    cat(sprintf(" ⚠ %s: no data — skipping\n", var_label)); return(invisible(NULL))
  }
  panels <- lapply(1:12, function(m) {
    sub <- df[df$month == m, ]
    n_obs <- nrow(sub)
    y_lbl <- if (m %in% c(1,5,9)) y_unit else NULL
    p <- ggplot2::ggplot(sub, ggplot2::aes(x=year, y=value)) +
      ggplot2::geom_point(colour=color_line, size=0.9, alpha=0.65) +
      ggplot2::geom_line(colour=color_line, linewidth=0.45, alpha=0.5)
    sig_label <- ""
    if (n_obs >= 10) {
      fit <- tryCatch(lm(value ~ year, data=sub), error=function(e) NULL)
      if (!is.null(fit)) {
        coef_sum <- summary(fit)$coefficients
        if ("year" %in% rownames(coef_sum)) {
          pval <- coef_sum["year", 4]
          slope <- coef(fit)["year"]
          pred <- data.frame(year=range(sub$year, na.rm=TRUE))
          pred$y <- predict(fit, newdata=pred)
          p <- p +
            ggplot2::geom_line(data=pred,
                               ggplot2::aes(x=year, y=y),
                               colour=color_ols, linewidth=1.1, linetype="solid",
                               inherit.aes=FALSE)
          sig_label <- if (!is.na(pval) && pval < 0.05)
            sprintf("%+.3f/yr *", slope) else sprintf("%+.3f/yr", slope)
        }
      }
    }
    p +
      ggplot2::annotate("text", x=-Inf, y=Inf, label=sig_label,
                        hjust=-0.1, vjust=1.4, size=2.5,
                        colour=ifelse(grepl("\\*", sig_label), color_ols, "grey45"),
                        fontface="bold") +
      ggplot2::scale_x_continuous(breaks=seq(1960, 2020, by=20)) +
      ggplot2::labs(title=month.abb[m], x=NULL, y=y_lbl,
                    subtitle=sprintf("n=%d", n_obs)) +
      ggplot2::theme_classic(base_size=8.5) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face="bold", size=9, hjust=0.5),
        plot.subtitle = ggplot2::element_text(size=7, colour="grey50", hjust=0.5),
        axis.title.y = ggplot2::element_text(size=8),
        axis.text = ggplot2::element_text(size=7),
        panel.grid.major = ggplot2::element_line(colour="grey94", linewidth=0.25),
        plot.margin = ggplot2::margin(2,3,2,3)
      )
  })
  fig <- patchwork::wrap_plots(panels, ncol=4) +
    patchwork::plot_annotation(
      title = sprintf("Calendar-Monthly Temporal Changes — %s — Nechako Basin", var_label),
      subtitle = paste0(
        "Each panel = one calendar month; each dot = one year's value. ",
        "Solid line = OLS trend. * p < 0.05. ",
        "Slope in ", y_unit, "/yr."),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face="bold", size=11, hjust=0.5),
        plot.subtitle = ggplot2::element_text(size=8, colour="grey35", hjust=0)))
  out_pdf <- file.path(out_dir, paste0(pdf_stem, ".pdf"))
  out_png <- file.path(out_dir, paste0(pdf_stem, ".png"))
  tryCatch(ggplot2::ggsave(out_pdf, fig, width=13, height=9, units="in", device="pdf"),
           error=function(e) cat(" ⚠ PDF: ", e$message, "\n"))
  tryCatch(ggplot2::ggsave(out_png, fig, width=13, height=9, units="in",
                           dpi=300, device="png"),
           error=function(e) cat(" ⚠ PNG: ", e$message, "\n"))
  cat(sprintf(" ✓ %s: %s / %s\n", var_label, basename(out_pdf), basename(out_png)))
  invisible(fig)
}
create_all_calmonth_temporal_changes <- function() {
  cat("\n── Function 13: Calendar-monthly temporal change plots ──\n")
  create_calmonth_temporal_changes(
    var_label = "Precipitation",
    monthly_vec = basin_avg_monthly$precip_mm_month,
    dates = basin_avg_monthly$date,
    y_unit = "mm/month",
    color_line = "#4575b4",
    color_ols = "#08306b",
    pdf_stem = "calmonth_temporal_Precipitation")
  create_calmonth_temporal_changes(
    var_label = "PET — Penman-Monteith",
    monthly_vec = basin_avg_monthly$pet_mm_month,
    dates = basin_avg_monthly$date,
    y_unit = "mm/month",
    color_line = "#d73027",
    color_ols = "#67000d",
    pdf_stem = "calmonth_temporal_PET_PM")
  if ("pet_thw_mm_month" %in% names(basin_avg_monthly)) {
    create_calmonth_temporal_changes(
      var_label = "PET — Thornthwaite",
      monthly_vec = basin_avg_monthly$pet_thw_mm_month,
      dates = basin_avg_monthly$date,
      y_unit = "mm/month",
      color_line = "#e08214",
      color_ols = "#7f2704",
      pdf_stem = "calmonth_temporal_PET_Thw")
  } else {
    cat(" ⚠ 13c: pet_thw_mm_month not found — skipped\n")
  }
  if ("tair_degC_month" %in% names(basin_avg_monthly)) {
    create_calmonth_temporal_changes(
      var_label = "Temperature (2-m)",
      monthly_vec = basin_avg_monthly$tair_degC_month,
      dates = basin_avg_monthly$date,
      y_unit = "°C",
      color_line = "#6a3d9a",
      color_ols = "#3d0072",
      pdf_stem = "calmonth_temporal_Temperature")
  } else {
    cat(" ⚠ 13d: tair_degC_month not found — skipped\n")
  }
  cat(" 13e: combined per-month slope profile...\n")
  .slope_profile_df <- function(dates, vals, var_label, unit) {
    do.call(rbind, lapply(1:12, function(m) {
      yr <- as.integer(format(dates, "%Y"))
      mo <- as.integer(format(dates, "%m"))
      sub <- data.frame(year=yr[mo==m], value=vals[mo==m])
      sub <- sub[!is.na(sub$value), ]
      if (nrow(sub) < 10) return(data.frame(
        month=m, month_abb=month.abb[m], slope=NA_real_, pval=NA_real_,
        sig=FALSE, variable=var_label, unit=unit, stringsAsFactors=FALSE))
      fit <- tryCatch(lm(value ~ year, data=sub), error=function(e) NULL)
      if (is.null(fit)) return(data.frame(
        month=m, month_abb=month.abb[m], slope=NA_real_, pval=NA_real_,
        sig=FALSE, variable=var_label, unit=unit, stringsAsFactors=FALSE))
      data.frame(
        month = m,
        month_abb = month.abb[m],
        slope = coef(fit)["year"],
        pval = summary(fit)$coefficients["year",4],
        sig = summary(fit)$coefficients["year",4] < 0.05,
        variable = var_label, unit=unit, stringsAsFactors=FALSE)
    }))
  }
  sp_list <- list(
    .slope_profile_df(basin_avg_monthly$date,
                      basin_avg_monthly$precip_mm_month, "Precipitation", "mm/month"),
    .slope_profile_df(basin_avg_monthly$date,
                      basin_avg_monthly$pet_mm_month, "PET_PM", "mm/month")
  )
  if ("pet_thw_mm_month" %in% names(basin_avg_monthly))
    sp_list[[length(sp_list)+1]] <- .slope_profile_df(
      basin_avg_monthly$date, basin_avg_monthly$pet_thw_mm_month, "PET_Thw", "mm/month")
  if ("tair_degC_month" %in% names(basin_avg_monthly))
    sp_list[[length(sp_list)+1]] <- .slope_profile_df(
      basin_avg_monthly$date, basin_avg_monthly$tair_degC_month, "Temperature", "°C")
  sp_df <- do.call(rbind, sp_list)
  sp_df$month_abb <- factor(sp_df$month_abb, levels=month.abb)
  sp_df$variable <- factor(sp_df$variable,
                           levels=c("Precipitation","PET_PM","PET_Thw","Temperature"))
  fwrite(as.data.table(sp_df),
         file.path(out_dir, "calmonth_ols_slope_profiles.csv"))
  var_pal <- c(Precipitation="#4575b4", PET_PM="#d73027",
               PET_Thw="#e08214", Temperature="#6a3d9a")
  p13e <- ggplot2::ggplot(sp_df,
                          ggplot2::aes(x=month_abb, y=slope, fill=variable,
                                       alpha=sig)) +
    ggplot2::geom_col(position="dodge", colour="white", linewidth=0.2, width=0.8) +
    ggplot2::geom_point(
      data=sp_df[sp_df$sig, ],
      ggplot2::aes(x=month_abb, y=slope + sign(slope)*max(abs(sp_df$slope),na.rm=TRUE)*0.07,
                   colour=variable, group=variable),
      shape=8, size=2.0, position=ggplot2::position_dodge(0.8),
      inherit.aes=FALSE, show.legend=FALSE) +
    ggplot2::geom_hline(yintercept=0, linewidth=0.45) +
    ggplot2::facet_wrap(~variable, ncol=1, scales="free_y",
                        labeller=ggplot2::labeller(variable=c(
                          Precipitation="Precipitation (mm/month/yr)",
                          PET_PM="PET — Penman-Monteith (mm/month/yr)",
                          PET_Thw="PET — Thornthwaite (mm/month/yr)",
                          Temperature="Temperature (°C/yr)"))) +
    ggplot2::scale_fill_manual(values=var_pal, guide="none") +
    ggplot2::scale_colour_manual(values=var_pal, guide="none") +
    ggplot2::scale_alpha_manual(values=c("TRUE"=1.0, "FALSE"=0.40), guide="none") +
    ggplot2::labs(
      title = "Per-Calendar-Month OLS Slope Profile — All Variables",
      subtitle = "Solid bars = significant (p < 0.05); faded = non-significant. ★ = p < 0.05 marker.",
      x = "Calendar month",
      y = "OLS slope (unit/yr)") +
    ggplot2::theme_classic(base_size=10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face="bold", size=12, hjust=0.5),
      strip.text = ggplot2::element_text(face="bold", size=9),
      axis.title.x = ggplot2::element_text(size=9),
      panel.grid.major.y = ggplot2::element_line(colour="grey92", linewidth=0.3),
      plot.margin = ggplot2::margin(4,10,4,4))
  out_pdf13e <- file.path(out_dir, "calmonth_slope_profile_all_vars.pdf")
  out_png13e <- file.path(out_dir, "calmonth_slope_profile_all_vars.png")
  tryCatch(ggplot2::ggsave(out_pdf13e, p13e, width=12, height=10, units="in", device="pdf"),
           error=function(e) cat(" ⚠ 13e PDF: ", e$message, "\n"))
  tryCatch(ggplot2::ggsave(out_png13e, p13e, width=12, height=10, units="in",
                           dpi=300, device="png"),
           error=function(e) cat(" ⚠ 13e PNG: ", e$message, "\n"))
  cat(sprintf(" ✓ 13e combined slope profile: %s\n", basename(out_png13e)))
  cat("── Function 13 complete.\n")
  invisible(NULL)
}
####################################################################################
# FUNCTION 14: PET_PM vs PET_Thw COMPARISON  [NEW]
####################################################################################
create_pet_comparison_plots <- function() {
  if (!has_pet_thw) {
    cat("\n⚠️ PET_Thw not present — skipping Function 14.\n")
    return(invisible(NULL))
  }
  if (!"pet_thw_mm_month" %in% names(basin_avg_monthly)) {
    cat("\n⚠️ pet_thw_mm_month not in basin_avg_monthly — skipping Function 14.\n")
    return(invisible(NULL))
  }
  cat("\n── Function 14: PET_PM vs PET_Thw comparison plots ──\n")
  pet_pm <- basin_avg_monthly$pet_mm_month
  pet_thw <- basin_avg_monthly$pet_thw_mm_month
  dates <- basin_avg_monthly$date
  yr <- as.integer(format(dates, "%Y"))
  mo <- as.integer(format(dates, "%m"))
  theme_cmp <- ggplot2::theme_classic(base_size=10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face="bold", size=11, hjust=0),
      plot.subtitle = ggplot2::element_text(size=8.5, colour="grey40", hjust=0),
      panel.grid.major = ggplot2::element_line(colour="grey93", linewidth=0.3),
      plot.margin = ggplot2::margin(4,8,4,4))
  cat(" 14a: time series overlay...\n")
  df_ts <- data.frame(
    date = dates,
    PET_PM = pet_pm,
    PET_Thw = pet_thw,
    pm_roll = rollmean(pet_pm, k=12, fill=NA, align="center"),
    thw_roll = rollmean(pet_thw, k=12, fill=NA, align="center"))
  p14a <- ggplot2::ggplot(df_ts, ggplot2::aes(x=date)) +
    ggplot2::geom_line(ggplot2::aes(y=PET_PM), colour=rgb(0.84,0.19,0.15,0.25), linewidth=0.45) +
    ggplot2::geom_line(ggplot2::aes(y=PET_Thw), colour=rgb(0.88,0.51,0.08,0.25), linewidth=0.45) +
    ggplot2::geom_line(ggplot2::aes(y=pm_roll), colour="#d73027", linewidth=1.2, na.rm=TRUE) +
    ggplot2::geom_line(ggplot2::aes(y=thw_roll), colour="#e08214", linewidth=1.2, linetype="dashed", na.rm=TRUE) +
    ggplot2::scale_x_date(date_breaks="10 years", date_labels="%Y",
                          expand=ggplot2::expansion(mult=c(0.01,0.01))) +
    ggplot2::labs(title="Basin-Mean PET — Penman-Monteith vs Thornthwaite (1950–2025)",
                  subtitle="Thick lines = 12-month rolling mean. Red = PM, Orange-dashed = Thw.",
                  x="Year", y="PET (mm/month)") +
    theme_cmp
  tryCatch(ggplot2::ggsave(file.path(out_dir, "pet_comparison_14a_timeseries.pdf"),
                           p14a, width=12, height=5, units="in", device="pdf"),
           error=function(e) cat(" ⚠ 14a PDF: ", e$message, "\n"))
  tryCatch(ggplot2::ggsave(file.path(out_dir, "pet_comparison_14a_timeseries.png"),
                           p14a, width=12, height=5, units="in", dpi=300, device="png"),
           error=function(e) cat(" ⚠ 14a PNG: ", e$message, "\n"))
  cat(" ✓ 14a saved\n")
  cat(" 14b: per-month slope comparison...\n")
  slopes_b <- do.call(rbind, lapply(1:12, function(m) {
    idx <- mo == m
    .fit <- function(x) {
      sub <- data.frame(yr=yr[idx], v=x[idx]); sub <- sub[!is.na(sub$v),]
      if (nrow(sub) <10) return(c(NA,NA))
      f <- tryCatch(lm(v~yr, data=sub), error=function(e) NULL)
      if (is.null(f)) return(c(NA,NA))
      c(coef(f)["yr"], summary(f)$coefficients["yr",4])
    }
    spm <- .fit(pet_pm); sth <- .fit(pet_thw)
    data.frame(month=m, month_abb=month.abb[m],
               slope_PM=spm[1], pval_PM=spm[2], sig_PM=!is.na(spm[2]) & spm[2] <0.05,
               slope_Thw=sth[1], pval_Thw=sth[2], sig_Thw=!is.na(sth[2]) & sth[2] <0.05)
  }))
  slopes_b$month_abb <- factor(slopes_b$month_abb, levels=month.abb)
  slopes_long <- data.table::melt(as.data.table(slopes_b),
                                  id.vars=c("month", "month_abb"),
                                  measure.vars=c("slope_PM", "slope_Thw"),
                                  variable.name="method", value.name="slope")
  slopes_long[, sig := ifelse(method=="slope_PM",
                              slopes_b$sig_PM[match(month, slopes_b$month)],
                              slopes_b$sig_Thw[match(month, slopes_b$month)])]
  slopes_long[, method_lbl := ifelse(method=="slope_PM", "Penman-Monteith", "Thornthwaite")]
  p14b <- ggplot2::ggplot(slopes_long,
                          ggplot2::aes(x=month_abb, y=slope, fill=method_lbl,
                                       alpha=sig)) +
    ggplot2::geom_col(position="dodge", colour="white", linewidth=0.2, width=0.75) +
    ggplot2::geom_hline(yintercept=0, linewidth=0.5) +
    ggplot2::scale_fill_manual(values=c("Penman-Monteith"="#d73027",
                                        "Thornthwaite"="#e08214"),
                               name=NULL) +
    ggplot2::scale_alpha_manual(values=c("TRUE"=1.0, "FALSE"=0.35), guide="none") +
    ggplot2::labs(
      title = "Per-Month OLS Slope: PET_PM vs PET_Thw",
      subtitle = "Solid bars = p < 0.05; faded = non-significant. (mm/month per year)",
      x="Calendar month", y="OLS slope (mm/month/yr)") +
    theme_cmp + ggplot2::theme(legend.position="top")
  tryCatch(ggplot2::ggsave(file.path(out_dir, "pet_comparison_14b_monthly_slopes.pdf"),
                           p14b, width=11, height=5, units="in", device="pdf"),
           error=function(e) cat(" ⚠ 14b PDF: ", e$message, "\n"))
  tryCatch(ggplot2::ggsave(file.path(out_dir, "pet_comparison_14b_monthly_slopes.png"),
                           p14b, width=11, height=5, units="in", dpi=300, device="png"),
           error=function(e) cat(" ⚠ 14b PNG: ", e$message, "\n"))
  cat(" ✓ 14b saved\n")
  cat(" 14c: basin trend maps (Sen's slope, PM vs Thw)...\n")
  result <- .pdf_safe(file.path(out_dir, "pet_comparison_14c_trendsmap.pdf"),
                      width=10, height=5, expr_fn=function() {
                        par(mfrow=c(1,2), mar=c(2,1.5,2.5,1), oma=c(1,1,3,1))
                        if (!is.null(r_pet))
                          plot_raster_panel(r_pet$mag_vc, "PET_PM: Sen's Slope (mm/yr)", c(-0.5,3),
                                            hcl.colors(101, "YlOrRd"), legend_title="mm/yr")
                        if (!is.null(r_pet_thw))
                          plot_raster_panel(r_pet_thw$mag_vc, "PET_Thw: Sen's Slope (mm/yr)", c(-0.5,3),
                                            hcl.colors(101, "YlOrRd"), legend_title="mm/yr")
                        mtext("PET Trend Comparison — VC-MK Sen's Slope | Penman-Monteith vs Thornthwaite",
                              outer=TRUE, cex=1.1, font=2, line=0.5)
                      })
  if (!is.null(result)) cat(" ✓ 14c saved\n")
  cat(" 14d: monthly PET bias (PM − Thw)...\n")
  bias_df <- do.call(rbind, lapply(1:12, function(m) {
    idx <- mo == m
    mu_pm <- mean(pet_pm[idx], na.rm=TRUE)
    mu_thw <- mean(pet_thw[idx], na.rm=TRUE)
    data.frame(month=m, month_abb=month.abb[m],
               bias=mu_pm - mu_thw, mu_pm=mu_pm, mu_thw=mu_thw)
  }))
  bias_df$month_abb <- factor(bias_df$month_abb, levels=month.abb)
  bias_df$sign_col <- ifelse(bias_df$bias > 0, "PM > Thw", "Thw > PM")
  p14d <- ggplot2::ggplot(bias_df,
                          ggplot2::aes(x=month_abb, y=bias, fill=sign_col)) +
    ggplot2::geom_col(colour="white", linewidth=0.25) +
    ggplot2::geom_hline(yintercept=0, linewidth=0.5) +
    ggplot2::scale_fill_manual(values=c("PM > Thw"="#d73027", "Thw > PM"="#e08214"),
                               name=NULL) +
    ggplot2::labs(
      title = "Monthly Mean PET Bias: Penman-Monteith minus Thornthwaite",
      subtitle = "Positive = PM exceeds Thw (radiation+wind component active)",
      x="Calendar month", y="PET bias (mm/month)") +
    theme_cmp + ggplot2::theme(legend.position="top")
  tryCatch(ggplot2::ggsave(file.path(out_dir, "pet_comparison_14d_monthly_bias.pdf"),
                           p14d, width=10, height=5, units="in", device="pdf"),
           error=function(e) cat(" ⚠ 14d PDF: ", e$message, "\n"))
  tryCatch(ggplot2::ggsave(file.path(out_dir, "pet_comparison_14d_monthly_bias.png"),
                           p14d, width=10, height=5, units="in", dpi=300, device="png"),
           error=function(e) cat(" ⚠ 14d PNG: ", e$message, "\n"))
  cat(" ✓ 14d saved\n")
  cat(" 14e: per-month scatter PM vs Thw...\n")
  sc_df <- data.frame(date=dates, pm=pet_pm, thw=pet_thw,
                      month_abb=factor(month.abb[mo], levels=month.abb),
                      decade=factor(paste0(floor(yr/10)*10, "s")))
  sc_df <- sc_df[!is.na(sc_df$pm) & !is.na(sc_df$thw), ]
  p14e <- ggplot2::ggplot(sc_df,
                          ggplot2::aes(x=thw, y=pm, colour=decade)) +
    ggplot2::geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey40") +
    ggplot2::geom_point(size=0.8, alpha=0.5) +
    ggplot2::facet_wrap(~month_abb, ncol=4, scales="free") +
    ggplot2::scale_colour_brewer(palette="RdYlBu", direction=-1, name="Decade") +
    ggplot2::labs(
      title = "Monthly PET Scatter: Penman-Monteith (y) vs Thornthwaite (x)",
      subtitle = "Dashed 1:1 line. Colour = decade. Divergence from 1:1 = radiation/wind component.",
      x="PET Thornthwaite (mm/month)", y="PET Penman-Monteith (mm/month)") +
    theme_cmp +
    ggplot2::theme(strip.text=ggplot2::element_text(face="bold", size=8),
                   legend.position="right")
  tryCatch(ggplot2::ggsave(file.path(out_dir, "pet_comparison_14e_scatter.pdf"),
                           p14e, width=13, height=9, units="in", device="pdf"),
           error=function(e) cat(" ⚠ 14e PDF: ", e$message, "\n"))
  tryCatch(ggplot2::ggsave(file.path(out_dir, "pet_comparison_14e_scatter.png"),
                           p14e, width=13, height=9, units="in", dpi=300, device="png"),
           error=function(e) cat(" ⚠ 14e PNG: ", e$message, "\n"))
  cat(" ✓ 14e saved\n")
  cat("── Function 14 complete.\n")
  invisible(NULL)
}
################################################################################
################################################################################
# FUNCTION 15: PET-PM SEASONAL TRENDS PLOT
################################################################################
create_pet_seasonal_trends <- function(
    df = basin_avg_monthly,
    out_dir_loc = out_dir,  # Renamed to avoid recursive promise error
    basin_label = "Nechako River Basin",
    year_start  = 1950L,
    year_end    = 2025L
) {
  cat("\n── Function 15: PET-PM seasonal trends plot ──\n")
  if (is.null(df) || !"date" %in% names(df)) {
    cat("  ⚠ Basin monthly data not available — skipping PET seasonal trends.\n")
    return(invisible(NULL))
  }
  
  pet_raw <- data.table::as.data.table(df) 
  pet_raw <- pet_raw[, .(date, pet_mm_month)]
  setnames(pet_raw, c("date", "pet"))
  pet_raw <- pet_raw[!is.na(pet)]
  
  pet_raw[, .date := as.Date(date)]
  pet_raw[, year  := as.integer(format(.date, "%Y"))]
  pet_raw[, month := as.integer(format(.date, "%m"))]
  pet_raw <- pet_raw[year >= year_start & year <= year_end]
  
  pet_raw[, season := data.table::fcase(
    month %in% c(12, 1, 2), "DJF (Dec\u2013Feb)",
    month %in% c(3, 4, 5),  "MAM (Mar\u2013May)",
    month %in% c(6, 7, 8),  "JJA (Jun\u2013Aug)",
    month %in% c(9, 10, 11),"SON (Sep\u2013Nov)"
  )]
  pet_raw[, season_year := data.table::fifelse(month == 12L, year + 1L, year)]
  
  SEASON_LEVELS <- c("DJF (Dec\u2013Feb)", "MAM (Mar\u2013May)",
                     "JJA (Jun\u2013Aug)", "SON (Sep\u2013Nov)")
  
  seas_means <- pet_raw[, .(pet_mean = mean(pet, na.rm = TRUE)), by = .(season_year, season)]
  seas_means <- seas_means[season_year >= year_start & season_year <= year_end]
  seas_means[, season := factor(season, levels = SEASON_LEVELS)]
  
  trend_stats <- seas_means[, {
    fit   <- lm(pet_mean ~ season_year)
    slope <- coef(fit)[["season_year"]]
    pval  <- summary(fit)$coefficients["season_year", "Pr(>|t|)"]
    list(slope = slope, p_value = pval)
  }, by = season]
  
  trend_label <- paste(
    sapply(SEASON_LEVELS, function(s) {
      row   <- trend_stats[season == s]
      sig   <- if (!is.na(row$p_value) && row$p_value < 0.05) " *" else ""
      sprintf("%s: %+.4f mm/yr%s", sub(" \\(.*", "", s), row$slope, sig)
    }),
    collapse = "   |   "
  )
  
  season_colours <- c(
    "DJF (Dec\u2013Feb)" = "#4472C4",
    "MAM (Mar\u2013May)" = "#70AD47",
    "JJA (Jun\u2013Aug)" = "#C00000",
    "SON (Sep\u2013Nov)" = "#ED7D31"
  )
  
  p_pet <- ggplot(seas_means, aes(x = season_year, y = pet_mean, colour = season, group = season)) +
    geom_line(linewidth = 0.6, alpha = 0.85) +
    geom_point(size = 1.2, alpha = 0.75) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linetype = "dashed", linewidth = 1.6) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50", linewidth = 0.4) +
    scale_colour_manual(values = season_colours, name = NULL,
                        guide = guide_legend(direction = "horizontal", nrow = 1,
                                             override.aes = list(linewidth = 1.4, size = 2))) +
    scale_x_continuous(breaks = seq(1950, 2025, by = 10), expand = expansion(mult = 0.01)) +
    labs(
      title = sprintf("Seasonal PET-PM Trends \u2014 %s (%d\u2013%d)", basin_label, year_start, year_end),
      subtitle = paste0("Annual seasonal means with OLS trend lines (dashed).  * p < 0.05\n", trend_label),
      x = "Year", y = "Seasonal mean PET-PM (mm/month)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13, margin = margin(b = 2)),
      plot.subtitle = element_text(size = 9.5, colour = "grey30", margin = margin(b = 8)),
      legend.position = "top", legend.text = element_text(size = 10),
      legend.key.width = unit(1.6, "cm"), panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
      axis.title = element_text(size = 11), axis.text = element_text(size = 10)
    )
  
  out_png <- file.path(out_dir_loc, "FigS_PET_Seasonal_Trends.png")
  out_pdf <- file.path(out_dir_loc, "FigS_PET_Seasonal_Trends.pdf")
  tryCatch({ ggsave(out_png, p_pet, width = 14, height = 7, units = "in", dpi = 300)
    cat(sprintf("  \u2713 PET seasonal trends saved: %s\n", basename(out_png))) },
    error = function(e) cat(sprintf("  \u26a0 PET PNG: %s\n", e$message)))
  tryCatch({ ggsave(out_pdf, p_pet, width = 14, height = 7, units = "in", device = "pdf")
    cat(sprintf("  \u2713 PET seasonal trends saved: %s\n", basename(out_pdf))) },
    error = function(e) cat(sprintf("  \u26a0 PET PDF: %s\n", e$message)))
  invisible(p_pet)
}
####################################################################################
# FUNCTION: Seasonal Climatology P vs PET (Bar Chart Style - Fig3.png)
####################################################################################
create_seasonality_p_pet_plot <- function() {
  cat("\n── Creating Seasonality P vs PET plot (FigS_Seasonality_P_PET.png) ──\n")
  
  # 1. Data Preparation
  df <- data.frame(
    date = basin_avg_monthly$date,
    P = basin_avg_monthly$precip_mm_month,
    PET_PM = basin_avg_monthly$pet_mm_month
  )
  df$Month <- factor(month.abb[as.integer(format(df$date, "%m"))], levels=month.abb)
  
  # Calculate Statistics
  stats_p <- aggregate(P ~ Month, data=df, FUN=function(x) c(mean=mean(x, na.rm=TRUE), sd=sd(x, na.rm=TRUE)))
  stats_pet <- aggregate(PET_PM ~ Month, data=df, FUN=function(x) c(mean=mean(x, na.rm=TRUE), sd=sd(x, na.rm=TRUE)))
  
  # Combine into a long format for plotting
  clim_df <- data.frame(
    Month = stats_p$Month,
    P_Mean = stats_p$P[, "mean"],
    P_SD = stats_p$P[, "sd"],
    PET_Mean = stats_pet$PET_PM[, "mean"],
    PET_SD = stats_pet$PET_PM[, "sd"],
    stringsAsFactors = FALSE
  )
  clim_df$Month <- factor(clim_df$Month, levels=month.abb)
  
  # Water Balance
  clim_df$Balance <- clim_df$P_Mean - clim_df$PET_Mean
  
  # Annual totals for subtitle
  ann_P <- sum(clim_df$P_Mean)
  ann_PET <- sum(clim_df$PET_Mean)
  subtitle_txt <- sprintf("Each bar = mean of 76 years (1950–2025). Annual totals: P = %.0f mm; PET(PM) = %.0f mm; aridity index (PET/P) = %.2f.",
                          ann_P, ann_PET, ann_PET/ann_P)
  
  # Long format for Panel (a)
  plot_data_a <- data.frame(
    Month = rep(clim_df$Month, 2),
    Variable = rep(c("Precipitation (P)", "PET(PM)"), each=nrow(clim_df)),
    Mean = c(clim_df$P_Mean, clim_df$PET_Mean),
    SD = c(clim_df$P_SD, clim_df$PET_SD)
  )
  plot_data_a$Variable <- factor(plot_data_a$Variable, levels=c("Precipitation (P)", "PET(PM)"))
  
  # 2. Panel (a)
  p_a <- ggplot(plot_data_a, aes(x=Month, y=Mean, fill=Variable)) +
    geom_col(position=position_dodge(0.8), width=0.7, color="white", linewidth=0.2) +
    geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD),
                  position=position_dodge(0.8), width=0.2, linewidth=0.6) +
    scale_fill_manual(values=c("Precipitation (P)"="#4575b4", "PET(PM)"="#d73027")) +
    labs(y="mm month⁻¹",
         title="(a) Mean monthly precipitation and PET(PM)",
         subtitle="76-year mean (1950–2025) | Error bars = ±1 SD") +
    theme_classic(base_size=12) +
    theme(
      plot.title = element_text(face="bold", size=11, hjust=0),
      plot.subtitle = element_text(size=8, color="grey40", hjust=0),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle=0, hjust=0.5),
      legend.position = "top",
      legend.justification = "right",
      panel.grid.major.y = element_line(color="grey90")
    )
  
  # 3. Panel (b)
  # Color logic for balance
  clim_df$BalanceCol <- ifelse(clim_df$Balance > 0, "#4575b4", "#d73027")
  
  p_b <- ggplot(clim_df, aes(x=Month, y=Balance, fill=BalanceCol)) +
    geom_col(width=0.7, color="white", linewidth=0.2, show.legend=FALSE) +
    geom_hline(yintercept=0, linewidth=0.5, color="grey50") +
    scale_fill_identity() +
    labs(x="Calendar month", y="P − PET(PM) (mm month⁻¹)",
         title="(b) Climatic water balance P − PET(PM)",
         subtitle="Blue = monthly surplus (P > PET); red = monthly deficit (PET > P)") +
    theme_classic(base_size=12) +
    theme(
      plot.title = element_text(face="bold", size=11, hjust=0),
      plot.subtitle = element_text(size=8, color="grey40", hjust=0),
      axis.text.x = element_text(angle=0, hjust=0.5),
      panel.grid.major.y = element_line(color="grey90")
    )
  
  # Annotation for P surplus (if applicable)
  p_b <- p_b + annotate("text", x=6.5, y=max(clim_df$Balance, na.rm=TRUE)*0.8,
                        label="P surplus", color="#4575b4", size=4, fontface="italic")
  
  # 4. Combine and Save
  combined <- (p_a / p_b) +
    plot_annotation(
      title = "Nechako River Basin — Seasonal cycle of precipitation and PET(PM)",
      subtitle = subtitle_txt,
      theme = theme(plot.title = element_text(face="bold", size=14, hjust=0),
                    plot.subtitle = element_text(size=10, color="grey40", hjust=0))
    )
  
  out_path <- file.path(out_dir, "FigS_Seasonality_P_PET.png")
  ggsave(out_path, combined, width=10, height=9, dpi=300)
  cat("✓ Saved:", out_path, "\n")
}
####################################################################################
# FUNCTION 9–11: Specific-point plots (updated for 4 variables)
####################################################################################
create_point_timeseries_plots <- function() {
  if (is.null(point_monthly_ts)) {
    cat("⚠️ point_monthly_timeseries.csv not available — skipping.\n"); return(NULL)
  }
  cat("📈 Creating specific-point time series plots...\n")
  point_monthly_ts[, date := as.Date(date)]
  for (pid in unique(point_monthly_ts$point_id)) {
    sub <- point_monthly_ts[point_id == pid]
    lon <- round(sub$lon_wgs84[1], 4); lat <- round(sub$lat_wgs84[1], 4)
    title_str <- sprintf("Point #%d (%.4f°N, %.4f°W)", pid, lat, abs(lon))
    make_var_panel <- function(var_name, y_label, color_line, color_smooth) {
      df <- sub[variable == var_name]
      if (!nrow(df)) return(ggplot() +
                              annotate("text",x=0.5,y=0.5,label=paste(var_name, "not available")) +
                              theme_void())
      df <- df[order(date)]
      df[, smooth := rollmean(value, k=min(12,nrow(df)), fill=NA, align="center")]
      ggplot(df, aes(x=date, y=value)) +
        geom_line(color=color_line, alpha=0.45, linewidth=0.5) +
        geom_line(aes(y=smooth), color=color_smooth, linewidth=1.1, na.rm=TRUE) +
        labs(title=var_name, x=NULL, y=y_label) +
        theme_minimal(base_size=10) +
        theme(plot.title=element_text(face="bold"), panel.grid.minor=element_blank(),
              axis.text.x=element_text(angle=30, hjust=1))
    }
    panels <- list(
      make_var_panel("Precipitation", "mm/month", "steelblue", "blue"),
      make_var_panel("PET", "mm/month", "tomato", "red3"),
      make_var_panel("PET_Thw", "mm/month", "orange", "#e08214"),
      make_var_panel("Temperature", "°C", "orchid", "purple4")
    )
    combined <- (panels[[1]] / panels[[2]] / panels[[3]] / panels[[4]]) +
      plot_annotation(
        title = title_str,
        subtitle = "ERA5-Land monthly | coloured line = 12-month rolling mean",
        theme = theme(
          plot.title = element_text(size=12, face="bold", hjust=0.5),
          plot.subtitle = element_text(size=9, colour="grey30", hjust=0.5)))
    out_path <- file.path(out_dir,
                          sprintf("point_%d_%.4fN_%.4fW_ts_pr_pet_temp.png",
                                  pid, lat, abs(lon)))
    ggsave(out_path, combined, width=12, height=13, dpi=150)
    cat(sprintf(" Saved: %s\n", basename(out_path)))
  }
  invisible(NULL)
}
create_point_monthly_trend_plots <- function() {
  if (is.null(point_trend_stats)) {
    cat("⚠️ point_trend_stats.csv not available — skipping.\n"); return(NULL)
  }
  cat("📊 Creating specific-point monthly trend bar charts...\n")
  var_units <- c(Precipitation="Sen slope (mm/yr per month)",
                 PET="Sen slope (mm/yr per month)",
                 PET_Thw="Sen slope (mm/yr per month)",
                 Temperature="Sen slope (°C/yr)")
  var_colors <- c(Precipitation="#4575b4", PET="#d73027",
                  PET_Thw="#e08214", Temperature="#7b2d8b")
  for (pid in unique(point_trend_stats$point_id)) {
    sub <- point_trend_stats[point_id == pid & period == "monthly"]
    if (!nrow(sub)) next
    lon <- round(sub$lon_wgs84[1], 4); lat <- round(sub$lat_wgs84[1], 4)
    sub[, month_abb := factor(month_abb, levels=month.abb)]
    sub[, sig_label := ifelse(!is.na(sig_vc) & sig_vc, "*", "")]
    sub[, label_y := sl_vc + ifelse(!is.na(sl_vc) & sl_vc >=0,
                                    0.03*max(abs(sl_vc),na.rm=TRUE),
                                    -0.03*max(abs(sl_vc),na.rm=TRUE))]
    panels <- lapply(c("Precipitation", "PET", "PET_Thw", "Temperature"), function(var) {
      d <- sub[variable == var]
      if (!nrow(d)) return(ggplot()+theme_void()+labs(title=paste(var, "— no data")))
      ggplot(d, aes(x=month_abb, y=sl_vc)) +
        geom_bar(stat="identity", colour="black", linewidth=0.3, width=0.7,
                 fill=var_colors[var]) +
        geom_hline(yintercept=0, linewidth=0.5) +
        geom_text(aes(y=label_y, label=sig_label), size=4.5, fontface="bold") +
        scale_x_discrete(labels=substr(month.abb,1,1)) +
        labs(title=var, x="Month", y=var_units[var]) +
        theme_classic(base_size=10) +
        theme(plot.title=element_text(face="bold",hjust=0.5))
    })
    combined <- (panels[[1]] | panels[[2]] | panels[[3]] | panels[[4]]) +
      plot_annotation(
        title = sprintf("Per-Month Sen's Slope — Point #%d (%.4f°N, %.4f°W)",
                        pid, lat, abs(lon)),
        subtitle = "* = significant by VC Mann-Kendall (p < 0.05)",
        theme = theme(
          plot.title = element_text(size=12, face="bold", hjust=0.5),
          plot.subtitle = element_text(size=9, colour="grey30", hjust=0.5)))
    out_path <- file.path(out_dir,
                          sprintf("point_%d_%.4fN_%.4fW_monthly_slopes.png",
                                  pid, lat, abs(lon)))
    ggsave(out_path, combined, width=16, height=5, dpi=150)
    cat(sprintf(" Saved: %s\n", basename(out_path)))
  }
  invisible(NULL)
}
create_point_summary_table <- function() {
  if (is.null(point_trend_stats)) {
    cat("⚠️ point_trend_stats.csv not available.\n"); return(NULL)
  }
  cat("\n══════════════════════════════════════════════════\n")
  cat(" SPECIFIC POINT TREND STATISTICS SUMMARY\n")
  cat("══════════════════════════════════════════════════\n")
  annual <- point_trend_stats[period == "annual"]
  if (nrow(annual)) {
    cat("Annual trends:\n")
    print(as.data.frame(annual[, .(
      point_id, lat_wgs84, lon_wgs84, variable,
      n_obs, year_min, year_max,
      tau_vc=round(tau_vc,3), p_vc=round(p_value_vc,4), sig_vc,
      slope_vc=round(sl_vc,5),
      tau_tfpw=round(tau_tfpw,3), p_tfpw=round(p_value_tfpw,4), sig_tfpw,
      changepoint_yr=first_changepoint_year
    )]))
  }
  sig_monthly <- point_trend_stats[period=="monthly" & !is.na(sig_vc) & sig_vc==TRUE]
  if (nrow(sig_monthly)) {
    cat("\nSignificant monthly trends (VC-MK p < 0.05):\n")
    print(as.data.frame(sig_monthly[, .(
      point_id, variable, month_abb,
      tau_vc=round(tau_vc,3), p_vc=round(p_value_vc,4),
      slope_vc=round(sl_vc,5), sig_tfpw
    )]))
  } else {
    cat("\n (No significant monthly trends at any specific point)\n")
  }
  invisible(annual)
}
####################################################################################
# FUNCTION 12: Temperature dedicated analysis (unchanged from original)
####################################################################################
create_temperature_dedicated_plots <- function() {
  if (!"tair_degC_month" %in% names(basin_avg_monthly)) {
    cat("\u26a0\ufe0f basin_avg_monthly$tair_degC_month not found — skipping Function 12.\n")
    return(invisible(NULL))
  }
  tair <- basin_avg_monthly$tair_degC_month
  dates <- basin_avg_monthly$date
  yr <- as.integer(format(dates, "%Y"))
  mo <- as.integer(format(dates, "%m"))
  yr_frac <- yr + (mo - 0.5) / 12
  .ols <- function(x, y) {
    keep <- !is.na(x) & !is.na(y); if (sum(keep) < 5) return(NULL)
    lm(y[keep] ~ x[keep])
  }
  theme_tair <- ggplot2::theme_classic(base_size=10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face="bold", size=11, hjust=0),
      plot.subtitle = ggplot2::element_text(size=8.5, colour="grey40", hjust=0),
      panel.grid.major = ggplot2::element_line(colour="grey94", linewidth=0.3),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(4, 8, 4, 4))
  cat("\n\u2500\u2500 Function 12: Temperature dedicated plots \u2500\u2500\n")
  tair_roll <- zoo::rollmean(tair, k=12, fill=NA, align="center")
  fit_all <- .ols(yr_frac, tair)
  pred_df <- NULL; ann_slope <- "trend not available"
  if (!is.null(fit_all)) {
    pred_df <- data.frame(date=dates[!is.na(tair)],
                          tfit=coef(fit_all)[1]+coef(fit_all)[2]*yr_frac[!is.na(tair)])
    ann_slope <- sprintf("%+.4f \u00b0C/yr (%+.3f \u00b0C/decade)",
                         coef(fit_all)[2], coef(fit_all)[2]*10)
  }
  ann_df <- NULL
  if (!is.null(basin_avg_annual) && "tair_degC_year" %in% names(basin_avg_annual))
    ann_df <- data.frame(date=as.Date(paste0(basin_avg_annual$year, "-07-01")),
                         tair=basin_avg_annual$tair_degC_year)
  df_ts <- data.frame(date=dates, tair=tair, troll=tair_roll)
  p12a <- ggplot2::ggplot(df_ts, ggplot2::aes(x=date)) +
    ggplot2::geom_line(ggplot2::aes(y=tair), colour=rgb(0.6,0.2,0.6,0.3), linewidth=0.45, na.rm=TRUE) +
    ggplot2::geom_line(ggplot2::aes(y=troll), colour="purple", linewidth=1.6, na.rm=TRUE)
  if (!is.null(ann_df))
    p12a <- p12a + ggplot2::geom_point(data=ann_df, ggplot2::aes(x=date,y=tair),
                                       colour=rgb(0.4,0,0.4,0.8), size=1.4, shape=19)
  if (!is.null(pred_df))
    p12a <- p12a + ggplot2::geom_line(data=pred_df, ggplot2::aes(x=date,y=tfit),
                                      colour="#d73027", linewidth=1.2, linetype="dashed")
  p12a <- p12a +
    ggplot2::geom_hline(yintercept=0, colour="grey50", linetype="dotted", linewidth=0.5) +
    ggplot2::scale_x_date(date_breaks="10 years", date_labels="%Y",
                          expand=ggplot2::expansion(mult=c(0.01,0.01))) +
    ggplot2::labs(title="Basin-Mean 2-m Air Temperature — Nechako River Basin (1950\u20132025)",
                  subtitle=paste0("Monthly (faded), 12-mo mean (purple), annual mean (dots), OLS (red dashed)\n",
                                  "OLS slope: ", ann_slope),
                  x="Year", y="Temperature (\u00b0C)") + theme_tair
  tryCatch(ggplot2::ggsave(file.path(out_dir, "temperature_basin_timeseries.pdf"),
                           p12a, width=10, height=5, units="in", device="pdf"),
           error=function(e) cat(" \u26a0 12a PDF: ", e$message, "\n"))
  tryCatch(ggplot2::ggsave(file.path(out_dir, "temperature_basin_timeseries.png"),
                           p12a, width=10, height=5, units="in", dpi=300, device="png"),
           error=function(e) cat(" \u26a0 12a PNG: ", e$message, "\n"))
  cat(" \u2713 12a saved\n")
  ann_tair <- tapply(tair, yr, mean, na.rm=TRUE)
  ann_years <- as.integer(names(ann_tair)); ann_tair <- as.numeric(ann_tair)
  clim_idx <- ann_years >=1991 & ann_years <=2020
  clim_mean <- if (sum(clim_idx) >0) mean(ann_tair[clim_idx], na.rm=TRUE) else mean(ann_tair, na.rm=TRUE)
  df_anom <- data.frame(year=ann_years, anomaly=ann_tair - clim_mean)
  df_anom$sign <- ifelse(df_anom$anomaly >=0, "Warm", "Cool")
  df_anom$run10 <- zoo::rollmean(df_anom$anomaly, k=10, fill=NA, align="center")
  fit_anom <- .ols(df_anom$year, df_anom$anomaly)
  anom_str <- if (!is.null(fit_anom)) sprintf("OLS: %+.4f \u00b0C/yr", coef(fit_anom)[2]) else ""
  p12b <- ggplot2::ggplot(df_anom, ggplot2::aes(x=year, y=anomaly, fill=sign)) +
    ggplot2::geom_bar(stat="identity", width=0.85, colour="white", linewidth=0.1) +
    #ggplot2::geom_line(ggplot2::aes(y=run10), colour="#1f78b4", linewidth=1.2, linetype="dashed", na.rm=TRUE) +
    ggplot2::scale_fill_manual(values=c(Warm="#d73027", Cool="#4575b4"), name=NULL) +
    ggplot2::geom_hline(yintercept=0, colour="grey50", linetype="dotted", linewidth=0.5) +
    ggplot2::scale_x_continuous(breaks=seq(1950,2025,by=10)) +
    ggplot2::labs(title="Annual Temperature Anomaly — Nechako River Basin",
                  subtitle=sprintf("Baseline: WMO 1991\u20132020 (%.2f°C). Blue dashed = 10-yr mean. %s",
                                   clim_mean, anom_str),
                  x="Year", y="Temperature anomaly (\u00b0C)") + theme_tair +
    ggplot2::theme(legend.position="top")
  p12b <- p12b 
  # + ggplot2::annotate("text", x=min(df_anom$year)+2, y=max(df_anom$anomaly)*0.85,
                                   # label="10-yr rolling mean", colour="#1f78b4", hjust=0, size=3, fontface="italic")
  tryCatch(ggplot2::ggsave(file.path(out_dir, "temperature_annual_anomaly.pdf"),
                           p12b, width=10, height=5, units="in", device="pdf"),
           error=function(e) cat(" \u26a0 12b PDF: ", e$message, "\n"))
  tryCatch(ggplot2::ggsave(file.path(out_dir, "temperature_annual_anomaly.png"),
                           p12b, width=10, height=5, units="in", dpi=300, device="png"),
           error=function(e) cat(" \u26a0 12b PNG: ", e$message, "\n"))
  cat(" \u2713 12b saved\n")
  cat("\u2500\u2500 Function 12 complete (12a, 12b).\n")
  invisible(NULL)
}
####################################################################################
# ── GENERATE ALL OUTPUTS
####################################################################################
cat("\n========================================\n")
cat("\U0001f3a8 CREATING VISUALIZATIONS\n")
cat("========================================\n\n")
for (var in c("Precipitation","PET","PET_Thw","Temperature")) {
  for (meth in c("vc","tfpw")) create_publication_figure(var, meth)
}
create_comparison_figure()
create_regime_shift_maps()
create_spectral_analysis_maps()
create_spatial_pattern_analysis()
create_method_comparison()
create_statistics_table()
create_basin_timeseries_plot()
create_all_calmonth_temporal_changes()
create_pet_comparison_plots()
create_pet_seasonal_trends()
create_seasonality_p_pet_plot()
create_temperature_dedicated_plots()
create_point_timeseries_plots()
create_point_monthly_trend_plots()
create_point_summary_table()
cat("\n========================================\n")
cat("\u2705 ALL VISUALIZATIONS CREATED\n")
cat("========================================\n")
cat(sprintf("Output directory: %s\n", normalizePath(out_dir)))
cat("\n[NEW] CALENDAR-MONTHLY TEMPORAL CHANGE OUTPUTS:\n")
cat(" calmonth_temporal_Precipitation.pdf/.png — 12-panel: Jan–Dec Pr time series\n")
cat(" calmonth_temporal_PET_PM.pdf/.png — 12-panel: Jan–Dec PET (PM)\n")
cat(" calmonth_temporal_PET_Thw.pdf/.png — 12-panel: Jan–Dec PET (Thw)\n")
cat(" calmonth_temporal_Temperature.pdf/.png — 12-panel: Jan–Dec Tair\n")
cat(" calmonth_slope_profile_all_vars.pdf/.png — combined per-month slope profile\n")
cat(" calmonth_ols_slope_profiles.csv — slope + significance table\n")
cat("\n[NEW] PET COMPARISON OUTPUTS:\n")
cat(" pet_comparison_14a_timeseries.pdf/.png — PM vs Thw monthly time series\n")
cat(" pet_comparison_14b_monthly_slopes.pdf/.png — per-month slope bar chart\n")
cat(" pet_comparison_14c_trendsmap.pdf — side-by-side Sen's slope maps\n")
cat(" pet_comparison_14d_monthly_bias.pdf/.png — monthly bias (PM − Thw)\n")
cat(" pet_comparison_14e_scatter.pdf/.png — per-month PM vs Thw scatter\n")
cat("\n[NEW] SEASONALITY OUTPUT:\n")
cat(" FigS_Seasonality_P_PET.png — seasonal cycle P vs PET (two-panel bar chart)\n")
cat("\nSPATIAL MAPS: all 4 variables (Pr, PET_PM, PET_Thw, Temperature)\n")
cat("SPECIFIC-POINT: 4-panel TS + 4-panel slope bars at each point\n")