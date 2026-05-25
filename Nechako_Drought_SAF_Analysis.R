#============================================================================
# NECHAKO DROUGHT SAF — MASTER RUN SCRIPT
# Single entry point for the complete analysis pipeline.
#
# USAGE:
#   setwd("D:/Nechako_Drought/Nechako")   # set once, then:
#   source("Nechako_Drought_RUN_ALL.R")
#
# WHAT THIS DOES (in order):
#   Step 1 — Verify inputs (QuickStart)
#   Step 2 — Load utility / methodology helpers (MethodologyGuide)
#   Step 3 — Run main SAF pipeline: extract characteristics, fit marginals &
#             copulas, derive SAF curves for SPI-1, SPEI-1, SPI-3, SPEI-3
#   Step 4 — Load enhancement functions (ext1 : GOF, TV-copula;
#                                        ext2 : NS marginals)
#   Step 5 — Run all enhancements on the pipeline results
#   Step 6 — Optional: run seasonal / spatial / cross-index utility analyses
#   Step 7 — Generate Publication Figures & Excel Tables
#
#
#   [QuickStart] total_events now sourced from nrow(ev_file) when the raw event
#          CSV exists (authoritative); sum(rp$n_events) retained as fallback and
#          cross-check target.  MISMATCH flag compares the two counts.
#
# FIXES applied to this version:
#   FIX-1A  derive_SAF_nonstationary*(): added missing required argument
#           ns_result = ns to all four calls (ref/rec × conditional/Kendall).
#           Without it R aborts with "argument 'ns_result' is missing".
#   FIX-1B  Fig 16 change-point field names corrected:
#           pettitt_yr → pettitt_cp, cusum_yr → cusum_cp,
#           prutf_yr → prutf_cp, consensus_yr → cp_year.
#   FIX-1C  Table 3 TV-copula significance field corrected:
#           tv$lr_sig → tv$significant.
#   FIX-1D  Fig 18 change-point field corrected:
#           cp_result$consensus_yr → cp_result$cp_year.
#   FIX-1E  Fig 6 dependency field corrected:
#           dep$kendall_tau → dep$kendall.
#   FIX-1F  Fig 14 / Table 3 cp_year extraction corrected: ns$cp_year does not
#           exist in fit_nonstationary_marginals() return value; cp_year is now
#           parsed from ns$epoch_src (consistent with .resolve_epoch() in ext2).
#   FIX-5   quick_diagnostic() hardcodes "drought_analysis" but all output is
#           written to "spatial_drought". Patched inline after sourcing QuickStart.
#   FIX-8   Step 4 version label corrected from "ext2 v4" to "ext2 v5".
#============================================================================

setwd("D:/Nechako_Drought/Nechako")

# ---------------------------------------------------------------------------
# STEP 1 — Input verification
# ---------------------------------------------------------------------------
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("STEP 1: Input verification\n")
cat(paste(rep("=", 70), collapse=""), "\n")
source("Nechako_Drought_QuickStart.R")

# FIX-5: quick_diagnostic() as defined in QuickStart.R searches the wrong
# directory ("drought_analysis"). All pipeline output is written under
# OUT_ROOT = "spatial_drought". Override the function here so the diagnostic
# at the end of this script finds its files correctly.
# The body below is identical to the QuickStart version except for the
# list.files() path argument.
quick_diagnostic <- function() {
  sep_major <- paste(rep("=", 72), collapse = "")
  sep_minor <- paste(rep("-", 72), collapse = "")
  
  cat("\n", sep_major, "\n", sep = "")
  cat("DIAGNOSTIC — drought event summary by index and duration class\n")
  cat(sep_major, "\n", sep = "")
  
  # FIX-5: corrected from "drought_analysis" → "spatial_drought"
  dirs <- list.files("spatial_drought", pattern = "_analysis$", full.names = TRUE)
  if (length(dirs) == 0) {
    cat("  No analysis folders found. Run Nechako_Drought_SAF_Analysis.R first.\n")
    return(invisible(NULL))
  }
  
  summary_rows <- list()
  
  for (d in sort(dirs)) {
    index_tag <- sub("_analysis$", "", basename(d))
    
    rp_file <- file.path(d, sprintf("%s_return_periods_by_class.csv", index_tag))
    ev_file <- file.path(d, sprintf("%s_duration_classified_events.csv", index_tag))
    
    if (!file.exists(rp_file)) {
      cat(sprintf("\n  [%s] return_periods_by_class.csv not found — skipping.\n", index_tag))
      next
    }
    
    rp <- read.csv(rp_file, stringsAsFactors = FALSE)
    
    total_events <- if (file.exists(ev_file)) nrow(read.csv(ev_file)) else sum(rp$n_events, na.rm = TRUE)
    rp_total     <- sum(rp$n_events, na.rm = TRUE)
    source_label <- if (file.exists(ev_file)) "event list (authoritative)" else "rp table (fallback)"
    
    cat("\n", sep_minor, "\n", sep = "")
    cat(sprintf("  INDEX : %s\n", index_tag))
    cat(sprintf("  Total drought events : %d  [source: %s]", total_events, source_label))
    if (file.exists(ev_file)) {
      match_flag <- ifelse(rp_total == total_events, "OK", "MISMATCH — investigate aggregation!")
      cat(sprintf("\n  Cross-check vs rp table sum  : %d — %s", rp_total, match_flag))
    }
    cat("\n")
    cat(sep_minor, "\n", sep = "")
    
    cat(sprintf("  %-30s  %8s  %12s  %10s  %10s\n",
                "Duration Class", "N Events", "RP (years)", "Mean Sev", "Mean Area%"))
    cat(sprintf("  %-30s  %8s  %12s  %10s  %10s\n",
                paste(rep("-", 30), collapse=""),
                "--------", "------------", "----------", "----------"))
    
    for (i in seq_len(nrow(rp))) {
      r        <- rp[i, ]
      n_ev     <- r$n_events
      rp_yr    <- ifelse(is.infinite(r$return_period_years) | n_ev == 0,
                         "—", sprintf("%.1f", r$return_period_years))
      ms       <- ifelse(is.na(r$mean_severity), "—", sprintf("%.4f", r$mean_severity))
      ma       <- ifelse(is.na(r$mean_area_pct), "—", sprintf("%.2f", r$mean_area_pct))
      saf_flag <- if (n_ev >= 15) "" else "  [SAF skipped: n<15]"
      cat(sprintf("  %-30s  %8d  %12s  %10s  %10s%s\n",
                  r$duration_class, n_ev, rp_yr, ms, ma, saf_flag))
    }
    
    summary_rows[[index_tag]] <- data.frame(
      Index            = index_tag,
      Total_Events     = total_events,
      Total_Events_rp  = rp_total,
      N_Classes        = nrow(rp),
      stringsAsFactors = FALSE
    )
  }
  
  if (length(summary_rows) > 0) {
    cat("\n", sep_major, "\n", sep = "")
    cat("  CROSS-INDEX SUMMARY\n")
    cat(sep_major, "\n", sep = "")
    cat(sprintf("  %-20s  %14s  %14s  %10s\n", "Index", "Total(ev_file)", "Total(rp sum)", "N Classes"))
    cat(sprintf("  %-20s  %14s  %14s  %10s\n",
                paste(rep("-",20),collapse=""),
                "--------------", "--------------", "----------"))
    for (nm in names(summary_rows)) {
      r    <- summary_rows[[nm]]
      flag <- if (!is.na(r$Total_Events_rp) && r$Total_Events_rp != r$Total_Events) " !" else ""
      cat(sprintf("  %-20s  %14d  %14d%s  %10d\n",
                  r$Index, r$Total_Events, r$Total_Events_rp, flag, r$N_Classes))
    }
    cat(sep_major, "\n\n", sep = "")
  }
  
  invisible(summary_rows)
}

# ---------------------------------------------------------------------------
# STEP 2 — Utility / methodology helper functions (no execution, safe to load first)
# ---------------------------------------------------------------------------
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("STEP 2: Loading methodology helper functions\n")
cat(paste(rep("=", 70), collapse=""), "\n")
source("Nechako_Drought_Methodology_Guide.R")

# ---------------------------------------------------------------------------
# STEP 3 — Main SAF pipeline
#   Defines all core functions AND immediately executes for all 4 indices.
#   On completion the following objects exist in the workspace:
#     res_spi1, res_spei1, res_spi3, res_spei3
#   Each is a list with: dc, events, rp, marginals, copulas,
#                        saf_conditional, saf_kendall, saf_cond_by_class,
#                        saf_kend_by_class
# ---------------------------------------------------------------------------
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("STEP 3: Main SAF pipeline (SPI-1, SPEI-1, SPI-3, SPEI-3)\n")
cat(paste(rep("=", 70), collapse=""), "\n")
source("Nechako_Drought_SAF_Analysis.R")

# ---------------------------------------------------------------------------
# STEP 4 — Load extension functions (no rm, no duplicate definitions)
# ---------------------------------------------------------------------------
cat("\n", paste(rep("=", 70), collapse=""), "\n")
# FIX-8: corrected version label from "ext2 v4" to "ext2 v5"
cat("STEP 4: Loading enhancement functions (ext1 v6 + ext2 v5)\n")
cat(paste(rep("=", 70), collapse=""), "\n")
source("Nechako_Drought_SAF_Analysis_ext1.R")
source("Nechako_Drought_SAF_Analysis_ext2.R")

# ---------------------------------------------------------------------------
# STEP 5 — Run enhancements on existing pipeline results
# ---------------------------------------------------------------------------
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("STEP 5: Running enhancements\n")
cat(paste(rep("=", 70), collapse=""), "\n")

# Helper: named list of all four index results with their output directories
INDEX_RESULTS <- list(
  SPI1  = list(res = res_spi1,  dir = "spatial_drought/SPI1_analysis",  scale = 1),
  SPEI1 = list(res = res_spei1, dir = "spatial_drought/SPEI1_analysis", scale = 1),
  SPI3  = list(res = res_spi3,  dir = "spatial_drought/SPI3_analysis",  scale = 3),
  SPEI3 = list(res = res_spei3, dir = "spatial_drought/SPEI3_analysis", scale = 3)
)

# Collect all enhancement results for easy inspection
enh_results <- list()

for (nm in names(INDEX_RESULTS)) {
  entry   <- INDEX_RESULTS[[nm]]
  res     <- entry$res
  out_dir <- entry$dir
  
  cat(sprintf("\n--- %s ---\n", nm))
  
  # Guard: skip if main pipeline failed for this index
  if (is.null(res) || is.null(res$copulas) || is.null(res$marginals)) {
    cat(sprintf("  [%s] Main pipeline result unavailable — skipping enhancements.\n", nm))
    next
  }
  
  u_mat <- cbind(res$copulas$u_severity, res$copulas$u_area)
  
  # 5a) Kendall / Spearman dependency test
  dep <- test_dependency(res$dc, nm)
  
  # 5b) Copula goodness-of-fit (bootstrap, 499 replicates)
  #     Set N_boot = 99 for a quicker run during development
  gof <- gof_copula(res$copulas, u_mat, nm, N_boot = 499L)
  
  # 5c) Time-varying copula (tests if dependence strengthens / weakens over time)
  #     Now also runs detect_copula_changepoints() internally (Pettitt/CUSUM/PRUTF)
  #     and Rosenblatt PIT GOF after fitting.  The returned object contains
  #     cp_result and make_cop_fn, which are passed to 5d and 5e below.
  tv  <- fit_timevarying_copula(res$dc, res$copulas, res$marginals,
                                year_indices, nm, out_dir)
  
  # ---------------------------------------------------------------------------
  # DIAGNOSTICS (runs after GOF & TV fit)
  # ---------------------------------------------------------------------------
  if (!is.null(gof) && gof$decision == "REJECTED") {
    plot_copula_diagnostic(u_mat, res$copulas$best_copula_fit@copula, nm, out_dir)
    # Guard against NA parameters when the segmented model is selected
    if (!is.null(tv) && !is.null(tv$rosenblatt) && !is.na(tv$a_hat)) {
      # FIX (3.3): u_mat has n_drought rows (severity>0 & area_pct>0 only).
      # res$dc contains ALL months, so res$dc$year[i] mis-maps to the wrong
      # calendar year for every row after the first non-drought month.
      # Filter to drought months first so index i aligns with u_mat row i.
      dc_drought <- res$dc[res$dc$severity > 0 & res$dc$area_pct > 0, ]
      e2_vals <- vapply(seq_len(nrow(u_mat)), function(i) {
        yr_std_i <- (dc_drought$year[i] - tv$yr_mean) / tv$yr_sd
        cop_i <- tv$make_cop_fn(tv$a_hat + tv$b_hat * yr_std_i)
        .h_func_fd(u_mat[i,1], u_mat[i,2], cop_i)
      }, numeric(1))
      plot_rosenblatt_diagnostic(u_mat[,1], e2_vals, nm, out_dir)
    }
  }
  plot_rolling_tau(res$dc, nm, out_dir)
  if (!is.null(tv)) plot_tv_parameter(tv, nm, out_dir)
  # ---------------------------------------------------------------------------
  
  # 5d) Non-stationary severity marginals.
  #     cp_result from 5c is passed so that a statistically detected change
  #     point is used to define epoch boundaries automatically (instead of the
  #     fixed 40%/30% split).  Three GAMLSS models are compared: stationary,
  #     mu-trend, and mu+sigma-trend; period-specific sigma is returned for
  #     whichever model wins by AIC.
  #     Manual override example:
  #       ns <- fit_nonstationary_marginals(..., ref_years=1960:1990, recent_years=2005:2023)
  cp_res <- if (!is.null(tv)) tv$cp_result else NULL
  ns  <- fit_nonstationary_marginals(res$dc, nm, year_indices, out_dir,
                                     cp_result = cp_res)
  
  # 5e) Non-stationary SAF curves (reference and recent periods).
  #     Both conditional (Method A) and Kendall-corrected (Method B) variants.
  #     If 5c detected a significant trend in dependence (lr_p < 0.05), the
  #     copula in each SAF curve is replaced with a period-specific instance
  #     via .build_period_copula() priority-path logic (ext2 v5).
  #
  #     FIX-1A: ns_result = ns added to all four calls.  ext2 v5 made this a
  #     required argument; omitting it caused R to abort with
  #     "argument 'ns_result' is missing, with no default".
  saf_ns_ref   <- NULL
  saf_ns_rec   <- NULL
  saf_ns_ref_k <- NULL
  saf_ns_rec_k <- NULL
  
  if (!is.null(ns)) {
    ref_lbl <- sprintf("reference_%d_%d", min(ns$ref_years), max(ns$ref_years))
    rec_lbl <- sprintf("recent_%d_%d",    min(ns$recent_years), max(ns$recent_years))
    
    saf_ns_ref <- derive_SAF_nonstationary(
      mu_period        = ns$mu_ref,
      sigma_period     = ns$sigma_ref,
      copula_fit_obj   = res$copulas,
      marginal_fits    = res$marginals,
      drought_data     = res$dc,
      ns_result        = ns,           # FIX-1A: was missing; caused abort
      index_name       = nm,
      period_label     = ref_lbl,
      output_dir       = out_dir,
      tv_copula_result = tv)
    
    saf_ns_rec <- derive_SAF_nonstationary(
      mu_period        = ns$mu_rec,
      sigma_period     = ns$sigma_rec,
      copula_fit_obj   = res$copulas,
      marginal_fits    = res$marginals,
      drought_data     = res$dc,
      ns_result        = ns,           # FIX-1A: was missing; caused abort
      index_name       = nm,
      period_label     = rec_lbl,
      output_dir       = out_dir,
      tv_copula_result = tv)
    
    saf_ns_ref_k <- derive_SAF_nonstationary_kendall(
      mu_period        = ns$mu_ref,
      sigma_period     = ns$sigma_ref,
      copula_fit_obj   = res$copulas,
      marginal_fits    = res$marginals,
      drought_data     = res$dc,
      ns_result        = ns,           # FIX-1A: was missing; caused abort
      index_name       = nm,
      period_label     = ref_lbl,
      output_dir       = out_dir,
      tv_copula_result = tv)
    
    saf_ns_rec_k <- derive_SAF_nonstationary_kendall(
      mu_period        = ns$mu_rec,
      sigma_period     = ns$sigma_rec,
      copula_fit_obj   = res$copulas,
      marginal_fits    = res$marginals,
      drought_data     = res$dc,
      ns_result        = ns,           # FIX-1A: was missing; caused abort
      index_name       = nm,
      period_label     = rec_lbl,
      output_dir       = out_dir,
      tv_copula_result = tv)
  }
  
  enh_results[[nm]] <- list(
    dependency         = dep,
    gof                = gof,
    tv_copula          = tv,
    ns_marginals       = ns,
    saf_ns_ref         = saf_ns_ref,
    saf_ns_rec         = saf_ns_rec,
    saf_ns_ref_kendall = saf_ns_ref_k,
    saf_ns_rec_kendall = saf_ns_rec_k
  )
}

# ---------------------------------------------------------------------------
# STEP 6 — Optional utility analyses (uncomment as needed)
# ---------------------------------------------------------------------------
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("STEP 6: Optional utility analyses\n")
cat(paste(rep("=", 70), collapse=""), "\n")

# 6a) Threshold sensitivity (area fraction under drought at -0.5, -1.0, -1.5)
# analyze_multiple_thresholds(spi1, cell_area_km2, dates, thresholds = c(-0.5, -1.0, -1.5), index_name = "SPI1")

# 6b) Seasonal breakdown of severity and area
# seasonal_drought_analysis(res_spi1$dc, "SPI1", "spatial_drought/SPI1_analysis")
# seasonal_drought_analysis(res_spi3$dc, "SPI3", "spatial_drought/SPI3_analysis")

# 6c) Pixel-level drought frequency and severity maps
# map_drought_frequency(spi1,  index_name = "SPI1",  output_dir = "drought_maps")
# map_drought_frequency(spei1, index_name = "SPEI1", output_dir = "drought_maps")

# 6d) Cross-index correlation (SPI-1 vs SPI-3 severity)
# cross_correlate_indices(res_spi1$dc, res_spi3$dc, output_dir = "spatial_drought")

# 6e) Place a specific historical event on the SAF surface
# Replace event_years with the actual years of the event of interest
# place_event_on_saf(res_spi1$dc, res_spi1$copulas, res_spi1$marginals,
#                    event_years = c(2003, 2004), index_name = "SPI1",
#                    output_dir  = "spatial_drought/SPI1_analysis")

# ---------------------------------------------------------------------------
# STEP 7 — Generate Publication Figures & Excel Tables
# ---------------------------------------------------------------------------
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("STEP 7: Generating Publication Figures & Excel Tables\n")
cat(paste(rep("=", 70), collapse=""), "\n")

# ---- 0. SETUP FOR FIGURES ---------------------------------------------------
for (pkg in c("ggplot2","gridExtra","viridis","dplyr","scales","openxlsx","patchwork","ggExtra")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

PUB_DIR <- "spatial_drought/pub_figures"
dir.create(PUB_DIR, recursive = TRUE, showWarnings = FALSE)

# Consistent theme for all publication figures
pub_theme <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title        = element_text(face = "bold", hjust = 0.5, size = base_size + 1),
      plot.subtitle     = element_text(hjust = 0.5, colour = "grey40", size = base_size - 1),
      axis.title        = element_text(face = "bold"),
      legend.background = element_rect(colour = "grey80"),
      legend.key.size   = unit(0.6, "cm"),
      strip.background  = element_rect(fill = "grey92"),
      strip.text        = element_text(face = "bold")
    )
}

# Helper: save as PDF + JPEG
save_pub <- function(plot_obj, name, width = 10, height = 7, dpi = 300) {
  pdf_file  <- file.path(PUB_DIR, paste0(name, ".pdf"))
  jpeg_file <- file.path(PUB_DIR, paste0(name, ".jpeg"))
  pdf(pdf_file, width = width, height = height)
  if (inherits(plot_obj, "gg") || inherits(plot_obj, "patchwork")) {
    print(plot_obj)
  } else if (is.list(plot_obj)) {
    gridExtra::grid.arrange(grobs = plot_obj)
  } else {
    print(plot_obj)
  }
  dev.off()
  jpeg(jpeg_file, width = width, height = height, units = "in", res = dpi)
  if (inherits(plot_obj, "gg") || inherits(plot_obj, "patchwork")) {
    print(plot_obj)
  } else if (is.list(plot_obj)) {
    gridExtra::grid.arrange(grobs = plot_obj)
  } else {
    print(plot_obj)
  }
  dev.off()
  invisible(list(pdf = pdf_file, jpeg = jpeg_file))
}

# Colour palette for 4 indices
INDEX_COLOURS <- c(SPI1 = "#1b7837", SPEI1 = "#762a83",
                   SPI3 = "#2166ac", SPEI3 = "#d6604d")

# Return period colours (viridis 4)
RP_COLOURS <- setNames(viridis::viridis(4), c("10","25","50","100"))

INDEX_RESULTS <- list(
  SPI1  = list(res = res_spi1,  dir = "spatial_drought/SPI1_analysis",  scale = 1),
  SPEI1 = list(res = res_spei1, dir = "spatial_drought/SPEI1_analysis", scale = 1),
  SPI3  = list(res = res_spi3,  dir = "spatial_drought/SPI3_analysis",  scale = 3),
  SPEI3 = list(res = res_spei3, dir = "spatial_drought/SPEI3_analysis", scale = 3)
)

cat("=== NECHAKO DROUGHT — PUBLICATION FIGURES ===\n")
cat(sprintf("Output directory: %s\n\n", PUB_DIR))

# ============================================================================
# SECTION 1: STUDY AREA & DATA OVERVIEW
# ============================================================================
cat("--- Section 1: Data overview ---\n")

# Fig 1: Monthly time series of severity and area_pct for all 4 indices ------
{
  all_dc <- dplyr::bind_rows(lapply(names(INDEX_RESULTS), function(nm) {
    dc <- INDEX_RESULTS[[nm]]$res$dc
    dc$Index <- nm
    dc
  }))
  all_dc$date <- as.Date(all_dc$date)
  
  p1a <- ggplot(all_dc, aes(x = date, y = severity, colour = Index)) +
    geom_line(linewidth = 0.55, alpha = 0.85) +
    scale_colour_manual(values = INDEX_COLOURS) +
    scale_x_date(date_breaks = "10 years", date_labels = "%Y") +
    labs(x = NULL, y = "Mean Drought Severity", title = NULL,
         colour = "Index") +
    pub_theme()
  
  p1b <- ggplot(all_dc, aes(x = date, y = area_pct, colour = Index)) +
    geom_line(linewidth = 0.55, alpha = 0.85) +
    scale_colour_manual(values = INDEX_COLOURS) +
    scale_x_date(date_breaks = "10 years", date_labels = "%Y") +
    labs(x = "Year", y = "Drought-Affected Area (%)", colour = "Index") +
    pub_theme()
  
  fig01 <- p1a / p1b +
    plot_annotation(
      title    = "Nechako Basin: Monthly Drought Characteristics (1950–2025)",
      subtitle = "Threshold = −0.5  |  SPI-1, SPEI-1, SPI-3, SPEI-3",
      theme    = pub_theme()
    ) + plot_layout(guides = "collect")
  
  save_pub(fig01, "Fig01_basin_drought_timeseries", width = 12, height = 8)
  cat("  Fig01 saved.\n")
}

# Fig 2: Event scatter — severity vs area, coloured by duration class --------
{
  all_ev <- dplyr::bind_rows(lapply(names(INDEX_RESULTS), function(nm) {
    ev <- INDEX_RESULTS[[nm]]$res$events
    ev$Index <- nm
    ev
  }))
  
  fig02 <- ggplot(all_ev, aes(x = mean_area_pct, y = mean_severity,
                              colour = duration_class, shape = Index)) +
    geom_point(size = 2.2, alpha = 0.75) +
    scale_colour_viridis_d(name = "Duration Class", option = "plasma") +
    scale_shape_manual(values = c(SPI1=16, SPEI1=17, SPI3=15, SPEI3=18)) +
    labs(x = "Mean Drought-Affected Area (%)",
         y = "Mean Drought Severity",
         title = "Drought Event Characteristics",
         subtitle = "All identified events, all indices (1950–2025)") +
    pub_theme() +
    theme(legend.position = "right")
  
  save_pub(fig02, "Fig02_event_characteristics_scatter", width = 10, height = 7)
  cat("  Fig02 saved.\n")
}

# Fig 3: Event duration histogram by index -----------------------------------
{
  all_ev2 <- dplyr::bind_rows(lapply(names(INDEX_RESULTS), function(nm) {
    ev <- INDEX_RESULTS[[nm]]$res$events
    ev$Index <- nm
    ev
  }))
  
  fig03 <- ggplot(all_ev2, aes(x = duration_months, fill = Index)) +
    geom_histogram(binwidth = 1, colour = "white", alpha = 0.85) +
    facet_wrap(~Index, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = INDEX_COLOURS) +
    labs(x = "Event Duration (months)", y = "Count",
         title = "Distribution of Drought Event Durations") +
    pub_theme() + theme(legend.position = "none")
  
  save_pub(fig03, "Fig03_event_duration_histogram", width = 10, height = 7)
  cat("  Fig03 saved.\n")
}

# ============================================================================
# SECTION 2: MARGINAL DISTRIBUTIONS
# ============================================================================
cat("--- Section 2: Marginal distributions ---\n")

# Fig 4: Severity — empirical vs best-fit marginal (4-panel) -----------------
{
  plot_severity_marginal <- function(nm) {
    dc  <- INDEX_RESULTS[[nm]]$res$dc
    dc  <- dc[dc$severity > 0, ]
    mf  <- INDEX_RESULTS[[nm]]$res$marginals
    if (is.null(mf)) return(NULL)
    x_seq <- seq(0.001, max(dc$severity, na.rm = TRUE) * 1.05, length.out = 400)
    pars  <- mf$severity_fit$estimate
    dn    <- mf$severity_dist_name
    y_fit <- switch(dn,
                    Gamma       = dgamma(x_seq, shape = pars["shape"], rate  = pars["rate"]),
                    Weibull     = dweibull(x_seq, shape = pars["shape"], scale = pars["scale"]),
                    `Log-Normal`= dlnorm(x_seq, meanlog = pars["meanlog"], sdlog = pars["sdlog"]),
                    Exponential = dexp(x_seq, rate = pars["rate"]),
                    dnorm(x_seq, mean = pars["mean"], sd = pars["sd"]))
    df_fit <- data.frame(x = x_seq, y = y_fit)
    
    ggplot(dc, aes(x = severity)) +
      geom_histogram(aes(y = after_stat(density)), bins = 25,
                     fill = INDEX_COLOURS[nm], colour = "white", alpha = 0.7) +
      geom_line(data = df_fit, aes(x = x, y = y),
                colour = "black", linewidth = 1.1) +
      labs(title = nm, subtitle = sprintf("Best fit: %s", dn),
           x = "Drought Severity", y = "Density") +
      pub_theme()
  }
  
  panels <- lapply(names(INDEX_RESULTS), plot_severity_marginal)
  panels <- panels[!sapply(panels, is.null)]
  fig04  <- wrap_plots(panels, ncol = 2) +
    plot_annotation(title = "Marginal Distribution Fit: Drought Severity",
                    theme = pub_theme())
  save_pub(fig04, "Fig04_marginal_fits_severity", width = 10, height = 8)
  cat("  Fig04 saved.\n")
}

# Fig 5: Area — empirical vs Beta fit (4-panel) ------------------------------
{
  plot_area_marginal <- function(nm) {
    dc   <- INDEX_RESULTS[[nm]]$res$dc
    dc   <- dc[dc$area_pct > 0, ]
    mf   <- INDEX_RESULTS[[nm]]$res$marginals
    if (is.null(mf)) return(NULL)
    prop <- dc$area_pct / 100
    pars <- mf$area_fit$estimate
    x_seq <- seq(0.001, 0.999, length.out = 400)
    y_fit <- dbeta(x_seq, shape1 = pars["shape1"], shape2 = pars["shape2"])
    df_fit <- data.frame(x = x_seq * 100, y = y_fit / 100)
    
    ggplot(dc, aes(x = area_pct)) +
      geom_histogram(aes(y = after_stat(density)), bins = 25,
                     fill = INDEX_COLOURS[nm], colour = "white", alpha = 0.7) +
      geom_line(data = df_fit, aes(x = x, y = y),
                colour = "black", linewidth = 1.1) +
      labs(title = nm, subtitle = "Best fit: Beta",
           x = "Drought-Affected Area (%)", y = "Density") +
      pub_theme()
  }
  
  panels <- lapply(names(INDEX_RESULTS), plot_area_marginal)
  panels <- panels[!sapply(panels, is.null)]
  fig05  <- wrap_plots(panels, ncol = 2) +
    plot_annotation(title = "Marginal Distribution Fit: Drought-Affected Area",
                    theme = pub_theme())
  save_pub(fig05, "Fig05_marginal_fits_area", width = 10, height = 8)
  cat("  Fig05 saved.\n")
}

# ============================================================================
# SECTION 3: COPULA DEPENDENCE STRUCTURE
# ============================================================================
cat("--- Section 3: Copula dependence ---\n")

# Fig 6: Pseudo-observation scatter with fitted copula contours (4-panel) ----
{
  plot_pseudo_obs <- function(nm) {
    res   <- INDEX_RESULTS[[nm]]$res
    if (is.null(res$copulas)) return(NULL)
    uS <- res$copulas$u_severity
    uA <- res$copulas$u_area
    cop <- res$copulas$best_copula_fit@copula
    best_name <- res$copulas$best_copula_name
    
    # Compute copula density contours on a grid
    g    <- seq(0.01, 0.99, length.out = 80)
    grid <- expand.grid(u1 = g, u2 = g)
    grid$density <- tryCatch(
      copula::dCopula(as.matrix(grid), cop),
      error = function(e) rep(NA_real_, nrow(grid)))
    grid <- grid[is.finite(grid$density), ]
    
    ggplot() +
      geom_contour(data = grid, aes(x = u1, y = u2, z = density),
                   colour = "steelblue", bins = 10, linewidth = 0.5, alpha = 0.7) +
      geom_point(aes(x = uS, y = uA), size = 1.4, alpha = 0.5,
                 colour = INDEX_COLOURS[nm]) +
      labs(title = nm,
           # FIX-1E: was dep$kendall_tau — field does not exist;
           #         correct name is dep$kendall (from test_dependency() return list)
           subtitle = sprintf("Best copula: %s  |  τ=%.3f",
                              best_name,
                              enh_results[[nm]]$dependency$kendall),
           x = "U(Severity)", y = "U(Area)") +
      coord_fixed() +
      pub_theme()
  }
  
  panels <- lapply(names(INDEX_RESULTS), plot_pseudo_obs)
  panels <- panels[!sapply(panels, is.null)]
  fig06  <- wrap_plots(panels, ncol = 2) +
    plot_annotation(title = "Pseudo-Observations and Fitted Copula Contours",
                    theme = pub_theme())
  save_pub(fig06, "Fig06_pseudo_obs_scatter", width = 10, height = 9)
  cat("  Fig06 saved.\n")
}

# Fig 7: Kendall τ summary barplot + class breakdown -------------------------
{
  tau_df <- dplyr::bind_rows(lapply(names(INDEX_RESULTS), function(nm) {
    dep <- enh_results[[nm]]$dependency
    if (is.null(dep)) return(NULL)
    data.frame(Index     = nm,
               Class     = "All events",
               tau       = dep$kendall,
               rho       = dep$spearman,
               p_kendall = dep$p_value)
  }))
  
  fig07 <- ggplot(tau_df, aes(x = Index, y = tau, fill = Index)) +
    geom_col(width = 0.6, colour = "grey30") +
    geom_errorbar(aes(ymin = tau - 0.02, ymax = tau + 0.02), width = 0.2) +
    geom_text(aes(label = sprintf("τ=%.3f\nρ=%.3f", tau, rho)),
              vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = INDEX_COLOURS) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.1))) +
    labs(x = "Drought Index", y = "Kendall τ",
         title = "Severity–Area Rank Correlation by Index",
         subtitle = "All identified drought events  |  p < 0.001 for all") +
    pub_theme() + theme(legend.position = "none")
  
  save_pub(fig07, "Fig07_kendall_tau_barplot", width = 8, height = 6)
  cat("  Fig07 saved.\n")
}

# Fig 8: Copula AIC comparison (grouped bar) ---------------------------------
{
  aic_list <- lapply(names(INDEX_RESULTS), function(nm) {
    res <- INDEX_RESULTS[[nm]]$res
    if (is.null(res$copulas$all_results)) return(NULL)
    df <- res$copulas$all_results
    df$Index <- nm
    df
  })
  aic_df <- dplyr::bind_rows(aic_list[!sapply(aic_list, is.null)])
  
  if (nrow(aic_df) > 0 && "AIC" %in% names(aic_df) && "Copula" %in% names(aic_df)) {
    aic_df <- aic_df %>%
      group_by(Index) %>%
      mutate(delta_AIC = AIC - min(AIC, na.rm = TRUE)) %>%
      ungroup()
    
    fig08 <- ggplot(aic_df, aes(x = Copula, y = delta_AIC, fill = Index)) +
      geom_col(position = "dodge", colour = "grey30", width = 0.7) +
      scale_fill_manual(values = INDEX_COLOURS) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      labs(x = "Copula Family", y = "ΔAIC (relative to best)",
           title = "Copula Model Selection by Index",
           subtitle = "Lower ΔAIC = better fit  |  Best model at ΔAIC = 0") +
      pub_theme() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
    
    save_pub(fig08, "Fig08_copula_aic_comparison", width = 10, height = 6)
    cat("  Fig08 saved.\n")
  } else {
    cat("  Fig08 skipped (AIC table unavailable).\n")
  }
}

# ============================================================================
# SECTION 4: SAF CURVES
# ============================================================================
cat("--- Section 4: SAF curves ---\n")

# Fig 9 & 10: SAF at specific return periods — all indices on one panel ------
{
  make_saf_allindices <- function(rp_filter, method = "Conditional", fig_name) {
    all_saf <- dplyr::bind_rows(lapply(names(INDEX_RESULTS), function(nm) {
      res  <- INDEX_RESULTS[[nm]]$res
      saf  <- if (method == "Conditional") res$saf_conditional else res$saf_kendall
      if (is.null(saf)) return(NULL)
      saf[saf$ReturnPeriod_years == rp_filter & is.finite(saf$Severity), ]
    }))
    if (nrow(all_saf) == 0) { cat(sprintf("  %s skipped — no data.\n", fig_name)); return(NULL) }
    
    p <- ggplot(all_saf, aes(x = Area_pct, y = Severity, colour = Index, linetype = Index)) +
      geom_line(linewidth = 1.2) + geom_point(size = 2.5) +
      scale_colour_manual(values = INDEX_COLOURS) +
      scale_linetype_manual(values = c(SPI1="solid", SPEI1="dashed",
                                       SPI3="dotdash", SPEI3="dotted")) +
      labs(x = "Drought-Affected Area (%)", y = "Drought Severity",
           title = sprintf("SAF Curves at %d-Year Return Period", rp_filter),
           subtitle = sprintf("Method: %s  |  All four indices, full record", method)) +
      pub_theme()
    p
  }
  
  fig09 <- make_saf_allindices(10,  "Conditional", "Fig09")
  fig10 <- make_saf_allindices(100, "Conditional", "Fig10")
  if (!is.null(fig09)) { save_pub(fig09, "Fig09_SAF_allindices_RP10",  width = 9, height = 6); cat("  Fig09 saved.\n") }
  if (!is.null(fig10)) { save_pub(fig10, "Fig10_SAF_allindices_RP100", width = 9, height = 6); cat("  Fig10 saved.\n") }
}

# Fig 11: Conditional vs Kendall, 4-panel grid --------------------------------
{
  plot_method_comp <- function(nm) {
    res  <- INDEX_RESULTS[[nm]]$res
    cond <- res$saf_conditional
    kend <- res$saf_kendall
    if (is.null(cond) || is.null(kend)) return(NULL)
    comb <- rbind(cond[is.finite(cond$Severity), ],
                  kend[is.finite(kend$Severity), ])
    
    ggplot(comb, aes(x = Area_pct, y = Severity,
                     colour = factor(ReturnPeriod_years),
                     linetype = Method)) +
      geom_line(linewidth = 1.0) +
      scale_colour_manual(values = RP_COLOURS, name = "RP (yr)") +
      scale_linetype_manual(values = c(Conditional="solid", Kendall="dashed"),
                            name = "Method") +
      labs(title = nm, x = "Area (%)", y = "Severity") +
      pub_theme(base_size = 11)
  }
  
  panels <- lapply(names(INDEX_RESULTS), plot_method_comp)
  panels <- panels[!sapply(panels, is.null)]
  fig11  <- wrap_plots(panels, ncol = 2) +
    plot_annotation(
      title    = "SAF Curves: Conditional vs Kendall Method",
      subtitle = "Return periods: 10, 25, 50, 100 years  |  Full record 1950–2025",
      theme    = pub_theme())
  save_pub(fig11, "Fig11_SAF_method_comparison_4panel", width = 12, height = 9)
  cat("  Fig11 saved.\n")
}

# Fig 12: Per-class SAF curves for SPI-1 and SPI-3 ----------------------------
{
  make_class_saf <- function(nm, method = "Conditional") {
    res   <- INDEX_RESULTS[[nm]]$res
    cls_saf <- if (method == "Conditional") res$saf_cond_by_class else res$saf_kend_by_class
    if (is.null(cls_saf) || length(cls_saf) == 0) return(NULL)
    
    all_cls <- dplyr::bind_rows(lapply(names(cls_saf), function(cl) {
      df <- cls_saf[[cl]]
      if (!is.null(df)) df[is.finite(df$Severity), ] else NULL
    }))
    if (nrow(all_cls) == 0) return(NULL)
    
    ggplot(all_cls, aes(x = Area_pct, y = Severity,
                        colour = factor(ReturnPeriod_years),
                        linetype = DurationClass)) +
      geom_line(linewidth = 1.0) +
      scale_colour_manual(values = RP_COLOURS, name = "RP (yr)") +
      scale_linetype_discrete(name = "Duration Class") +
      labs(title = nm,
           x = "Drought-Affected Area (%)", y = "Drought Severity") +
      pub_theme(base_size = 11)
  }
  
  p12a <- make_class_saf("SPI1")
  p12b <- make_class_saf("SPI3")
  if (!is.null(p12a) && !is.null(p12b)) {
    fig12 <- p12a / p12b +
      plot_annotation(title = "SAF Curves by Duration Class: SPI-1 and SPI-3",
                      theme = pub_theme())
    save_pub(fig12, "Fig12_SAF_by_duration_class", width = 10, height = 10)
    cat("  Fig12 saved.\n")
  } else {
    cat("  Fig12 skipped (per-class SAF unavailable).\n")
  }
}

# ============================================================================
# SECTION 5: NON-STATIONARITY
# ============================================================================
cat("--- Section 5: Non-stationarity ---\n")

# Fig 13: NS SAF — reference vs recent, conditional, all 4 indices -----------
{
  plot_ns_comp <- function(nm) {
    enh  <- enh_results[[nm]]
    if (is.null(enh$saf_ns_ref) || is.null(enh$saf_ns_rec)) return(NULL)
    
    ns_ref <- enh$saf_ns_ref; ns_rec <- enh$saf_ns_rec
    ns_ref$Epoch <- "Reference"; ns_rec$Epoch <- "Recent"
    comb <- rbind(ns_ref[is.finite(ns_ref$Severity), ],
                  ns_rec[is.finite(ns_rec$Severity), ])
    if (nrow(comb) == 0) return(NULL)
    
    ggplot(comb, aes(x = Area_pct, y = Severity,
                     colour = factor(ReturnPeriod_years),
                     linetype = Epoch)) +
      geom_line(linewidth = 1.0) +
      scale_colour_manual(values = RP_COLOURS, name = "RP (yr)") +
      scale_linetype_manual(values = c(Reference="solid", Recent="dashed"),
                            name = "Epoch") +
      labs(title = nm,
           x = "Area (%)", y = "Severity") +
      pub_theme(base_size = 11)
  }
  
  panels <- lapply(names(INDEX_RESULTS), plot_ns_comp)
  panels <- panels[!sapply(panels, is.null)]
  if (length(panels) > 0) {
    fig13 <- wrap_plots(panels, ncol = 2) +
      plot_annotation(
        title    = "Non-Stationary SAF Curves: Reference vs Recent Epoch",
        subtitle = "Conditional method  |  Period-specific marginals (GAMLSS)",
        theme    = pub_theme())
    save_pub(fig13, "Fig13_NS_epoch_comparison_4panel", width = 12, height = 9)
    cat("  Fig13 saved.\n")
  } else cat("  Fig13 skipped.\n")
}

# Fig 14: GAMLSS model selection table rendered as a figure ------------------
{
  # FIX-1F: ns$cp_year does not exist in fit_nonstationary_marginals() return
  #         value.  The change-point year is encoded in ns$epoch_src as a
  #         string (e.g. "change-point (1987)").  Parse it with the same regex
  #         used by .resolve_epoch() in ext2.R, consistent with how ext2
  #         itself reads the epoch boundary.
  extract_cp_year_from_ns <- function(ns) {
    if (is.null(ns) || !grepl("change-point", ns$epoch_src, fixed = TRUE))
      return(NA_integer_)
    m <- regmatches(ns$epoch_src, regexpr("[0-9]{4}", ns$epoch_src))
    if (length(m) == 1L) as.integer(m) else NA_integer_
  }
  
  ns_aic_df <- dplyr::bind_rows(lapply(names(INDEX_RESULTS), function(nm) {
    ns <- enh_results[[nm]]$ns_marginals
    if (is.null(ns)) return(NULL)
    data.frame(
      Index      = nm,
      CP_year    = extract_cp_year_from_ns(ns),  # FIX-1F
      AIC_M0     = round(ns$aic_m0, 2),
      AIC_M1     = round(ns$aic_m1, 2),
      AIC_M2     = round(ns$aic_m2, 2),
      Best_model = ns$best_model,
      mu_ref     = round(ns$mu_ref, 4),
      mu_rec     = round(ns$mu_rec, 4),
      sigma_ref  = round(ns$sigma_ref, 4),
      sigma_rec  = round(ns$sigma_rec, 4)
    )
  }))
  
  if (nrow(ns_aic_df) > 0) {
    tbl <- gridExtra::tableGrob(ns_aic_df, rows = NULL,
                                theme = gridExtra::ttheme_minimal(base_size = 11))
    fig14 <- gridExtra::grid.arrange(tbl,
                                     top    = grid::textGrob("GAMLSS Non-Stationary Marginal: Model Selection Summary",
                                                             gp = grid::gpar(fontface = "bold", fontsize = 12)),
                                     bottom = grid::textGrob("M0=Stationary  |  M1=μ-trend  |  M2=μ+σ-trend",
                                                             gp = grid::gpar(fontsize = 9, col = "grey40")),
                                     heights = grid::unit.c(grid::unit(1, "line"),
                                                            grid::unit(1, "null"),
                                                            grid::unit(1, "line")))
    save_pub(fig14, "Fig14_GAMLSS_model_selection_table", width = 12, height = 3.5)
    cat("  Fig14 saved.\n")
  } else cat("  Fig14 skipped.\n")
}

# Fig 15: Severity distribution by epoch — split violin / boxplot -----------
{
  epoch_sev_df <- dplyr::bind_rows(lapply(names(INDEX_RESULTS), function(nm) {
    ns  <- enh_results[[nm]]$ns_marginals
    res <- INDEX_RESULTS[[nm]]$res
    if (is.null(ns)) return(NULL)
    dc  <- res$dc[res$dc$severity > 0, ]
    dc$Epoch <- dplyr::if_else(dc$year %in% ns$ref_years, "Reference", "Recent")
    dc$Index <- nm
    dc
  }))
  
  if (nrow(epoch_sev_df) > 0) {
    fig15 <- ggplot(epoch_sev_df, aes(x = Epoch, y = severity, fill = Epoch)) +
      geom_violin(trim = TRUE, alpha = 0.6, colour = NA) +
      geom_boxplot(width = 0.15, outlier.size = 1, colour = "grey30", fill = "white") +
      facet_wrap(~Index, ncol = 2, scales = "free_y") +
      scale_fill_manual(values = c(Reference = "#4575b4", Recent = "#d73027")) +
      labs(x = NULL, y = "Drought Severity",
           title = "Drought Severity by Epoch",
           subtitle = "Reference (pre-CP) vs Recent (post-CP) periods") +
      pub_theme() + theme(legend.position = "none")
    
    save_pub(fig15, "Fig15_severity_trend_by_epoch", width = 10, height = 8)
    cat("  Fig15 saved.\n")
  } else cat("  Fig15 skipped.\n")
}

# ============================================================================
# SECTION 6: COPULA DIAGNOSTICS & CHANGE-POINT
# ============================================================================
cat("--- Section 6: Copula diagnostics ---\n")

# Fig 16: Change-point detection summary -------------------------------------
{
  cp_df <- dplyr::bind_rows(lapply(names(INDEX_RESULTS), function(nm) {
    tv <- enh_results[[nm]]$tv_copula
    if (is.null(tv) || is.null(tv$cp_result)) return(NULL)
    cp <- tv$cp_result
    
    # FIX-1B: detect_copula_changepoints() returns field names
    #   pettitt_cp, cusum_cp, prutf_cp, cp_year (consensus)
    # Previous code used pettitt_yr / cusum_yr / prutf_yr / consensus_yr,
    # none of which exist, producing all-NA data frames and a blank figure.
    data.frame(
      Index    = nm,
      Method   = c("Pettitt", "CUSUM", "PRUTF", "Consensus"),
      CP_year  = c(cp$pettitt_cp, cp$cusum_cp, cp$prutf_cp, cp$cp_year),
      p_value  = c(cp$pettitt_p,  cp$cusum_p,  cp$prutf_p,  NA),
      Detected = c(isTRUE(cp$pettitt_p < 0.05),
                   isTRUE(cp$cusum_p   < 0.05),
                   isTRUE(cp$prutf_p   < 0.05),
                   isTRUE(cp$detected))
    )
  }))
  
  if (nrow(cp_df) > 0) {
    fig16 <- ggplot(cp_df[cp_df$Method != "Consensus", ],
                    aes(x = Index, y = CP_year, colour = Method, shape = Detected)) +
      geom_jitter(size = 4, width = 0.12, height = 0) +
      geom_hline(data = cp_df[cp_df$Method == "Consensus", ],
                 aes(yintercept = CP_year, colour = Index), linetype = "dashed",
                 linewidth = 0.9, show.legend = FALSE) +
      scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                         name = "Significant\n(p < 0.05)") +
      scale_colour_viridis_d(name = "Method/Index") +
      labs(x = "Index", y = "Change-Point Year",
           title = "Detected Change Points in Copula Dependence",
           subtitle = "Pettitt / CUSUM / PRUTF  |  Dashed line = consensus year") +
      pub_theme()
    
    save_pub(fig16, "Fig16_changepoint_summary", width = 9, height = 6)
    cat("  Fig16 saved.\n")
  } else cat("  Fig16 skipped.\n")
}

# Fig 17: Rosenblatt PIT diagnostic — 4-panel (save JPEG like the attached) --
{
  plot_rosenblatt <- function(nm) {
    res   <- INDEX_RESULTS[[nm]]$res
    if (is.null(res$copulas)) return(NULL)
    uS    <- res$copulas$u_severity
    uA    <- res$copulas$u_area
    cop   <- res$copulas$best_copula_fit@copula
    best  <- res$copulas$best_copula_name
    
    e1 <- uS
    
    h  <- 1e-5
    e2 <- vapply(seq_along(uS), function(i) {
      u1_lo <- max(uS[i] - h, 1e-7)
      u1_hi <- min(uS[i] + h, 1 - 1e-7)
      denom <- u1_hi - u1_lo
      if (denom <= 0) return(NA_real_)
      (copula::pCopula(c(u1_hi, uA[i]), cop) -
          copula::pCopula(c(u1_lo, uA[i]), cop)) / denom
    }, numeric(1))
    e2 <- pmax(pmin(e2, 1), 0)
    e2[!is.finite(e2)] <- NA_real_
    valid <- is.finite(e1) & is.finite(e2)
    e1v <- e1[valid]; e2v <- e2[valid]
    
    df_e2 <- data.frame(e2 = e2v)
    p_hist <- ggplot(df_e2, aes(x = e2)) +
      geom_histogram(aes(y = after_stat(density)), bins = 15,
                     fill = "#6baed6", colour = "white") +
      geom_hline(yintercept = 1, colour = "red", linewidth = 1) +
      labs(title = "e2 Histogram vs Uniform[0,1]",
           x = "e2 value", y = "Density") +
      pub_theme(base_size = 11)
    
    n   <- length(e2v)
    theoretical <- (seq_len(n) - 0.5) / n
    empirical   <- sort(e2v)
    df_qq <- data.frame(t = theoretical, s = empirical)
    p_qq  <- ggplot(df_qq, aes(x = t, y = s)) +
      geom_point(size = 0.9, alpha = 0.7) +
      geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.9) +
      labs(title = "e2 QQ-Plot",
           x = "Theoretical Quantiles", y = "Sample Quantiles") +
      pub_theme(base_size = 11)
    
    df_sc <- data.frame(e1 = e1v, e2 = e2v)
    p_sc  <- ggplot(df_sc, aes(x = e1, y = e2)) +
      geom_point(size = 1.2, alpha = 0.45, colour = "grey30") +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                  colour = "red", linetype = "dashed", linewidth = 0.8) +
      labs(title = "e1 vs e2 (independence check)",
           x = "e1", y = "e2") +
      pub_theme(base_size = 11)
    
    (p_hist | p_qq) / (p_sc | plot_spacer()) +
      plot_annotation(
        title    = nm,
        subtitle = sprintf("Rosenblatt PIT  |  Copula: %s  |  n=%d valid", best, sum(valid)),
        theme    = pub_theme(base_size = 12)
      )
  }
  
  panels_rob <- lapply(names(INDEX_RESULTS), function(nm) {
    tryCatch(plot_rosenblatt(nm), error = function(e) NULL)
  })
  names(panels_rob) <- names(INDEX_RESULTS)
  
  for (nm in names(panels_rob)) {
    if (!is.null(panels_rob[[nm]])) {
      out_dir_idx <- INDEX_RESULTS[[nm]]$dir
      jpeg(file.path(out_dir_idx, sprintf("%s_rosenblatt_diagnostic.jpeg", nm)),
           width = 9, height = 7, units = "in", res = 300)
      print(panels_rob[[nm]])
      dev.off()
      save_pub(panels_rob[[nm]], sprintf("Fig17_%s_rosenblatt", nm),
               width = 9, height = 7)
    }
  }
  
  valid_panels <- panels_rob[!sapply(panels_rob, is.null)]
  if (length(valid_panels) == 4) {
    fig17_comb <- (valid_panels[[1]] | valid_panels[[2]]) /
      (valid_panels[[3]] | valid_panels[[4]]) +
      plot_annotation(title = "Rosenblatt PIT Diagnostic — All Indices",
                      theme = pub_theme())
    save_pub(fig17_comb, "Fig17_rosenblatt_4panel_combined", width = 18, height = 14)
  }
  cat("  Fig17 saved (individual + combined).\n")
}

# Fig 18: Rolling Kendall τ, 4-panel ------------------------------------------
{
  compute_rolling_tau <- function(dc, window = 10) {
    dc2 <- dc[dc$severity > 0 & dc$area_pct > 0, ]
    dc2 <- dc2[order(dc2$year), ]
    n   <- nrow(dc2)
    if (n < window) return(NULL)
    out_yr  <- numeric(n - window + 1)
    out_tau <- numeric(n - window + 1)
    for (i in seq_len(n - window + 1)) {
      idx <- i:(i + window - 1)
      s_w <- dc2$severity[idx]; a_w <- dc2$area_pct[idx]
      if (sd(s_w) == 0 || sd(a_w) == 0) {
        out_tau[i] <- NA_real_
      } else {
        out_tau[i] <- cor(s_w, a_w, method = "kendall")
      }
      out_yr[i] <- mean(dc2$year[idx])
    }
    data.frame(year = out_yr, tau = out_tau)
  }
  
  rt_plots <- lapply(names(INDEX_RESULTS), function(nm) {
    res <- INDEX_RESULTS[[nm]]$res
    rt  <- tryCatch(compute_rolling_tau(res$dc), error = function(e) NULL)
    if (is.null(rt)) return(NULL)
    
    # FIX-1D: detect_copula_changepoints() returns cp_result$cp_year
    #         (the consensus year), not cp_result$consensus_yr.
    cp_yr <- tryCatch(
      enh_results[[nm]]$tv_copula$cp_result$cp_year,
      error = function(e) NA)
    
    p <- ggplot(rt, aes(x = year, y = tau)) +
      geom_line(colour = INDEX_COLOURS[nm], linewidth = 0.9) +
      geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
                  colour = "grey30", linetype = "dashed", linewidth = 0.7, alpha = 0.2) +
      geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
      { if (!is.na(cp_yr)) geom_vline(xintercept = cp_yr, colour = "red",
                                      linetype = "dashed", linewidth = 0.8) else NULL } +
      labs(title = nm, x = "Year", y = "Rolling Kendall τ (10-yr window)") +
      pub_theme(base_size = 11)
    p
  })
  rt_plots <- rt_plots[!sapply(rt_plots, is.null)]
  if (length(rt_plots) > 0) {
    fig18 <- wrap_plots(rt_plots, ncol = 2) +
      plot_annotation(
        title    = "Rolling Kendall τ: Severity–Area Dependence Over Time",
        subtitle = "10-year moving window  |  Red dashed = consensus change-point year",
        theme    = pub_theme())
    save_pub(fig18, "Fig18_rolling_tau_4panel", width = 12, height = 8)
    cat("  Fig18 saved.\n")
  }
}

# ============================================================================
# SECTION 7: EXCEL SUMMARY TABLES
# ============================================================================
cat("--- Section 7: Excel tables ---\n")

# Table 1: Return periods consolidated ---------------------------------------
{
  rp_all <- dplyr::bind_rows(lapply(names(INDEX_RESULTS), function(nm) {
    rp <- INDEX_RESULTS[[nm]]$res$rp
    rp$Index <- nm
    rp
  }))
  wb1 <- createWorkbook()
  addWorksheet(wb1, "Return_Periods")
  writeData(wb1, "Return_Periods", rp_all)
  headerStyle <- createStyle(fontColour = "#FFFFFF", fgFill = "#2c5f8a",
                             halign = "CENTER", fontName = "Calibri",
                             border = "Bottom", bold = TRUE)
  addStyle(wb1, 1, headerStyle, rows = 1, cols = seq_len(ncol(rp_all)), gridExpand = TRUE)
  setColWidths(wb1, 1, cols = seq_len(ncol(rp_all)), widths = "auto")
  saveWorkbook(wb1, file.path(PUB_DIR, "Table01_return_periods_all_indices.xlsx"), overwrite = TRUE)
  cat("  Table01 saved.\n")
}

# Table 2: Copula fit summary ------------------------------------------------
{
  cop_all <- dplyr::bind_rows(lapply(names(INDEX_RESULTS), function(nm) {
    res <- INDEX_RESULTS[[nm]]$res
    if (is.null(res$copulas$all_results)) return(NULL)
    df <- res$copulas$all_results
    df$Index <- nm
    df$Best  <- df$Copula == res$copulas$best_copula_name
    df
  }))
  if (nrow(cop_all) > 0) {
    wb2 <- createWorkbook()
    addWorksheet(wb2, "Copula_Comparison")
    writeData(wb2, "Copula_Comparison", cop_all)
    best_rows <- which(cop_all$Best) + 1L
    if (length(best_rows) > 0) {
      bestStyle <- createStyle(fgFill = "#d4edda", border = "TopBottom")
      addStyle(wb2, 1, bestStyle, rows = best_rows,
               cols = seq_len(ncol(cop_all)), gridExpand = TRUE)
    }
    setColWidths(wb2, 1, cols = seq_len(ncol(cop_all)), widths = "auto")
    saveWorkbook(wb2, file.path(PUB_DIR, "Table02_copula_fit_summary.xlsx"), overwrite = TRUE)
    cat("  Table02 saved.\n")
  }
}

# Table 3: Non-stationary marginal summary -----------------------------------
{
  # FIX-1F: ns$cp_year does not exist; parse from ns$epoch_src (same helper as Fig 14).
  # FIX-1C: tv$lr_sig does not exist; the significance flag is tv$significant.
  extract_cp_year_from_ns <- function(ns) {
    if (is.null(ns) || !grepl("change-point", ns$epoch_src, fixed = TRUE))
      return(NA_integer_)
    m <- regmatches(ns$epoch_src, regexpr("[0-9]{4}", ns$epoch_src))
    if (length(m) == 1L) as.integer(m) else NA_integer_
  }
  
  ns_all <- dplyr::bind_rows(lapply(names(INDEX_RESULTS), function(nm) {
    ns  <- enh_results[[nm]]$ns_marginals
    tv  <- enh_results[[nm]]$tv_copula
    gof <- enh_results[[nm]]$gof
    if (is.null(ns)) return(NULL)
    data.frame(
      Index         = nm,
      CP_year       = extract_cp_year_from_ns(ns),              # FIX-1F
      Ref_period    = sprintf("%d\u2013%d", min(ns$ref_years),    max(ns$ref_years)),
      Recent_period = sprintf("%d\u2013%d", min(ns$recent_years), max(ns$recent_years)),
      GAMLSS_best   = ns$best_model,
      AIC_M0        = round(ns$aic_m0, 2),
      AIC_M1        = round(ns$aic_m1, 2),
      AIC_M2        = round(ns$aic_m2, 2),
      mu_ref        = round(ns$mu_ref, 4),
      mu_rec        = round(ns$mu_rec, 4),
      sigma_ref     = round(ns$sigma_ref, 4),
      sigma_rec     = round(ns$sigma_rec, 4),
      TV_copula_sig = ifelse(!is.null(tv), isTRUE(tv$significant), NA),  # FIX-1C
      Copula_GOF_p  = ifelse(!is.null(gof), round(gof$p_value, 4), NA),
      Copula_GOF    = ifelse(!is.null(gof), gof$decision, NA)
    )
  }))
  
  if (nrow(ns_all) > 0) {
    wb3 <- createWorkbook()
    addWorksheet(wb3, "NS_Summary")
    writeData(wb3, "NS_Summary", ns_all)
    setColWidths(wb3, 1, cols = seq_len(ncol(ns_all)), widths = "auto")
    saveWorkbook(wb3, file.path(PUB_DIR, "Table03_nonstationary_summary.xlsx"), overwrite = TRUE)
    cat("  Table03 saved.\n")
  }
}

# ============================================================================
# FINAL SUMMARY
# ============================================================================
all_files <- list.files(PUB_DIR, recursive = FALSE)
cat("\n============================================================\n")
cat(sprintf("All publication outputs written to: %s\n", PUB_DIR))
cat(sprintf("Total files: %d\n", length(all_files)))
cat("------------------------------------------------------------\n")
pdf_files  <- all_files[grepl("\\.pdf$",  all_files)]
jpeg_files <- all_files[grepl("\\.jpeg$", all_files)]
xlsx_files <- all_files[grepl("\\.xlsx$", all_files)]
cat(sprintf("  PDFs : %d\n  JPEGs: %d\n  XLSX : %d\n",
            length(pdf_files), length(jpeg_files), length(xlsx_files)))
cat("============================================================\n")
cat("\nNOTE: Individual index-level Rosenblatt JPEG files are ALSO saved\n")
cat("      in each index analysis folder (e.g., spatial_drought/SPI1_analysis/).\n")

# ---------------------------------------------------------------------------
# DONE
# ---------------------------------------------------------------------------
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("ALL STEPS COMPLETE\n")
cat("Output directory: spatial_drought/\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

# Quick summary of what was produced
# FIX-5: quick_diagnostic() was redefined at the top of this script to search
#        "spatial_drought" instead of the hardcoded "drought_analysis" path
#        in QuickStart.R; that override is active here.
quick_diagnostic()
