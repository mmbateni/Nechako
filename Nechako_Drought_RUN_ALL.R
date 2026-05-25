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
#   Step 4 — Load enhancement functions (ext1 v4: GOF, TV-copula;
#                                        ext2 v4: NS marginals)
#   Step 5 — Run all enhancements on the pipeline results
#   Step 6 — Optional: run seasonal / spatial / cross-index utility analyses
#
# v6 BUG FIXES (ext1):
#   [ext1] detect_copula_changepoints(): rolling Kendall tau loop now pre-checks
#          sd(severity) and sd(area_pct) before calling cor(); zero-variance windows
#          are silently set to NA and counted in a single aggregate log message
#          instead of emitting ≈9 "standard deviation is zero" warnings for SPI-3.
#
#   [QuickStart] total_events now sourced from nrow(ev_file) when the raw event
#          CSV exists (authoritative); sum(rp$n_events) retained as fallback and
#          cross-check target.  MISMATCH flag compares the two counts.
#
# v4 BUG FIXES (ext1 + ext2):
#   [ext1] rosenblatt_pit_gof(): e2 was computed as C(u_sev | u_area)
#          (wrong direction) via cCopula(..., indices=1), producing e2 ≡ e1,
#          τ = 1.000, and 0 % valid values for Plackett.  Fixed: e2 is now
#          ∂C(u1,u2)/∂u1 = C(u_area | u_severity) via finite-difference of
#          pCopula (.h_func_fd helper), which works for all copula families.
#   [ext2] fit_nonstationary_marginals(): yr_match was year_vec[1:nrow(dc)]
#          (first N months of monthly record) instead of dc$year (actual
#          observation years).  This truncated recent_years to ~1998 and
#          falsely classified 1996 as a near-edge change point for SPEI-1.
#          Fixed: yr_match <- dc$year.  Recent epochs now extend to 2025.
#============================================================================

setwd("D:/Nechako_Drought/Nechako")

# ---------------------------------------------------------------------------
# STEP 1 — Input verification
# ---------------------------------------------------------------------------
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("STEP 1: Input verification\n")
cat(paste(rep("=", 70), collapse=""), "\n")
source("Nechako_Drought_QuickStart.R")

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
cat("STEP 4: Loading enhancement functions (ext1 v6 + ext2 v4)\n")
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
  SPI1  = list(res = res_spi1,  dir = "drought_analysis/SPI1_analysis",  scale = 1),
  SPEI1 = list(res = res_spei1, dir = "drought_analysis/SPEI1_analysis", scale = 1),
  SPI3  = list(res = res_spi3,  dir = "drought_analysis/SPI3_analysis",  scale = 3),
  SPEI3 = list(res = res_spei3, dir = "drought_analysis/SPEI3_analysis", scale = 3)
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
  #     evaluated at the epoch midpoint year using the fitted TV copula.
  saf_ns_ref <- NULL; saf_ns_rec <- NULL
  saf_ns_ref_k <- NULL; saf_ns_rec_k <- NULL
  if (!is.null(ns)) {
    ref_lbl <- sprintf("reference_%d_%d", min(ns$ref_years), max(ns$ref_years))
    rec_lbl <- sprintf("recent_%d_%d",    min(ns$recent_years), max(ns$recent_years))
    
    saf_ns_ref <- derive_SAF_nonstationary(
      mu_period        = ns$mu_ref,
      sigma_period     = ns$sigma_ref,
      copula_fit_obj   = res$copulas,
      marginal_fits    = res$marginals,
      drought_data     = res$dc,
      index_name       = nm,
      period_label     = ref_lbl,
      output_dir       = out_dir,
      tv_copula_result = tv)           # period-specific copula if TV significant
    
    saf_ns_rec <- derive_SAF_nonstationary(
      mu_period        = ns$mu_rec,
      sigma_period     = ns$sigma_rec,
      copula_fit_obj   = res$copulas,
      marginal_fits    = res$marginals,
      drought_data     = res$dc,
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
      index_name       = nm,
      period_label     = rec_lbl,
      output_dir       = out_dir,
      tv_copula_result = tv)
  }
  
  enh_results[[nm]] <- list(
    dependency   = dep,
    gof          = gof,
    tv_copula    = tv,
    ns_marginals = ns,
    saf_ns_ref   = saf_ns_ref,
    saf_ns_rec   = saf_ns_rec,
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
# seasonal_drought_analysis(res_spi1$dc, "SPI1", "drought_analysis/SPI1_analysis")
# seasonal_drought_analysis(res_spi3$dc, "SPI3", "drought_analysis/SPI3_analysis")

# 6c) Pixel-level drought frequency and severity maps
# map_drought_frequency(spi1,  index_name = "SPI1",  output_dir = "drought_maps")
# map_drought_frequency(spei1, index_name = "SPEI1", output_dir = "drought_maps")

# 6d) Cross-index correlation (SPI-1 vs SPI-3 severity)
# cross_correlate_indices(res_spi1$dc, res_spi3$dc, output_dir = "drought_analysis")

# 6e) Place a specific historical event on the SAF surface
# Replace event_years with the actual years of the event of interest
# place_event_on_saf(res_spi1$dc, res_spi1$copulas, res_spi1$marginals,
#                    event_years = c(2003, 2004), index_name = "SPI1",
#                    output_dir  = "drought_analysis/SPI1_analysis")

# ---------------------------------------------------------------------------
# DONE
# ---------------------------------------------------------------------------
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("ALL STEPS COMPLETE\n")
cat("Output directory: drought_analysis/\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

# Quick summary of what was produced
quick_diagnostic()