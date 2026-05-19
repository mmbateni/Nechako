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
#   Step 4 — Load enhancement functions (ext1: GOF, TV-copula; ext2: NS marginals)
#   Step 5 — Run all enhancements on the pipeline results
#   Step 6 — Optional: run seasonal / spatial / cross-index utility analyses
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
cat("STEP 4: Loading enhancement functions (ext1 + ext2)\n")
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
  tv  <- fit_timevarying_copula(res$dc, res$copulas, res$marginals,
                                year_indices, nm, out_dir)
  
  # 5d) Non-stationary severity marginals
  #     ref_years / recent_years = NULL => auto-derived from record length.
  #     Override example:
  #       ns <- fit_nonstationary_marginals(..., ref_years=1960:1990, recent_years=2005:2023)
  ns  <- fit_nonstationary_marginals(res$dc, nm, year_indices, out_dir)
  
  # 5e) Non-stationary SAF curves (reference and recent periods)
  #     Produced for BOTH the conditional (Method A) and Kendall-corrected
  #     (Method B, consistent with De Michele et al. 2026) variants.
  saf_ns_ref <- NULL; saf_ns_rec <- NULL
  saf_ns_ref_k <- NULL; saf_ns_rec_k <- NULL
  if (!is.null(ns)) {
    ref_lbl <- sprintf("reference_%d_%d", min(ns$ref_years), max(ns$ref_years))
    rec_lbl <- sprintf("recent_%d_%d",    min(ns$recent_years), max(ns$recent_years))
    
    saf_ns_ref <- derive_SAF_nonstationary(
      mu_period      = ns$mu_ref,
      sigma_ns       = ns$sigma,
      copula_fit_obj = res$copulas,
      marginal_fits  = res$marginals,
      drought_data   = res$dc,
      index_name     = nm,
      period_label   = ref_lbl,
      output_dir     = out_dir)
    
    saf_ns_rec <- derive_SAF_nonstationary(
      mu_period      = ns$mu_rec,
      sigma_ns       = ns$sigma,
      copula_fit_obj = res$copulas,
      marginal_fits  = res$marginals,
      drought_data   = res$dc,
      index_name     = nm,
      period_label   = rec_lbl,
      output_dir     = out_dir)
    
    # Kendall-corrected non-stationary SAF (Issue 5 fix)
    saf_ns_ref_k <- derive_SAF_nonstationary_kendall(
      mu_period      = ns$mu_ref,
      sigma_ns       = ns$sigma,
      copula_fit_obj = res$copulas,
      marginal_fits  = res$marginals,
      drought_data   = res$dc,
      index_name     = nm,
      period_label   = ref_lbl,
      output_dir     = out_dir)
    
    saf_ns_rec_k <- derive_SAF_nonstationary_kendall(
      mu_period      = ns$mu_rec,
      sigma_ns       = ns$sigma,
      copula_fit_obj = res$copulas,
      marginal_fits  = res$marginals,
      drought_data   = res$dc,
      index_name     = nm,
      period_label   = rec_lbl,
      output_dir     = out_dir)
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
# cat("\n", paste(rep("=", 70), collapse=""), "\n")
# cat("STEP 6: Optional utility analyses\n")
# cat(paste(rep("=", 70), collapse=""), "\n")

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
#     Replace event_years with the actual years of the event of interest
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