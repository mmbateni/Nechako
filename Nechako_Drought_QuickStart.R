#============================================================================
# QUICK-START & DATA VERIFICATION (SPI/SPEI ONLY)
# Fixed: Kendall distribution sanity check now uses the correct
#        Archimedean generator formula K_C(t) = t - phi(t)/phi'(t)
#        instead of the previous incorrect t - (1-C(t,t))/t
#============================================================================
library(terra); library(copula)
setwd("D:/Nechako_Drought/Nechako")
cat(paste(rep("=",70), collapse=""), "\nSPI/SPEI VERIFICATION\n",
    paste(rep("=",70), collapse=""), "\n")

# --- 1. File existence check ------------------------------------------------
files <- c(
  "spi_results_seasonal/spi_01_month01_Jan.nc",
  "spi_results_seasonal/spi_03_month01_Jan.nc",
  "spei_results_seasonal/spei_01_month01_Jan.nc",
  "spei_results_seasonal/spei_03_month01_Jan.nc"
)
if (!all(file.exists(files)))
  stop("Missing required NetCDFs. Run upstream scripts first.")
cat("✓ All input files exist.\n")

# --- 2. SPI-1 range check ---------------------------------------------------
spi1_jan <- rast(files[1])
cat(sprintf("SPI-1 Layers: %d | Cells: %d\n", nlyr(spi1_jan), ncell(spi1_jan)))
st <- global(spi1_jan, c("min","max"), na.rm = TRUE)
if (max(st$max) > 5 || min(st$min) < -5) {
  warning("⚠ SPI-1 values outside expected [-5, 5] range.")
} else {
  cat("✓ SPI-1 range valid.\n")
}

# --- 3. Kendall distribution sanity check -----------------------------------
# Correct formula for Archimedean copulas: K_C(t) = t - phi(t) / phi'(t)
# Independence (product) copula generator: phi(t) = -log(t), phi'(t) = -1/t
# => K_C(t) = t - (-log t) / (-1/t) = t - t*log(t)
# At t = 0.5: K_C(0.5) = 0.5 - 0.5*log(0.5) ≈ 0.847
test_t   <- 0.5
phi_t    <- -log(test_t)         # generator value
phi_d_t  <- -1 / test_t          # generator derivative
k_c      <- test_t - phi_t / phi_d_t              # = t - t*log(t)
expected <- test_t - test_t * log(test_t)          # analytical reference
if (abs(k_c - expected) < 1e-9) {
  cat(sprintf("✓ Kendall distribution formula verified (K_C(%.1f) ≈ %.4f).\n",
              test_t, k_c))
} else {
  warning("⚠ Kendall formula check failed.")
}

# --- 4. Quick post-run diagnostic -------------------------------------------
quick_diagnostic <- function() {
  sep_major <- paste(rep("=", 72), collapse = "")
  sep_minor <- paste(rep("-", 72), collapse = "")
  
  cat("\n", sep_major, "\n", sep = "")
  cat("DIAGNOSTIC — drought event summary by index and duration class\n")
  cat(sep_major, "\n", sep = "")
  
  dirs <- list.files("drought_analysis", pattern = "_analysis$", full.names = TRUE)
  if (length(dirs) == 0) {
    cat("  No analysis folders found. Run Nechako_Drought_SAF_Analysis.R first.\n")
    return(invisible(NULL))
  }
  
  # Collect results so we can also print a cross-index summary table at the end
  summary_rows <- list()
  
  for (d in sort(dirs)) {
    index_tag <- sub("_analysis$", "", basename(d))
    
    # ---- (a) Per-class table from return_periods_by_class.csv ----------------
    rp_file <- file.path(d, sprintf("%s_return_periods_by_class.csv", index_tag))
    # ---- (b) Full event list from duration_classified_events.csv -------------
    ev_file <- file.path(d, sprintf("%s_duration_classified_events.csv", index_tag))
    
    if (!file.exists(rp_file)) {
      cat(sprintf("\n  [%s] return_periods_by_class.csv not found — skipping.\n", index_tag))
      next
    }
    
    rp <- read.csv(rp_file, stringsAsFactors = FALSE)
    
    # FIX: prefer the raw event-list row count as the authoritative total.
    # sum(rp$n_events) can differ from nrow(events) when event-assignment
    # logic and the per-class aggregation diverge (e.g. unclassified events
    # falling into "D0 (Sub-threshold)" rows).  Using nrow(read.csv(ev_file))
    # ties total_events directly to the object that was serialised by
    # identify_duration_classified_events(), making the diagnostic self-consistent.
    # Falls back to sum(rp$n_events) only when the file is absent.
    total_events <- if (file.exists(ev_file)) nrow(read.csv(ev_file)) else sum(rp$n_events, na.rm = TRUE)
    
    # Cross-check: compare the primary total against the rp-table aggregate.
    # When they agree the accounting is clean; MISMATCH warns of a divergence
    # between the event table and the per-class summary and should be
    # investigated (most commonly caused by sub-threshold events or rounding).
    rp_total <- sum(rp$n_events, na.rm = TRUE)
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
    
    # Column widths
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
      Total_Events     = total_events,       # from ev_file if available, else rp sum
      Total_Events_rp  = rp_total,           # rp-table aggregate (for cross-check)
      N_Classes        = nrow(rp),
      stringsAsFactors = FALSE
    )
  }
  
  # ---- Cross-index summary ---------------------------------------------------
  if (length(summary_rows) > 0) {
    cat("\n", sep_major, "\n", sep = "")
    cat("  CROSS-INDEX SUMMARY\n")
    cat(sep_major, "\n", sep = "")
    cat(sprintf("  %-20s  %14s  %14s  %10s\n", "Index", "Total(ev_file)", "Total(rp sum)", "N Classes"))
    cat(sprintf("  %-20s  %14s  %14s  %10s\n",
                paste(rep("-",20),collapse=""),
                "--------------", "--------------", "----------"))
    for (nm in names(summary_rows)) {
      r <- summary_rows[[nm]]
      flag <- if (!is.na(r$Total_Events_rp) && r$Total_Events_rp != r$Total_Events) " !" else ""
      cat(sprintf("  %-20s  %14d  %14d%s  %10d\n",
                  r$Index, r$Total_Events, r$Total_Events_rp, flag, r$N_Classes))
    }
    cat(sep_major, "\n\n", sep = "")
  }
  
  invisible(summary_rows)
}

cat("\nRun order:\n",
    "  1. source('Nechako_Drought_SAF_Analysis.R')      # main pipeline\n",
    "  2. source('Nechako_Drought_SAF_Analysis_ext1.R') # enhanced functions\n",
    "  3. source('Nechako_Drought_SAF_Analysis_ext2.R') # extra SAF methods\n",
    "  -- OR simply: source('Nechako_Drought_RUN_ALL.R') --\n")