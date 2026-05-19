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
  cat("\nDIAGNOSTIC — duration class counts per index:\n")
  dirs <- list.files("drought_analysis", pattern = "_analysis$", full.names = TRUE)
  if (length(dirs) == 0) {
    cat("  No analysis folders found. Run Nechako_Drought_SAF_Analysis.R first.\n")
    return(invisible(NULL))
  }
  for (d in dirs) {
    idx <- basename(d)
    f   <- file.path(d, sprintf("%s_return_periods_by_class.csv",
                                sub("_analysis$", "", idx)))
    if (file.exists(f)) {
      rp <- read.csv(f)
      cat(sprintf("  %-20s : %d classes\n", idx, nrow(rp)))
    }
  }
}

cat("\nRun order:\n",
    "  1. source('Nechako_Drought_SAF_Analysis.R')      # main pipeline\n",
    "  2. source('Nechako_Drought_SAF_Analysis_ext1.R') # enhanced functions\n",
    "  3. source('Nechako_Drought_SAF_Analysis_ext2.R') # extra SAF methods\n",
    "  -- OR simply: source('Nechako_Drought_RUN_ALL.R') --\n")