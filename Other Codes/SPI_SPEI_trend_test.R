####################################################################################
# SPI/SPEI SEASONAL TO TREND RESULTS CONVERTER -
# 
# FIXES:
# 1. Windows-compatible parallel processing (uses parLapply instead of mclapply)
# 2. Proper error handling to continue on failures
# 3. Better diagnostics for troubleshooting
#
# INPUT: spi_XX_monthYY_Mon.csv files (seasonal results)
# OUTPUT: spi_XX_results.csv files (trend analysis results)
####################################################################################

library(data.table)
library(Kendall)

# Parallel processing setup (Windows compatible)
USE_PARALLEL <- TRUE
if (USE_PARALLEL && requireNamespace("parallel", quietly = TRUE)) {
  library(parallel)
  
  # Detect OS
  is_windows <- .Platform$OS.type == "windows"
  
  if (is_windows) {
    N_CORES <- max(1, detectCores() - 1)
    cl <- makeCluster(N_CORES)
    clusterEvalQ(cl, library(Kendall))
    cat(sprintf("✓ Parallel processing enabled (Windows): %d cores\n", N_CORES))
  } else {
    N_CORES <- max(1, detectCores() - 1)
    cl <- NULL
    cat(sprintf("✓ Parallel processing enabled (Unix): %d cores\n", N_CORES))
  }
} else {
  USE_PARALLEL <- FALSE
  N_CORES <- 1
  cl <- NULL
  cat("ℹ Parallel processing disabled\n")
}

setwd("D:/Nechako_Drought/Nechako/")

# Configuration
timescales <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 24)
indices <- c("spi", "spei")
months <- sprintf("%02d", 1:12)
month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

output_dir <- "./temporal_spi_spei/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

####################################################################################
# VECTORIZED HELPER FUNCTIONS
####################################################################################

#' Vectorized Mann-Kendall test for multiple time series
vectorized_mann_kendall <- function(ts_matrix) {
  n_pixels <- nrow(ts_matrix)
  
  # Initialize results
  results <- data.table(
    tau_vc = rep(NA_real_, n_pixels),
    p_value_vc = rep(NA_real_, n_pixels),
    sl_vc = rep(NA_real_, n_pixels),
    filtered_vc = rep(TRUE, n_pixels)
  )
  
  # Function to process one pixel
  process_pixel <- function(i) {
    x <- ts_matrix[i, ]
    x_clean <- x[!is.na(x)]
    
    # Check minimum data and variance
    if (length(x_clean) < 10) {
      return(list(tau = NA, pval = NA, slope = NA, filtered = TRUE))
    }
    
    var_x <- var(x_clean, na.rm = TRUE)
    if (is.na(var_x) || var_x < 1e-6) {
      return(list(tau = NA, pval = NA, slope = NA, filtered = TRUE))
    }
    
    # Mann-Kendall test
    tryCatch({
      mk_result <- MannKendall(x_clean)
      
      # Sen's slope
      slope <- NA_real_
      if (length(x_clean) >= 3) {
        n <- length(x_clean)
        slopes <- numeric((n * (n - 1)) / 2)
        k <- 1
        for (i in 1:(n-1)) {
          for (j in (i+1):n) {
            slopes[k] <- (x_clean[j] - x_clean[i]) / (j - i)
            k <- k + 1
          }
        }
        slope <- median(slopes, na.rm = TRUE)
      }
      
      list(
        tau = as.numeric(mk_result$tau),
        pval = as.numeric(mk_result$sl),
        slope = slope,
        filtered = FALSE
      )
    }, error = function(e) {
      list(tau = NA, pval = NA, slope = NA, filtered = TRUE)
    })
  }
  
  # Process in parallel (Windows-compatible)
  if (USE_PARALLEL && n_pixels > 100) {
    if (is_windows) {
      # Windows: use parLapply with cluster
      pixel_results <- parLapply(cl, 1:n_pixels, process_pixel)
    } else {
      # Unix: use mclapply
      pixel_results <- mclapply(1:n_pixels, process_pixel, mc.cores = N_CORES)
    }
  } else {
    pixel_results <- lapply(1:n_pixels, process_pixel)
  }
  
  # Extract results
  results$tau_vc <- sapply(pixel_results, function(x) x$tau)
  results$p_value_vc <- sapply(pixel_results, function(x) x$pval)
  results$sl_vc <- sapply(pixel_results, function(x) x$slope)
  results$filtered_vc <- sapply(pixel_results, function(x) x$filtered)
  
  return(results)
}

#' Vectorized event detection
vectorized_event_detection <- function(ts_matrix, threshold = -1.0) {
  n_pixels <- nrow(ts_matrix)
  
  results <- data.table(
    n_events = rep(NA_real_, n_pixels),
    mean_duration = rep(NA_real_, n_pixels),
    max_intensity = rep(NA_real_, n_pixels)
  )
  
  for (i in 1:n_pixels) {
    x <- ts_matrix[i, ]
    x_clean <- x[!is.na(x)]
    
    if (length(x_clean) == 0) next
    
    # Detect events
    in_event <- x_clean < threshold
    if (!any(in_event)) {
      results$n_events[i] <- 0
      results$mean_duration[i] <- 0
      results$max_intensity[i] <- 0
      next
    }
    
    event_starts <- which(diff(c(FALSE, in_event)) == 1)
    event_ends <- which(diff(c(in_event, FALSE)) == -1)
    
    n_events <- length(event_starts)
    durations <- event_ends - event_starts + 1
    
    intensities <- sapply(seq_along(event_starts), function(j) {
      min(x_clean[event_starts[j]:event_ends[j]])
    })
    
    results$n_events[i] <- n_events
    results$mean_duration[i] <- mean(durations)
    results$max_intensity[i] <- min(intensities)
  }
  
  return(results)
}

#' Vectorized regime shift detection
vectorized_regime_shift <- function(ts_matrix, year_vec) {
  n_pixels <- nrow(ts_matrix)
  regime_years <- rep(NA_real_, n_pixels)
  
  for (i in 1:n_pixels) {
    x <- ts_matrix[i, ]
    valid_idx <- !is.na(x)
    
    if (sum(valid_idx) < 20) next
    
    x_clean <- x[valid_idx]
    years_clean <- year_vec[valid_idx]
    
    # Change point detection
    x_centered <- x_clean - mean(x_clean)
    cumsum_x <- cumsum(x_centered)
    change_point <- which.max(abs(cumsum_x))
    
    regime_years[i] <- years_clean[change_point]
  }
  
  return(regime_years)
}

#' Vectorized runs test
vectorized_runs_test <- function(ts_matrix, threshold = 0) {
  n_pixels <- nrow(ts_matrix)
  
  results <- data.table(
    p_value_runs = rep(NA_real_, n_pixels),
    clustering = rep(NA_character_, n_pixels),
    filtered_runs = rep(TRUE, n_pixels)
  )
  
  for (i in 1:n_pixels) {
    x <- ts_matrix[i, ]
    x_clean <- x[!is.na(x)]
    
    if (length(x_clean) < 10) next
    
    # Binary series
    binary <- ifelse(x_clean < threshold, 0, 1)
    
    # Count runs
    runs <- rle(binary)
    n_runs <- length(runs$lengths)
    n <- length(binary)
    n1 <- sum(binary == 1)
    n0 <- sum(binary == 0)
    
    if (n1 == 0 || n0 == 0) next
    
    # Expected runs
    expected_runs <- (2 * n0 * n1) / n + 1
    var_runs <- (2 * n0 * n1 * (2 * n0 * n1 - n)) / (n^2 * (n - 1))
    
    if (var_runs <= 0) next
    
    # Z-score and p-value
    z <- (n_runs - expected_runs) / sqrt(var_runs)
    p_value <- 2 * pnorm(-abs(z))
    
    clustering <- if (n_runs < expected_runs) "clustered" else 
      if (n_runs > expected_runs) "dispersed" else "random"
    
    results$p_value_runs[i] <- p_value
    results$clustering[i] <- clustering
    results$filtered_runs[i] <- FALSE
  }
  
  return(results)
}

####################################################################################
# MAIN PROCESSING FUNCTION
####################################################################################

process_seasonal_to_trends_optimized <- function(index_type, timescale) {
  cat(sprintf("\n╔══════════════════════════════════════════════════╗\n"))
  cat(sprintf("║  Processing %s-%02d                              \n", toupper(index_type), timescale))
  cat(sprintf("╚══════════════════════════════════════════════════╝\n"))
  
  input_dir <- sprintf("./%s_results_seasonal/", index_type)  
  if (!dir.exists(input_dir)) {
    cat(sprintf("❌ Directory not found: %s\n", input_dir))
    return(NULL)
  }
  
  # STEP 1: Load all monthly files ONCE
  cat("\n[1/5] Loading monthly files...\n")
  all_months_data <- vector("list", 12)
  
  start_time <- Sys.time()
  
  for (m in 1:12) {
    month_file <- sprintf("%s%s_%02d_month%s_%s.csv", 
                          input_dir, index_type, timescale, 
                          months[m], month_names[m])
    
    if (!file.exists(month_file)) {
      cat(sprintf("  ⚠️  Missing: %s\n", month_names[m]))
      next
    }
    
    all_months_data[[m]] <- fread(month_file, showProgress = FALSE)
    cat(sprintf("  ✓ %s: %d rows\n", month_names[m], nrow(all_months_data[[m]])))
  }
  
  load_time <- Sys.time() - start_time
  cat(sprintf("  Time: %.1f seconds\n", as.numeric(load_time, units = "secs")))
  
  # Validate
  n_files_loaded <- sum(!sapply(all_months_data, is.null))
  if (n_files_loaded == 0) {
    cat("❌ No monthly files loaded\n")
    return(NULL)
  }
  
  # STEP 2: Extract metadata and setup
  cat("\n[2/5] Setting up spatial grid...\n")
  ref_data <- all_months_data[[which(!sapply(all_months_data, is.null))[1]]]
  
  year_cols <- grep("^[0-9]{4}$", names(ref_data), value = TRUE)
  years <- as.integer(year_cols)
  n_years <- length(years)
  n_pixels <- nrow(ref_data)
  n_months <- n_years * 12
  
  # IDENTIFY VALID PIXELS (not all NA across all months)
  cat("  Identifying valid pixels...\n")
  valid_pixels <- rep(FALSE, n_pixels)
  for (i in 1:n_pixels) {
    has_data <- FALSE
    for (m in 1:12) {
      if (!is.null(all_months_data[[m]])) {
        month_data <- all_months_data[[m]][i, ]
        year_values <- unlist(month_data[, ..year_cols])
        if (sum(!is.na(year_values)) > 0) {
          has_data <- TRUE
          break
        }
      }
    }
    valid_pixels[i] <- has_data
  }
  
  n_valid <- sum(valid_pixels)
  cat(sprintf("  Valid pixels: %d (%.1f%%)\n", n_valid, n_valid/n_pixels*100))
  
  # Filter to valid pixels only
  ref_data <- ref_data[valid_pixels, ]
  for (m in 1:12) {
    if (!is.null(all_months_data[[m]])) {
      all_months_data[[m]] <- all_months_data[[m]][valid_pixels, ]
    }
  }
  n_pixels <- n_valid  # Update pixel count
  
  cat(sprintf("  Pixels: %d\n", n_pixels))
  cat(sprintf("  Years: %d to %d (%d years)\n", min(years), max(years), n_years))
  cat(sprintf("  Total time steps: %d months\n", n_months))
  
  # STEP 3: Build time series matrix (VECTORIZED!)
  cat("\n[3/5] Building time series matrix...\n")
  start_time <- Sys.time()
  
  # Pre-allocate matrix: rows = pixels, cols = time steps
  ts_matrix <- matrix(NA_real_, nrow = n_pixels, ncol = n_months)
  
  # Fill matrix efficiently
  time_idx <- 1
  for (year in years) {
    year_str <- as.character(year)
    
    for (m in 1:12) {
      if (!is.null(all_months_data[[m]]) && year_str %in% names(all_months_data[[m]])) {
        # Extract entire column at once (vectorized!)
        ts_matrix[, time_idx] <- all_months_data[[m]][[year_str]]
      }
      time_idx <- time_idx + 1
    }
  }
  
  build_time <- Sys.time() - start_time
  cat(sprintf("  Matrix size: %d × %d\n", nrow(ts_matrix), ncol(ts_matrix)))
  cat(sprintf("  Non-NA values: %.1f%%\n", 
              sum(!is.na(ts_matrix)) / length(ts_matrix) * 100))
  cat(sprintf("  Time: %.1f seconds\n", as.numeric(build_time, units = "secs")))
  
  # STEP 4: Perform all analyses (VECTORIZED!)
  cat("\n[4/5] Computing statistics...\n")
  start_time <- Sys.time()
  
  # Initialize results with spatial coordinates
  results <- data.table(
    lon = ref_data$lon,
    lat = ref_data$lat
  )
  
  # Trend analysis (vectorized)
  cat("  → Mann-Kendall tests...")
  trend_results <- tryCatch({
    vectorized_mann_kendall(ts_matrix)
  }, error = function(e) {
    cat(sprintf("\n  ❌ ERROR: %s\n", e$message))
    data.table(
      tau_vc = rep(NA_real_, n_pixels),
      p_value_vc = rep(NA_real_, n_pixels),
      sl_vc = rep(NA_real_, n_pixels),
      filtered_vc = rep(TRUE, n_pixels)
    )
  })
  results <- cbind(results, trend_results)
  cat(" Done\n")
  
  # Event characteristics (vectorized)
  cat("  → Event detection...")
  event_results <- vectorized_event_detection(ts_matrix)
  results <- cbind(results, event_results)
  cat(" Done\n")
  
  # Regime shifts (vectorized)
  cat("  → Regime shift detection...")
  year_vec <- rep(years, each = 12) + (rep(1:12, times = n_years) - 1) / 12
  results$regime_shift_year <- vectorized_regime_shift(ts_matrix, year_vec)
  cat(" Done\n")
  
  # Runs test (vectorized)
  cat("  → Temporal clustering...")
  runs_results <- vectorized_runs_test(ts_matrix)
  results <- cbind(results, runs_results)
  cat(" Done\n")
  
  # Add spectral peaks placeholder
  results$n_spectral_peaks <- 0
  
  stats_time <- Sys.time() - start_time
  cat(sprintf("  Time: %.1f seconds\n", as.numeric(stats_time, units = "secs")))
  
  # STEP 5: Add time series and save
  cat("\n[5/5] Saving results...\n")
  start_time <- Sys.time()
  
  # Add time series columns
  for (t in 1:n_months) {
    results[[sprintf("time_%03d", t)]] <- ts_matrix[, t]
  }
  
  # Summary statistics
  cat("\n  📊 Results Summary:\n")
  cat("  ─────────────────────────────────────────────────\n")
  cat(sprintf("     Total pixels: %d\n", n_pixels))
  cat(sprintf("     Valid tau: %d (%.1f%%)\n", 
              sum(!is.na(results$tau_vc)),
              sum(!is.na(results$tau_vc)) / n_pixels * 100))
  cat(sprintf("     Valid p-values: %d (%.1f%%)\n", 
              sum(!is.na(results$p_value_vc)),
              sum(!is.na(results$p_value_vc)) / n_pixels * 100))
  cat(sprintf("     Passed variance check: %d (%.1f%%)\n", 
              sum(!results$filtered_vc),
              sum(!results$filtered_vc) / n_pixels * 100))
  
  if (sum(!is.na(results$tau_vc)) > 0) {
    cat(sprintf("     Median tau: %.4f\n", 
                median(results$tau_vc, na.rm = TRUE)))
    cat(sprintf("     Significant trends (p<0.05): %d (%.1f%%)\n",
                sum(results$p_value_vc < 0.05, na.rm = TRUE),
                sum(results$p_value_vc < 0.05, na.rm = TRUE) / sum(!is.na(results$p_value_vc)) * 100))
  }
  
  # Save
  output_file <- sprintf("%s%s_%02d_results.csv", output_dir, index_type, timescale)
  fwrite(results, output_file, showProgress = FALSE)
  
  file_size <- file.info(output_file)$size / 1024 / 1024
  save_time <- Sys.time() - start_time
  
  cat(sprintf("\n  ✅ Saved: %s (%.2f MB)\n", basename(output_file), file_size))
  cat(sprintf("  Time: %.1f seconds\n", as.numeric(save_time, units = "secs")))
  
  return(results)
}

####################################################################################
# RUN PROCESSING FOR ALL TIMESCALES
####################################################################################

cat("\n╔════════════════════════════════════════════════════════╗\n")
cat("║  SPI/SPEI SEASONAL → TREND CONVERTER (WINDOWS FIXED)   ║\n")
cat("╚════════════════════════════════════════════════════════╝\n\n")

cat("PERFORMANCE FEATURES:\n")
cat("  ✓ Each file read exactly once\n")
cat("  ✓ Vectorized matrix operations\n")
cat("  ✓ Efficient data.table processing\n")
if (USE_PARALLEL) {
  if (is_windows) {
    cat(sprintf("  ✓ Windows-compatible parallel processing (%d cores)\n", N_CORES))
  } else {
    cat(sprintf("  ✓ Parallel processing (%d cores)\n", N_CORES))
  }
}
cat("\n")

total_start <- Sys.time()
all_results <- list()
success_count <- 0
fail_count <- 0

for (idx in indices) {
  cat(sprintf("\n█████████████████████████████████████████████████████████\n"))
  cat(sprintf("█  %s\n", toupper(idx)))
  cat(sprintf("█████████████████████████████████████████████████████████\n"))
  
  for (scale in timescales) {
    scale_start <- Sys.time()
    
    result <- tryCatch({
      process_seasonal_to_trends_optimized(idx, scale)
    }, error = function(e) {
      cat(sprintf("\n  ❌ ERROR: %s\n", e$message))
      NULL
    })
    
    scale_time <- Sys.time() - scale_start
    
    if (!is.null(result)) {
      all_results[[sprintf("%s_%02d", idx, scale)]] <- result
      success_count <- success_count + 1
      cat(sprintf("\n  ✅ Total time for %s-%02d: %.1f seconds\n", 
                  toupper(idx), scale, as.numeric(scale_time, units = "secs")))
    } else {
      fail_count <- fail_count + 1
    }
  }
}

# Cleanup cluster if Windows
if (USE_PARALLEL && is_windows && !is.null(cl)) {
  stopCluster(cl)
  cat("\n✓ Cluster stopped\n")
}

total_time <- Sys.time() - total_start

cat("\n\n╔════════════════════════════════════════════════════════╗\n")
cat("║  PROCESSING COMPLETE                                   ║\n")
cat("╚════════════════════════════════════════════════════════╝\n\n")

cat("📊 SUMMARY:\n")
cat("═══════════════════════════════════════════════════════\n")
cat(sprintf("  Successful: %d\n", success_count))
cat(sprintf("  Failed: %d\n", fail_count))
cat(sprintf("  Total time: %.1f minutes\n", as.numeric(total_time, units = "mins")))
if (success_count > 0) {
  cat(sprintf("  Average per file: %.1f seconds\n", 
              as.numeric(total_time, units = "secs") / success_count))
}

cat("\n📁 OUTPUT FILES:\n")
cat("═══════════════════════════════════════════════════════\n")
output_files <- list.files(output_dir, pattern = "_results\\.csv$", full.names = FALSE)

if (length(output_files) > 0) {
  total_size <- 0
  for (f in sort(output_files)) {
    file_path <- file.path(output_dir, f)
    file_size <- file.info(file_path)$size / 1024 / 1024
    total_size <- total_size + file_size
    cat(sprintf("  ✓ %-30s (%.2f MB)\n", f, file_size))
  }
  cat(sprintf("\n  Total output size: %.2f MB\n", total_size))
} else {
  cat("  ⚠️  No output files created\n")
}

cat("\n✅ Ready for visualization!\n")
cat("   Next: Run SPI_SPEI_trends_visualization.R\n\n")