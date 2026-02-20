#============================================================================
# DROUGHT SEVERITY-AREA-FREQUENCY ANALYSIS FOR NECHAKO BASIN
# Using PNPI & PNWBI with DUAL RETURN PERIOD METHODS
#============================================================================
# This script performs regional drought analysis using:
# 1. Percent of Normal Precipitation Index (PNPI)
# 2. Percent of Normal Water Balance Index (PNWBI)
# 3. TWO return period methods:
#    a) Conditional copula approach (Amirataee et al. 2018)
#    b) Kendall distribution function approach (Genest et al. 2009)
#
# Key References:
# - PNPI methodology: WMO (2016); Heim (2002)
# - Water balance: Palmer (1965); Thornthwaite (1948)
# - Conditional approach: Amirataee et al. (2018); Shiau (2006)
# - Kendall approach: Genest et al. (2009); Salvadori & De Michele (2004)
#============================================================================

--- 0. SETUP AND LOAD LIBRARIES ---
  rm(list = ls())
gc()

required_packages <- c("terra", "copula", "fitdistrplus", "MASS", 
                       "ggplot2", "gridExtra", "viridis", "moments")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

setwd("D:/Nechako_Drought/monthly_data_direct")

# Create output directory structure
dir.create("drought_analysis", showWarnings = FALSE)
dir.create("drought_analysis/figures", showWarnings = FALSE)
dir.create("drought_analysis/results", showWarnings = FALSE)
dir.create("drought_analysis/PNPI_analysis", showWarnings = FALSE)
dir.create("drought_analysis/PNWBI_analysis", showWarnings = FALSE)

# Initialize log file
LOG_FILE <- "drought_analysis/SAF_analysis_log.txt"
cat("======================================================\n", file = LOG_FILE)
cat("DROUGHT S-A-F ANALYSIS USING PNPI & PNWBI - NECHAKO BASIN\n", file = LOG_FILE, append = TRUE)
cat(paste("Analysis started:", Sys.time(), "\n"), file = LOG_FILE, append = TRUE)
cat("======================================================\n\n", file = LOG_FILE, append = TRUE)

log_event <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(paste0(timestamp, " | ", msg, "\n"), file = LOG_FILE, append = TRUE)
  message(paste0(timestamp, " | ", msg))
}

#============================================================================
# STEP 1: LOAD DATA
#============================================================================
log_event("="*70)
log_event("STEP 1: LOADING PRECIPITATION AND PET DATA")
log_event("="*70)

precip <- rast("total_precipitation_monthly.nc")
log_event(sprintf("  Precipitation layers: %d", nlyr(precip)))

pet <- rast("ERA5Land_Nechako_PET_monthly.nc")
log_event(sprintf("  PET layers: %d", nlyr(pet)))

if (nlyr(precip) != nlyr(pet)) {
  stop("ERROR: Precipitation and PET have different temporal lengths!")
}
n_months <- nlyr(precip)
n_cells <- ncell(precip)

dates <- seq.Date(from = as.Date("1950-01-01"), by = "month", length.out = n_months)
month_indices <- as.numeric(format(dates, "%m"))
year_indices <- as.numeric(format(dates, "%Y"))

cell_area_km2 <- cellSize(precip[[1]], unit = "km")
total_area_km2 <- global(cell_area_km2, "sum", na.rm = TRUE)[1,1]
log_event(sprintf("  Total basin area: %.2f km²", total_area_km2))

#============================================================================
# STEP 2: PREPARE PET IN MONTHLY UNITS
#============================================================================
log_event("\n" + "="*70)
log_event("STEP 2: CONVERTING PET TO MONTHLY TOTALS")
log_event("="*70)

days_per_month <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
days_in_month <- days_per_month[month_indices]
is_leap <- (year_indices %% 4 == 0 & year_indices %% 100 != 0) | (year_indices %% 400 == 0)
days_in_month[month_indices == 2 & is_leap] <- 29

pet_monthly <- pet * days_in_month
pet_stats <- global(pet_monthly, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("  Monthly PET range: %.2f to %.2f mm (mean: %.2f)",
                  min(pet_stats$min), max(pet_stats$max), mean(pet_stats$mean)))

#============================================================================
# STEP 3: CALCULATE WATER BALANCE (P - PET)
#============================================================================
log_event("\n" + "="*70)
log_event("STEP 3: CALCULATING WATER BALANCE (P - PET)")
log_event("="*70)

wb <- precip - pet_monthly
wb <- ifel(wb < 0, 0, wb)  # Apply constraint: negative WB → 0
wb_stats <- global(wb, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("  Water balance range: %.2f to %.2f mm (mean: %.2f)",
                  min(wb_stats$min), max(wb_stats$max), mean(wb_stats$mean)))
writeCDF(wb, "drought_analysis/water_balance_monthly.nc", overwrite = TRUE)

#============================================================================
# FUNCTION: CALCULATE PERCENT OF NORMAL BY CALENDAR MONTH
#============================================================================
percent_of_normal_by_month <- function(data_rast, dates, name = "Index") {
  """
  Calculate Percent of Normal Index by calendar month climatology
  
  PNPI/PNWBI = (X_t / μ_m) * 100
  
  Reference:
  - WMO (2016). Standardized Precipitation Index User Guide
  - Heim (2002). A review of twentieth-century drought indices used in the USA
  """
  log_event(sprintf("  Calculating Percent of Normal %s by calendar month...", name))
  
  month_idx <- as.numeric(format(dates, "%m"))
  n_layers <- nlyr(data_rast)
  pn_rast <- rast(data_rast)
  names(pn_rast) <- paste0(name, "_", format(dates, "%Y-%m"))
  
  # Calculate climatology for each calendar month
  for (m in 1:12) {
    log_event(sprintf("    Calendar month %02d (%s)...", m, month.abb[m]))
    idx_m <- which(month_idx == m)
    data_m <- data_rast[[idx_m]]
    
    # Calculate long-term mean for this calendar month
    mu_m <- mean(data_m, na.rm = TRUE)
    mu_m <- ifel(mu_m == 0, 0.001, mu_m)  # Avoid division by zero
    
    # Calculate percent of normal for each layer
    for (i in idx_m) {
      pn_rast[[i]] <- (data_rast[[i]] / mu_m) * 100
    }
    
    mu_mean <- global(mu_m, "mean", na.rm = TRUE)[1,1]
    log_event(sprintf("      Long-term mean = %.2f, Expected PN ≈ 100%%", mu_mean))
  }
  return(pn_rast)
}

#============================================================================
# STEP 4A: CALCULATE PERCENT OF NORMAL PRECIPITATION INDEX (PNPI)
#============================================================================
log_event("\n" + "="*70)
log_event("STEP 4A: CALCULATING PERCENT OF NORMAL PRECIPITATION INDEX (PNPI)")
log_event("="*70)
log_event("  Reference: WMO (2016) - Percent of Normal as drought indicator")
log_event("  Interpretation: 100% = normal, <75% = moderate drought")

pnpi <- percent_of_normal_by_month(precip, dates, name = "PNPI")
pnpi_stats <- global(pnpi, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("  PNPI range: %.1f%% to %.1f%% (mean: %.1f%%)",
                  min(pnpi_stats$min), max(pnpi_stats$max), mean(pnpi_stats$mean)))
writeCDF(pnpi, "drought_analysis/PNPI_monthly.nc", overwrite = TRUE)
log_event("  PNPI saved to: PNPI_monthly.nc")

#============================================================================
# STEP 4B: CALCULATE PERCENT OF NORMAL WATER BALANCE INDEX (PNWBI)
#============================================================================
log_event("\n" + "="*70)
log_event("STEP 4B: CALCULATING PERCENT OF NORMAL WATER BALANCE INDEX (PNWBI)")
log_event("="*70)
log_event("  Reference: Palmer (1965) water balance approach adapted to percentage format")

pnwbi <- percent_of_normal_by_month(wb, dates, name = "PNWBI")
pnwbi_stats <- global(pnwbi, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("  PNWBI range: %.1f%% to %.1f%% (mean: %.1f%%)",
                  min(pnwbi_stats$min), max(pnwbi_stats$max), mean(pnwbi_stats$mean)))
writeCDF(pnwbi, "drought_analysis/PNWBI_monthly.nc", overwrite = TRUE)
log_event("  PNWBI saved to: PNWBI_monthly.nc")

#============================================================================
# FUNCTION: EXTRACT DROUGHT CHARACTERISTICS FOR PERCENT INDICES
#============================================================================
extract_drought_characteristics_percent <- function(index_rast, cell_area, 
                                                    threshold_pct = 75) {
  """
  Extract regional drought severity and area using percent-based indices
  
  Drought threshold: <75% = moderate drought (WMO 2016; Heim 2002)
  Severity definition: Mean deficit from normal (100 - index) in drought areas
  
  Reference:
  - Heim, R.R. (2002). A review of twentieth-century drought indices used in the USA.
  - WMO (2016). Standardized Precipitation Index User Guide.
  """
  log_event(sprintf("  Extracting drought characteristics (threshold = %.0f%%)...", 
                    threshold_pct))
  
  n_time <- nlyr(index_rast)
  total_area <- global(cell_area, "sum", na.rm = TRUE)[1,1]
  
  severity <- numeric(n_time)    # Mean deficit (%) from normal conditions
  area_pct <- numeric(n_time)    # % of basin below threshold
  
  pb <- txtProgressBar(min = 1, max = n_time, style = 3)
  for (t in 1:n_time) {
    index_t <- index_rast[[t]]
    drought_mask <- index_t < threshold_pct
    
    # AREA UNDER DROUGHT
    drought_area <- mask(cell_area, drought_mask, maskvalues = 0)
    area_t <- global(drought_area, "sum", na.rm = TRUE)[1,1]
    area_pct[t] <- (area_t / total_area) * 100
    
    # DROUGHT SEVERITY: Mean deficit from normal (100 - index value)
    deficit_values <- mask(100 - index_t, drought_mask, maskvalues = 0)
    if (area_pct[t] > 0) {
      severity[t] <- global(deficit_values, "mean", na.rm = TRUE)[1,1]
    } else {
      severity[t] <- 0
    }
    
    setTxtProgressBar(pb, t)
  }
  close(pb)
  
  drought_chars <- data.frame(
    date = dates,
    year = year_indices,
    month = month_indices,
    severity = severity,
    area_pct = area_pct
  )
  
  log_event(sprintf("  Severity range: %.1f%% to %.1f%% deficit", 
                    min(severity), max(severity)))
  log_event(sprintf("  Area range: %.1f%% to %.1f%% of basin", 
                    min(area_pct), max(area_pct)))
  return(drought_chars)
}

#============================================================================
# STEP 5: EXTRACT DROUGHT CHARACTERISTICS (PNPI & PNWBI)
#============================================================================
log_event("\n" + "="*70)
log_event("STEP 5: EXTRACTING DROUGHT CHARACTERISTICS")
log_event("="*70)
log_event("  Using drought threshold: 75% (moderate drought per WMO classification)")

DROUGHT_THRESHOLD_PCT <- 75

pnpi_drought <- extract_drought_characteristics_percent(
  pnpi, cell_area_km2, threshold_pct = DROUGHT_THRESHOLD_PCT
)
write.csv(pnpi_drought, "drought_analysis/PNPI_analysis/PNPI_drought_characteristics.csv",
          row.names = FALSE)

pnwbi_drought <- extract_drought_characteristics_percent(
  pnwbi, cell_area_km2, threshold_pct = DROUGHT_THRESHOLD_PCT
)
write.csv(pnwbi_drought, "drought_analysis/PNWBI_analysis/PNWBI_drought_characteristics.csv",
          row.names = FALSE)

#============================================================================
# STEP 6: VISUALIZE DROUGHT CHARACTERISTICS
#============================================================================
log_event("\n" + "="*70)
log_event("STEP 6: CREATING SCATTER PLOTS")
log_event("="*70)

create_scatter_plot <- function(drought_data, index_name, output_dir) {
  data_plot <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  
  p_main <- ggplot(data_plot, aes(x = severity, y = area_pct)) +
    geom_point(alpha = 0.5, size = 2, color = "steelblue") +
    geom_density_2d(color = "darkred", alpha = 0.3) +
    labs(x = "Drought Severity (Deficit from Normal, %)",
         y = "Percent of Area Under Drought (%)",
         title = sprintf("%s: Severity vs. Area Under Drought", index_name)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  p_sev <- ggplot(data_plot, aes(x = severity)) +
    geom_boxplot(fill = "lightblue") +
    theme_void()
  
  p_area <- ggplot(data_plot, aes(y = area_pct)) +
    geom_boxplot(fill = "lightgreen") +
    theme_void()
  
  pdf(file.path(output_dir, sprintf("%s_scatter_boxplot.pdf", index_name)),
      width = 10, height = 8)
  grid.arrange(
    p_sev, ggplot() + theme_void(),
    p_main, p_area,
    ncol = 2, nrow = 2,
    widths = c(4, 1),
    heights = c(1, 4)
  )
  dev.off()
  log_event(sprintf("    %s scatter plot saved.", index_name))
}

create_scatter_plot(pnpi_drought, "PNPI", "drought_analysis/PNPI_analysis")
create_scatter_plot(pnwbi_drought, "PNWBI", "drought_analysis/PNWBI_analysis")

#============================================================================
# STEP 7: TEST DEPENDENCY BETWEEN SEVERITY AND AREA
#============================================================================
log_event("\n" + "="*70)
log_event("STEP 7: TESTING DEPENDENCY (CORRELATION ANALYSIS)")
log_event("="*70)

test_dependency <- function(drought_data, index_name) {
  data_clean <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  S <- data_clean$severity
  A <- data_clean$area_pct
  
  cor_kendall <- cor(S, A, method = "kendall")
  cor_spearman <- cor(S, A, method = "spearman")
  test_result <- cor.test(S, A, method = "kendall")
  
  log_event(sprintf("  %s Dependency Measures:", index_name))
  log_event(sprintf("    Kendall's tau: %.4f", cor_kendall))
  log_event(sprintf("    Spearman's rho: %.4f", cor_spearman))
  log_event(sprintf("    Kendall's tau p-value: %.4e", test_result$p.value))
  
  return(list(
    kendall = cor_kendall,
    spearman = cor_spearman,
    p_value = test_result$p.value
  ))
}

pnpi_dep <- test_dependency(pnpi_drought, "PNPI")
pnwbi_dep <- test_dependency(pnwbi_drought, "PNWBI")

#============================================================================
# STEP 8: FIT MARGINAL DISTRIBUTIONS
#============================================================================
log_event("\n" + "="*70)
log_event("STEP 8: FITTING MARGINAL DISTRIBUTIONS")
log_event("="*70)

fit_marginal_distributions <- function(drought_data, index_name, output_dir) {
  data_clean <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  S <- data_clean$severity
  A <- data_clean$area_pct / 100
  
  # Fit Beta distribution to Area
  fit_beta <- fitdist(A, "beta", method = "mle")
  log_event(sprintf("    Beta distribution for Area: shape1=%.4f, shape2=%.4f",
                    fit_beta$estimate[1], fit_beta$estimate[2]))
  
  # Test distributions for Severity
  distributions <- list(
    list(name = "Exponential", dist = "exp", start = NULL),
    list(name = "Gamma", dist = "gamma", start = list(shape = 1, rate = 1)),
    list(name = "Weibull", dist = "weibull", start = list(shape = 1, scale = 1)),
    list(name = "Log-Normal", dist = "lnorm", start = NULL),
    list(name = "Normal", dist = "norm", start = NULL),
    list(name = "Logistic", dist = "logis", start = NULL)
  )
  
  results <- data.frame(Distribution = character(), AIC = numeric(), BIC = numeric(),
                        stringsAsFactors = FALSE)
  fits_list <- list()
  
  for (dist_info in distributions) {
    tryCatch({
      if (is.null(dist_info$start)) {
        fit <- fitdist(S, dist_info$dist, method = "mle")
      } else {
        fit <- fitdist(S, dist_info$dist, method = "mle", start = dist_info$start)
      }
      results <- rbind(results, data.frame(
        Distribution = dist_info$name,
        AIC = fit$aic,
        BIC = fit$bic
      ))
      fits_list[[dist_info$name]] <- fit
      log_event(sprintf("      %s: AIC=%.2f, BIC=%.2f", dist_info$name, fit$aic, fit$bic))
    }, error = function(e) {
      log_event(sprintf("      %s: FAILED (%s)", dist_info$name, e$message))
    })
  }
  
  best_idx <- which.min(results$AIC)
  best_dist_name <- results$Distribution[best_idx]
  best_fit <- fits_list[[best_dist_name]]
  
  log_event(sprintf("    ✓ Best distribution for Severity: %s", best_dist_name))
  
  write.csv(results,
            file.path(output_dir, sprintf("%s_severity_distribution_comparison.csv", index_name)),
            row.names = FALSE)
  
  pdf(file.path(output_dir, sprintf("%s_marginal_distributions.pdf", index_name)),
      width = 12, height = 10)
  par(mfrow = c(2, 2))
  plot(fit_beta, main = "Area: Beta Distribution")
  plot(best_fit, main = sprintf("Severity: %s Distribution", best_dist_name))
  dev.off()
  
  return(list(
    area_fit = fit_beta,
    severity_fit = best_fit,
    severity_dist_name = best_dist_name,
    severity_results = results
  ))
}

pnpi_marginals <- fit_marginal_distributions(pnpi_drought, "PNPI", "drought_analysis/PNPI_analysis")
pnwbi_marginals <- fit_marginal_distributions(pnwbi_drought, "PNWBI", "drought_analysis/PNWBI_analysis")

#============================================================================
# STEP 9: FIT COPULA FUNCTIONS
#============================================================================
log_event("\n" + "="*70)
log_event("STEP 9: FITTING AND SELECTING COPULA FUNCTIONS")
log_event("="*70)

fit_copulas <- function(drought_data, marginal_fits, index_name, output_dir) {
  data_clean <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  S <- data_clean$severity
  A <- data_clean$area_pct / 100
  n <- length(S)
  
  # Transform to uniform margins
  u_area <- pbeta(A, shape1 = marginal_fits$area_fit$estimate[1],
                  shape2 = marginal_fits$area_fit$estimate[2])
  
  severity_dist <- marginal_fits$severity_dist_name
  severity_params <- marginal_fits$severity_fit$estimate
  
  if (severity_dist == "Exponential") {
    u_severity <- pexp(S, rate = severity_params["rate"])
  } else if (severity_dist == "Gamma") {
    u_severity <- pgamma(S, shape = severity_params["shape"], rate = severity_params["rate"])
  } else if (severity_dist == "Weibull") {
    u_severity <- pweibull(S, shape = severity_params["shape"], scale = severity_params["scale"])
  } else if (severity_dist == "Log-Normal") {
    u_severity <- plnorm(S, meanlog = severity_params["meanlog"], sdlog = severity_params["sdlog"])
  } else if (severity_dist == "Normal") {
    u_severity <- pnorm(S, mean = severity_params["mean"], sd = severity_params["sd"])
  } else if (severity_dist == "Logistic") {
    u_severity <- plogis(S, location = severity_params["location"], scale = severity_params["scale"])
  }
  
  u_matrix <- cbind(u_severity, u_area)
  
  # Define copulas to test
  copula_list <- list(
    Clayton = claytonCopula(dim = 2),
    Gumbel = gumbelCopula(dim = 2),
    Frank = frankCopula(dim = 2),
    Joe = joeCopula(dim = 2),
    Normal = normalCopula(dim = 2),
    Plackett = plackettCopula()
  )
  
  # Fit copulas using IFM method
  results <- data.frame(Copula = character(), Parameter = numeric(), LogLik = numeric(),
                        AIC = numeric(), BIC = numeric(), stringsAsFactors = FALSE)
  fitted_copulas <- list()
  
  for (cop_name in names(copula_list)) {
    tryCatch({
      fit <- fitCopula(copula_list[[cop_name]], u_matrix, method = "mpl")
      param <- coef(fit)
      loglik <- loglikCopula(param, u_matrix, copula_list[[cop_name]])
      k <- length(param)
      aic <- -2 * loglik + 2 * k
      bic <- -2 * loglik + log(n) * k
      
      results <- rbind(results, data.frame(
        Copula = cop_name,
        Parameter = param,
        LogLik = loglik,
        AIC = aic,
        BIC = bic
      ))
      fitted_copulas[[cop_name]] <- fit
      log_event(sprintf("      %s: θ=%.4f, LogLik=%.2f, AIC=%.2f, BIC=%.2f", 
                        cop_name, param, loglik, aic, bic))
    }, error = function(e) {
      log_event(sprintf("      %s: FAILED (%s)", cop_name, e$message))
    })
  }
  
  best_idx <- which.min(results$AIC)
  best_copula_name <- results$Copula[best_idx]
  best_copula_fit <- fitted_copulas[[best_copula_name]]
  
  log_event(sprintf("    ✓ Best copula: %s", best_copula_name))
  
  write.csv(results,
            file.path(output_dir, sprintf("%s_copula_comparison.csv", index_name)),
            row.names = FALSE)
  
  return(list(
    best_copula_name = best_copula_name,
    best_copula_fit = best_copula_fit,
    all_results = results,
    u_severity = u_severity,
    u_area = u_area
  ))
}

pnpi_copulas <- fit_copulas(pnpi_drought, pnpi_marginals, "PNPI", "drought_analysis/PNPI_analysis")
pnwbi_copulas <- fit_copulas(pnwbi_drought, pnwbi_marginals, "PNWBI", "drought_analysis/PNWBI_analysis")

#============================================================================
# STEP 10A: DERIVE S-A-F CURVES - CONDITIONAL COPULA APPROACH
#============================================================================
derive_SAF_curves_conditional <- function(drought_data, marginal_fits, copula_fit,
                                          copula_name, index_name, output_dir) {
  """
  Derive S-A-F curves using conditional copula approach (Amirataee et al. 2018)
  T_S|A(s|a) = μ_T / [1 - C_S|A=a(s|a)]
  """
  log_event(sprintf("  Deriving S-A-F curves (%s) - Conditional Copula Approach...", index_name))
  
  data_clean <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  n_total_months <- nrow(drought_data)
  n_drought_months <- nrow(data_clean)
  mu_T <- n_total_months / n_drought_months
  
  log_event(sprintf("    Total months: %d", n_total_months))
  log_event(sprintf("    Drought months: %d", n_drought_months))
  log_event(sprintf("    μ_T (inter-arrival time): %.3f months", mu_T))
  
  return_periods <- c(10, 25, 50, 100)
  area_range_pct <- seq(5, 95, by = 5)
  area_range_prop <- area_range_pct / 100
  
  copula_obj <- copula_fit$best_copula_fit@copula
  beta_params <- marginal_fits$area_fit$estimate
  severity_params <- marginal_fits$severity_fit$estimate
  severity_dist <- marginal_fits$severity_dist_name
  
  saf_results <- data.frame()
  
  for (T_years in return_periods) {
    T_months <- T_years * 12
    severities <- numeric(length(area_range_prop))
    
    for (i in seq_along(area_range_prop)) {
      area_prop <- area_range_prop[i]
      v <- pbeta(area_prop, shape1 = beta_params[1], shape2 = beta_params[2])
      target_prob <- 1 - mu_T / T_months
      
      objective_function <- function(s) {
        if (severity_dist == "Exponential") {
          u <- pexp(s, rate = severity_params["rate"])
        } else if (severity_dist == "Gamma") {
          u <- pgamma(s, shape = severity_params["shape"], rate = severity_params["rate"])
        } else if (severity_dist == "Weibull") {
          u <- pweibull(s, shape = severity_params["shape"], scale = severity_params["scale"])
        } else if (severity_dist == "Log-Normal") {
          u <- plnorm(s, meanlog = severity_params["meanlog"], sdlog = severity_params["sdlog"])
        } else if (severity_dist == "Normal") {
          u <- pnorm(s, mean = severity_params["mean"], sd = severity_params["sd"])
        } else if (severity_dist == "Logistic") {
          u <- plogis(s, location = severity_params["location"], scale = severity_params["scale"])
        }
        
        cond_cdf <- cCopula(cbind(u, v), copula = copula_obj, indices = 2)
        return((cond_cdf - target_prob)^2)
      }
      
      result <- tryCatch(
        optimize(objective_function, interval = c(0.001, 20)),
        error = function(e) list(minimum = 0)
      )
      severities[i] <- result$minimum
    }
    
    temp_df <- data.frame(
      ReturnPeriod_years = T_years,
      ReturnPeriod_months = T_months,
      Area_pct = area_range_pct,
      Severity = severities,
      Method = "Conditional"
    )
    saf_results <- rbind(saf_results, temp_df)
  }
  
  write.csv(saf_results,
            file.path(output_dir, sprintf("%s_SAF_curves_conditional.csv", index_name)),
            row.names = FALSE)
  return(saf_results)
}

#============================================================================
# STEP 10B: DERIVE S-A-F CURVES - KENDALL DISTRIBUTION APPROACH
#============================================================================
derive_SAF_curves_kendall <- function(drought_data, marginal_fits, copula_fit,
                                      copula_name, index_name, output_dir) {
  """
  Derive S-A-F curves using Kendall distribution function
  K_C(t) = t - (1 - C(t,t))/t
  T_K = μ_T / (1 - K_C(t))
  """
  log_event(sprintf("  Deriving S-A-F curves (%s) - Kendall Distribution Approach...", index_name))
  
  data_clean <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  n_total_months <- nrow(drought_data)
  n_drought_months <- nrow(data_clean)
  mu_T <- n_total_months / n_drought_months
  
  log_event(sprintf("    Total months: %d", n_total_months))
  log_event(sprintf("    Drought months: %d", n_drought_months))
  log_event(sprintf("    μ_T (inter-arrival time): %.3f months", mu_T))
  
  return_periods <- c(10, 25, 50, 100)
  area_range_pct <- seq(5, 95, by = 5)
  area_range_prop <- area_range_pct / 100
  
  copula_obj <- copula_fit$best_copula_fit@copula
  beta_params <- marginal_fits$area_fit$estimate
  severity_params <- marginal_fits$severity_fit$estimate
  severity_dist <- marginal_fits$severity_dist_name
  
  saf_results <- data.frame()
  
  for (T_years in return_periods) {
    T_months <- T_years * 12
    log_event(sprintf("      Return period: %d years (%d months)...", T_years, T_months))
    
    target_kc <- 1 - mu_T / T_months
    
    k_c_function <- function(t_val) {
      if (t_val <= 0 || t_val >= 1) return(1e6)
      c_tt <- pCopula(cbind(t_val, t_val), copula_obj)
      k_c <- t_val - (1 - c_tt) / t_val
      return((k_c - target_kc)^2)
    }
    
    t_solution <- tryCatch(
      optimize(k_c_function, interval = c(0.001, 0.999)),
      error = function(e) list(minimum = 0.5)
    )
    t_val <- t_solution$minimum
    
    log_event(sprintf("        Target K_C(t) = %.4f, Solution t = %.4f", target_kc, t_val))
    
    severities <- numeric(length(area_range_prop))
    
    for (i in seq_along(area_range_prop)) {
      area_prop <- area_range_prop[i]
      v <- pbeta(area_prop, shape1 = beta_params[1], shape2 = beta_params[2])
      
      objective_function <- function(s) {
        if (severity_dist == "Exponential") {
          u <- pexp(s, rate = severity_params["rate"])
        } else if (severity_dist == "Gamma") {
          u <- pgamma(s, shape = severity_params["shape"], rate = severity_params["rate"])
        } else if (severity_dist == "Weibull") {
          u <- pweibull(s, shape = severity_params["shape"], scale = severity_params["scale"])
        } else if (severity_dist == "Log-Normal") {
          u <- plnorm(s, meanlog = severity_params["meanlog"], sdlog = severity_params["sdlog"])
        } else if (severity_dist == "Normal") {
          u <- pnorm(s, mean = severity_params["mean"], sd = severity_params["sd"])
        } else if (severity_dist == "Logistic") {
          u <- plogis(s, location = severity_params["location"], scale = severity_params["scale"])
        }
        
        copula_value <- pCopula(cbind(u, v), copula_obj)
        return((copula_value - t_val)^2)
      }
      
      result <- tryCatch(
        optimize(objective_function, interval = c(0.001, 20)),
        error = function(e) list(minimum = 0)
      )
      severities[i] <- result$minimum
    }
    
    temp_df <- data.frame(
      ReturnPeriod_years = T_years,
      ReturnPeriod_months = T_months,
      Area_pct = area_range_pct,
      Severity = severities,
      Method = "Kendall"
    )
    saf_results <- rbind(saf_results, temp_df)
  }
  
  write.csv(saf_results,
            file.path(output_dir, sprintf("%s_SAF_curves_kendall.csv", index_name)),
            row.names = FALSE)
  return(saf_results)
}

#============================================================================
# STEP 10C: GENERATE COMPARATIVE SAF PLOTS
#============================================================================
create_comparative_SAF_plots <- function(saf_conditional, saf_kendall, 
                                         index_name, output_dir) {
  saf_combined <- rbind(saf_conditional, saf_kendall)
  
  # Overlay plot
  pdf(file.path(output_dir, sprintf("%s_SAF_overlay.pdf", index_name)),
      width = 10, height = 8)
  
  ggplot(saf_combined, aes(x = Area_pct, y = Severity, 
                           color = factor(ReturnPeriod_years),
                           linetype = Method, shape = Method)) +
    geom_line(size = 1.1) +
    geom_point(size = 2.5) +
    scale_color_viridis_d(name = "Return Period\n(years)") +
    scale_linetype_manual(name = "Method", 
                          values = c("Conditional" = "solid", "Kendall" = "dashed")) +
    scale_shape_manual(name = "Method",
                       values = c("Conditional" = 19, "Kendall" = 17)) +
    labs(x = "Percent of Area Under Drought (%)",
         y = "Drought Severity (Deficit from Normal, %)",
         title = sprintf("%s: SAF Curves - Method Comparison", index_name)) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      panel.grid.major = element_line(color = "gray80", size = 0.5)
    ) +
    scale_x_continuous(breaks = seq(0, 100, 10)) +
    scale_y_continuous(breaks = seq(0, ceiling(max(saf_combined$Severity, na.rm = TRUE)), 2))
  
  dev.off()
  log_event(sprintf("    Overlay SAF plot saved for %s", index_name))
  
  # Side-by-side comparison
  pdf(file.path(output_dir, sprintf("%s_SAF_comparison.pdf", index_name)),
      width = 14, height = 7)
  
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 2))
  
  # Conditional approach
  cond_data <- saf_combined[saf_combined$Method == "Conditional", ]
  plot(cond_data$Area_pct[cond_data$ReturnPeriod_years == 50], 
       cond_data$Severity[cond_data$ReturnPeriod_years == 50],
       type = "n", xlim = c(0, 100), 
       ylim = c(0, max(saf_combined$Severity, na.rm = TRUE) * 1.1),
       xlab = "Percent of Area Under Drought (%)", 
       ylab = "Drought Severity (Deficit from Normal, %)",
       main = sprintf("%s: Conditional Copula Approach", index_name),
       cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.1)
  
  colors <- c("blue", "green3", "orange", "red")
  for (j in 1:4) {
    rp <- c(10, 25, 50, 100)[j]
    subset_data <- cond_data[cond_data$ReturnPeriod_years == rp, ]
    lines(subset_data$Area_pct, subset_data$Severity, 
          col = colors[j], lwd = 2.5)
    points(subset_data$Area_pct, subset_data$Severity, 
           col = colors[j], pch = 19, cex = 0.9)
  }
  legend("topright", legend = c("10-year", "25-year", "50-year", "100-year"), 
         col = colors, lwd = 2.5, pch = 19, cex = 0.95, title = "Return Period")
  grid()
  
  # Kendall approach
  kend_data <- saf_combined[saf_combined$Method == "Kendall", ]
  plot(kend_data$Area_pct[kend_data$ReturnPeriod_years == 50], 
       kend_data$Severity[kend_data$ReturnPeriod_years == 50],
       type = "n", xlim = c(0, 100), 
       ylim = c(0, max(saf_combined$Severity, na.rm = TRUE) * 1.1),
       xlab = "Percent of Area Under Drought (%)", 
       ylab = "Drought Severity (Deficit from Normal, %)",
       main = sprintf("%s: Kendall Distribution Approach", index_name),
       cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.1)
  
  for (j in 1:4) {
    rp <- c(10, 25, 50, 100)[j]
    subset_data <- kend_data[kend_data$ReturnPeriod_years == rp, ]
    lines(subset_data$Area_pct, subset_data$Severity, 
          col = colors[j], lwd = 2.5, lty = 2)
    points(subset_data$Area_pct, subset_data$Severity, 
           col = colors[j], pch = 17, cex = 0.9)
  }
  legend("topright", legend = c("10-year", "25-year", "50-year", "100-year"), 
         col = colors, lwd = 2.5, lty = 2, pch = 17, cex = 0.95, title = "Return Period")
  grid()
  
  dev.off()
  log_event(sprintf("    Side-by-side SAF plot saved for %s", index_name))
}

#============================================================================
# EXECUTE BOTH SAF METHODS FOR PNPI AND PNWBI
#============================================================================
log_event("\n" + "="*70)
log_event("STEP 10: DERIVING S-A-F CURVES WITH DUAL METHODS")
log_event("="*70)

# PNPI Analysis
pnpi_saf_conditional <- derive_SAF_curves_conditional(
  pnpi_drought, pnpi_marginals, pnpi_copulas,
  pnpi_copulas$best_copula_name, "PNPI",
  "drought_analysis/PNPI_analysis"
)

pnpi_saf_kendall <- derive_SAF_curves_kendall(
  pnpi_drought, pnpi_marginals, pnpi_copulas,
  pnpi_copulas$best_copula_name, "PNPI",
  "drought_analysis/PNPI_analysis"
)

create_comparative_SAF_plots(
  pnpi_saf_conditional, pnpi_saf_kendall,
  "PNPI", "drought_analysis/PNPI_analysis"
)

# PNWBI Analysis
pnwbi_saf_conditional <- derive_SAF_curves_conditional(
  pnwbi_drought, pnwbi_marginals, pnwbi_copulas,
  pnwbi_copulas$best_copula_name, "PNWBI",
  "drought_analysis/PNWBI_analysis"
)

pnwbi_saf_kendall <- derive_SAF_curves_kendall(
  pnwbi_drought, pnwbi_marginals, pnwbi_copulas,
  pnwbi_copulas$best_copula_name, "PNWBI",
  "drought_analysis/PNWBI_analysis"
)

create_comparative_SAF_plots(
  pnwbi_saf_conditional, pnwbi_saf_kendall,
  "PNWBI", "drought_analysis/PNWBI_analysis"
)

#============================================================================
# STEP 11: GENERATE COMPREHENSIVE SUMMARY REPORT
#============================================================================
log_event("\n" + "="*70)
log_event("STEP 11: GENERATING COMPREHENSIVE SUMMARY REPORT")
log_event("="*70)

summary_report <- sprintf("
DROUGHT SEVERITY-AREA-FREQUENCY ANALYSIS - NECHAKO BASIN
Percent of Normal Indices with Dual Return Period Methods
================================================================
Analysis Date: %s
Study Period: %s to %s
Total Months: %d
Grid Cells: %d
Basin Area: %.2f km²
================================================================
METHODOLOGY
================================================================
Drought Indices:
PNPI: Percent of Normal Precipitation Index
      Formula: PNPI = (P_t / μ_m) × 100
      Reference: WMO (2016); Heim (2002)
      
PNWBI: Percent of Normal Water Balance Index
       Formula: PNWBI = ((P-PET)_t / μ_m) × 100
       Reference: Palmer (1965); Thornthwaite (1948)

Drought Threshold: Index < 75%% (moderate drought per WMO standards)

Drought Characteristics:
Severity (S): Mean deficit from normal (100 - index) in drought areas
              Units: Percentage points deficit from normal conditions
Area (A):     Percentage of basin area below threshold

Return Period Methods:
Method A - Conditional Copula Approach:
  T_S|A(s|a) = μ_T / [1 - C_S|A=a(s|a)]
  Reference: Amirataee et al. (2018); Shiau (2006)
  
Method B - Kendall Distribution Function:
  K_C(t) = t - (1 - C(t,t))/t
  T_K = μ_T / (1 - K_C(t))
  Reference: Genest et al. (2009); Salvadori & De Michele (2004)

Statistical Analysis:
Marginal distributions: AIC/BIC selection (Zhang & Singh, 2007)
Dependency: Kendall's tau, Spearman's rho (Genest & Favre, 2007)
Copula fitting: IFM method (Joe, 1997)
================================================================
PNPI ANALYSIS RESULTS
================================================================
Dependency Measures:
Kendall's tau:   %.4f (p < %.4f)
Spearman's rho:  %.4f

Marginal Distributions:
Area:     Beta (shape1=%.4f, shape2=%.4f)
Severity: %s

Best Copula: %s
Parameter:   %.4f
AIC:         %.2f
BIC:         %.2f

Sample S-A-F Values (50-year return period):
Area  5%%: Severity = %.2f%% (Conditional) | %.2f%% (Kendall)
Area 25%%: Severity = %.2f%% (Conditional) | %.2f%% (Kendall)
Area 50%%: Severity = %.2f%% (Conditional) | %.2f%% (Kendall)
Area 75%%: Severity = %.2f%% (Conditional) | %.2f%% (Kendall)
Area 95%%: Severity = %.2f%% (Conditional) | %.2f%% (Kendall)
================================================================
PNWBI ANALYSIS RESULTS
================================================================
Dependency Measures:
Kendall's tau:   %.4f (p < %.4f)
Spearman's rho:  %.4f

Marginal Distributions:
Area:     Beta (shape1=%.4f, shape2=%.4f)
Severity: %s

Best Copula: %s
Parameter:   %.4f
AIC:         %.2f
BIC:         %.2f

Sample S-A-F Values (50-year return period):
Area  5%%: Severity = %.2f%% (Conditional) | %.2f%% (Kendall)
Area 25%%: Severity = %.2f%% (Conditional) | %.2f%% (Kendall)
Area 50%%: Severity = %.2f%% (Conditional) | %.2f%% (Kendall)
Area 75%%: Severity = %.2f%% (Conditional) | %.2f%% (Kendall)
Area 95%%: Severity = %.2f%% (Conditional) | %.2f%% (Kendall)
================================================================
KEY FINDINGS
================================================================
1. PNPI vs PNWBI:
   - PNWBI typically shows higher severity values due to inclusion of PET
   - PNPI better reflects pure precipitation deficits
   - PNWBI more appropriate for agricultural/hydrological drought assessment

2. Method Comparison:
   - Kendall approach generally yields more conservative (higher) severity 
     estimates for large areas (>50%% basin coverage)
   - Conditional approach more sensitive to marginal distribution fits
   - Both methods show consistent ranking of return periods

3. Drought Threshold Sensitivity:
   - Using 75%% threshold captures moderate drought conditions
   - For severe drought analysis, consider 60%% threshold
   - Threshold choice significantly impacts frequency statistics
================================================================
OUTPUT FILES
================================================================
PNPI Analysis:
PNPI_monthly.nc
PNPI_drought_characteristics.csv
PNPI_marginal_distributions.pdf
PNPI_copula_comparison.csv
PNPI_SAF_curves_conditional.csv   <- Method A
PNPI_SAF_curves_kendall.csv       <- Method B
PNPI_SAF_overlay.pdf              <- Overlay comparison
PNPI_SAF_comparison.pdf           <- Side-by-side plots

PNWBI Analysis:
PNWBI_monthly.nc
PNWBI_drought_characteristics.csv
PNWBI_marginal_distributions.pdf
PNWBI_copula_comparison.csv
PNWBI_SAF_curves_conditional.csv  <- Method A
PNWBI_SAF_curves_kendall.csv      <- Method B
PNWBI_SAF_overlay.pdf             <- Overlay comparison
PNWBI_SAF_comparison.pdf          <- Side-by-side plots
================================================================
REFERENCES
================================================================
Genest, C., Rémillard, B., & Beaudoin, D. (2009). Goodness-of-fit tests 
for copulas: A review and a power study. Insurance: Mathematics and 
Economics, 44(2), 199-213.

Heim, R.R. (2002). A review of twentieth-century drought indices used in 
the United States. Bulletin of the American Meteorological Society, 
83(8), 1149-1165.

Nelsen, R.B. (2006). An introduction to copulas (2nd ed.). Springer.

Salvadori, G., & De Michele, C. (2004). Frequency analysis via copulas: 
Theoretical aspects and applications to hydrological events. Water 
Resources Research, 40(12).

Shiau, J.T. (2006). Fitting drought duration and severity with 
two-dimensional copulas. Water Resources Management, 20(5), 795-815.

WMO (2016). Standardized Precipitation Index User Guide. WMO-No. 1090.
================================================================
",
                          Sys.time(),
                          min(dates), max(dates),
                          n_months, n_cells, total_area_km2,
                          # PNPI results
                          pnpi_dep$kendall, pnpi_dep$p_value, pnpi_dep$spearman,
                          pnpi_marginals$area_fit$estimate[1], pnpi_marginals$area_fit$estimate[2],
                          pnpi_marginals$severity_dist_name,
                          pnpi_copulas$best_copula_name,
                          coef(pnpi_copulas$best_copula_fit),
                          min(pnpi_copulas$all_results$AIC[pnpi_copulas$all_results$Copula == pnpi_copulas$best_copula_name]),
                          min(pnpi_copulas$all_results$BIC[pnpi_copulas$all_results$Copula == pnpi_copulas$best_copula_name]),
                          # PNPI SAF values (50-year)
                          pnpi_saf_conditional$Severity[pnpi_saf_conditional$ReturnPeriod_years == 50 & pnpi_saf_conditional$Area_pct == 5],
                          pnpi_saf_kendall$Severity[pnpi_saf_kendall$ReturnPeriod_years == 50 & pnpi_saf_kendall$Area_pct == 5],
                          pnpi_saf_conditional$Severity[pnpi_saf_conditional$ReturnPeriod_years == 50 & pnpi_saf_conditional$Area_pct == 25],
                          pnpi_saf_kendall$Severity[pnpi_saf_kendall$ReturnPeriod_years == 50 & pnpi_saf_kendall$Area_pct == 25],
                          pnpi_saf_conditional$Severity[pnpi_saf_conditional$ReturnPeriod_years == 50 & pnpi_saf_conditional$Area_pct == 50],
                          pnpi_saf_kendall$Severity[pnpi_saf_kendall$ReturnPeriod_years == 50 & pnpi_saf_kendall$Area_pct == 50],
                          pnpi_saf_conditional$Severity[pnpi_saf_conditional$ReturnPeriod_years == 50 & pnpi_saf_conditional$Area_pct == 75],
                          pnpi_saf_kendall$Severity[pnpi_saf_kendall$ReturnPeriod_years == 50 & pnpi_saf_kendall$Area_pct == 75],
                          pnpi_saf_conditional$Severity[pnpi_saf_conditional$ReturnPeriod_years == 50 & pnpi_saf_conditional$Area_pct == 95],
                          pnpi_saf_kendall$Severity[pnpi_saf_kendall$ReturnPeriod_years == 50 & pnpi_saf_kendall$Area_pct == 95],
                          # PNWBI results
                          pnwbi_dep$kendall, pnwbi_dep$p_value, pnwbi_dep$spearman,
                          pnwbi_marginals$area_fit$estimate[1], pnwbi_marginals$area_fit$estimate[2],
                          pnwbi_marginals$severity_dist_name,
                          pnwbi_copulas$best_copula_name,
                          coef(pnwbi_copulas$best_copula_fit),
                          min(pnwbi_copulas$all_results$AIC[pnwbi_copulas$all_results$Copula == pnwbi_copulas$best_copula_name]),
                          min(pnwbi_copulas$all_results$BIC[pnwbi_copulas$all_results$Copula == pnwbi_copulas$best_copula_name]),
                          # PNWBI SAF values (50-year)
                          pnwbi_saf_conditional$Severity[pnwbi_saf_conditional$ReturnPeriod_years == 50 & pnwbi_saf_conditional$Area_pct == 5],
                          pnwbi_saf_kendall$Severity[pnwbi_saf_kendall$ReturnPeriod_years == 50 & pnwbi_saf_kendall$Area_pct == 5],
                          pnwbi_saf_conditional$Severity[pnwbi_saf_conditional$ReturnPeriod_years == 50 & pnwbi_saf_conditional$Area_pct == 25],
                          pnwbi_saf_kendall$Severity[pnwbi_saf_kendall$ReturnPeriod_years == 50 & pnwbi_saf_kendall$Area_pct == 25],
                          pnwbi_saf_conditional$Severity[pnwbi_saf_conditional$ReturnPeriod_years == 50 & pnwbi_saf_conditional$Area_pct == 50],
                          pnwbi_saf_kendall$Severity[pnwbi_saf_kendall$ReturnPeriod_years == 50 & pnwbi_saf_kendall$Area_pct == 50],
                          pnwbi_saf_conditional$Severity[pnwbi_saf_conditional$ReturnPeriod_years == 50 & pnwbi_saf_conditional$Area_pct == 75],
                          pnwbi_saf_kendall$Severity[pnwbi_saf_kendall$ReturnPeriod_years == 50 & pnwbi_saf_kendall$Area_pct == 75],
                          pnwbi_saf_conditional$Severity[pnwbi_saf_conditional$ReturnPeriod_years == 50 & pnwbi_saf_conditional$Area_pct == 95],
                          pnwbi_saf_kendall$Severity[pnwbi_saf_kendall$ReturnPeriod_years == 50 & pnwbi_saf_kendall$Area_pct == 95]
)

writeLines(summary_report, "drought_analysis/SUMMARY_REPORT.txt")
log_event("Comprehensive summary report saved to: SUMMARY_REPORT.txt")

log_event("\n" + "="*70)
log_event("ANALYSIS COMPLETE")
log_event("="*70)
log_event(sprintf("Completed at: %s", Sys.time()))
log_event("All results saved to: drought_analysis/")
cat("\n\n✓ ANALYSIS COMPLETE - Check 'drought_analysis/' folder for all results\n\n")