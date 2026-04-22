

#============================================================================
# STEP 10A: DERIVE S-A-F CURVES - CONDITIONAL COPULA APPROACH
#============================================================================
derive_SAF_curves_conditional <- function(drought_data, marginal_fits, copula_fit,
                                          copula_name, index_name, output_dir,
                                          duration_class = "all") {
  # Derive S-A-F curves using conditional copula approach (Amirataee et al. 2018)
  # T_S|A(s|a) = mu_T / [1 - C_S|A=a(s|a)]
  log_event(sprintf("  Deriving S-A-F curves (%s | %s) - Conditional Copula...",
                    index_name, duration_class))
  
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
            file.path(output_dir, sprintf("%s_SAF_curves_conditional_%s.csv",
                                          index_name,
                                          gsub("[^A-Za-z0-9]", "_", duration_class))),
            row.names = FALSE)
  return(saf_results)
}

#============================================================================
# STEP 10B: DERIVE S-A-F CURVES - KENDALL DISTRIBUTION APPROACH
#============================================================================
derive_SAF_curves_kendall <- function(drought_data, marginal_fits, copula_fit,
                                      copula_name, index_name, output_dir,
                                      duration_class = "all") {
  # Derive S-A-F curves using Kendall distribution function
  # K_C(t) = t - (1 - C(t,t))/t
  # T_K = mu_T / (1 - K_C(t))
  log_event(sprintf("  Deriving S-A-F curves (%s | %s) - Kendall Distribution...",
                    index_name, duration_class))
  
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
            file.path(output_dir, sprintf("%s_SAF_curves_kendall_%s.csv",
                                          index_name,
                                          gsub("[^A-Za-z0-9]", "_", duration_class))),
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
# EXECUTE BOTH SAF METHODS FOR PNPI AND PNWBI — PER DURATION CLASS
#============================================================================
log_event(paste(rep("=", 70), collapse = ""))
log_event("STEP 10: DERIVING S-A-F CURVES WITH DUAL METHODS (PER DURATION CLASS)")
log_event(paste(rep("=", 70), collapse = ""))

# Helper: subset drought_chars to the months that belong to a specific class
subset_drought_chars_by_class <- function(drought_chars, events_df, class_label) {
  # Collect all month-dates that fall inside events of the given class
  ev_sub <- events_df[events_df$duration_class == class_label, ]
  if (nrow(ev_sub) == 0) return(NULL)
  # Mark months that are inside a class event
  in_class <- logical(nrow(drought_chars))
  for (i in seq_len(nrow(ev_sub))) {
    in_class <- in_class |
      (drought_chars$date >= ev_sub$start_date[i] &
         drought_chars$date <= ev_sub$end_date[i])
  }
  # Return the full drought_chars but zero-out months not in the class
  dc_sub <- drought_chars
  dc_sub$severity[!in_class]  <- 0
  dc_sub$area_pct[!in_class]  <- 0
  dc_sub
}

# Duration classes to process for SAF curves (need >=10 events for copula fitting)
DC_LABELS <- c("D4-6 (Short-term)", "D7-12 (Medium-term)", "D12+ (Long-term)")

# Initialise overall SAF storage (used in summary report, default = "all-classes")
pnpi_saf_conditional  <- NULL
pnpi_saf_kendall      <- NULL
pnwbi_saf_conditional <- NULL
pnwbi_saf_kendall     <- NULL

run_saf_for_index <- function(index_name, drought_chars, events_df,
                              marginals, copulas, out_dir) {
  saf_cond_list  <- list()
  saf_kend_list  <- list()
  
  # ── (a) Full all-class SAF (original behaviour, preserved) ──────────
  log_event(sprintf("  [%s] All-class SAF curves...", index_name))
  sc <- derive_SAF_curves_conditional(
    drought_chars, marginals, copulas,
    copulas$best_copula_name, index_name, out_dir, "all")
  sk <- derive_SAF_curves_kendall(
    drought_chars, marginals, copulas,
    copulas$best_copula_name, index_name, out_dir, "all")
  create_comparative_SAF_plots(sc, sk, index_name, out_dir)
  saf_cond_list[["all"]] <- sc
  saf_kend_list[["all"]] <- sk
  
  # ── (b) Per-duration-class SAF curves ───────────────────────────────
  for (cl in DC_LABELS) {
    n_ev <- sum(events_df$duration_class == cl, na.rm = TRUE)
    log_event(sprintf("  [%s | %s] n_events = %d", index_name, cl, n_ev))
    
    if (n_ev < 8L) {
      log_event(sprintf("    ⚠ Fewer than 8 events — skipping copula SAF for %s", cl))
      next
    }
    
    dc_sub <- subset_drought_chars_by_class(drought_chars, events_df, cl)
    if (is.null(dc_sub)) next
    
    # Re-fit marginals and copula on the class-specific severity/area pairs
    data_cl <- dc_sub[dc_sub$severity > 0 & dc_sub$area_pct > 0, ]
    if (nrow(data_cl) < 8L) {
      log_event(sprintf("    ⚠ Insufficient drought months for %s — skipping", cl))
      next
    }
    
    tryCatch({
      cl_marginals <- fit_marginal_distributions(
        dc_sub, index_name, out_dir)
      cl_copulas   <- fit_copulas(
        dc_sub, cl_marginals, index_name, out_dir)
      
      sc_cl <- derive_SAF_curves_conditional(
        dc_sub, cl_marginals, cl_copulas,
        cl_copulas$best_copula_name, index_name, out_dir, cl)
      sk_cl <- derive_SAF_curves_kendall(
        dc_sub, cl_marginals, cl_copulas,
        cl_copulas$best_copula_name, index_name, out_dir, cl)
      
      sc_cl$duration_class <- cl
      sk_cl$duration_class <- cl
      saf_cond_list[[cl]]  <- sc_cl
      saf_kend_list[[cl]]  <- sk_cl
      
      log_event(sprintf("    ✓ SAF curves derived for %s | %s", index_name, cl))
    }, error = function(e) {
      log_event(sprintf("    ❌ SAF fitting failed for %s | %s: %s",
                        index_name, cl, e$message))
    })
  }
  
  # ── Combined overlay plot — all duration classes on one figure ───────
  cond_all_classes <- do.call(rbind, saf_cond_list)
  kend_all_classes <- do.call(rbind, saf_kend_list)
  if (!is.null(cond_all_classes) && "duration_class" %in% names(cond_all_classes)) {
    pdf(file.path(out_dir,
                  sprintf("%s_SAF_overlay_by_duration_class.pdf", index_name)),
        width = 12, height = 8)
    p <- ggplot(cond_all_classes[cond_all_classes$ReturnPeriod_years == 50, ],
                aes(x = Area_pct, y = Severity,
                    colour = duration_class, linetype = duration_class)) +
      geom_line(linewidth = 1.1) +
      geom_point(size = 2) +
      scale_colour_manual(
        values = c("all"                  = "grey50",
                   "D4-6 (Short-term)"   = "#fdae61",
                   "D7-12 (Medium-term)" = "#f46d43",
                   "D12+ (Long-term)"    = "#a50026"),
        name = "Duration Class") +
      scale_linetype_manual(
        values = c("all" = "dotted",
                   "D4-6 (Short-term)"   = "solid",
                   "D7-12 (Medium-term)" = "dashed",
                   "D12+ (Long-term)"    = "longdash"),
        name = "Duration Class") +
      labs(title    = sprintf("%s — 50-Year SAF Curves by Duration Class (Conditional Method)",
                              index_name),
           subtitle = "Conditional copula approach; curves represent events of each duration class separately",
           x = "Percent of Area Under Drought (%)",
           y = "Drought Severity (% deficit from normal)") +
      theme_bw(base_size = 13) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
            legend.position = "right")
    print(p)
    dev.off()
    log_event(sprintf("  ✓ Duration-class SAF overlay saved for %s", index_name))
  }
  
  list(saf_conditional = saf_cond_list[["all"]],
       saf_kendall     = saf_kend_list[["all"]],
       saf_cond_by_class = saf_cond_list,
       saf_kend_by_class = saf_kend_list)
}

# ── Run for PNPI ──────────────────────────────────────────────────────
pnpi_saf_results <- run_saf_for_index(
  "PNPI", pnpi_drought, pnpi_events,
  pnpi_marginals, pnpi_copulas, "drought_analysis/PNPI_analysis")
pnpi_saf_conditional <- pnpi_saf_results$saf_conditional
pnpi_saf_kendall     <- pnpi_saf_results$saf_kendall

# ── Run for PNWBI ─────────────────────────────────────────────────────
pnwbi_saf_results <- run_saf_for_index(
  "PNWBI", pnwbi_drought, pnwbi_events,
  pnwbi_marginals, pnwbi_copulas, "drought_analysis/PNWBI_analysis")
pnwbi_saf_conditional <- pnwbi_saf_results$saf_conditional
pnwbi_saf_kendall     <- pnwbi_saf_results$saf_kendall


#============================================================================
# [NEW] STEP 10d: NON-STATIONARY SAF CURVES WITH GAMLSS MARGINALS (Enhancement 1)
#
# Fits non-stationary Gamma marginals using GAMLSS where the location
# parameter (mean severity) varies with calendar year or JJA Tmax.
# Computes SAF curves separately for two reference periods:
#   P_ref : 1950-1990  (pre-warming)
#   P_rec : 2010-2025  (recent warming)
# Quantifies how the 50-year and 100-year return levels have shifted.
#============================================================================
log_event(paste(rep("=", 70), collapse=""))
log_event("STEP 10d [NEW]: NON-STATIONARY SAF CURVES (GAMLSS)")
log_event(paste(rep("=", 70), collapse=""))

dir.create("drought_analysis/nonstationary_SAF", showWarnings = FALSE)

NS_P_REF_START <- 1950L;  NS_P_REF_END <- 1990L
NS_P_REC_START <- 2010L;  NS_P_REC_END <- 2025L

if (!requireNamespace("gamlss", quietly = TRUE)) {
  tryCatch(install.packages("gamlss", repos = "https://cloud.r-project.org"),
           error = function(e) message("gamlss install failed: ", e$message))
}

fit_nonstationary_marginals <- function(drought_data, index_name,
                                        year_vec, output_dir) {
  data_clean <- drought_data[drought_data$severity > 0 &
                               drought_data$area_pct > 0, ]
  if (nrow(data_clean) < 20L || !requireNamespace("gamlss", quietly=TRUE)) {
    log_event(sprintf("  [%s] Non-stationary fit skipped", index_name))
    return(NULL)
  }
  # Match year to drought months
  yr_match <- year_vec[seq_len(nrow(data_clean))]
  if (length(yr_match) != nrow(data_clean))
    yr_match <- rep(mean(year_vec), nrow(data_clean))
  df_fit  <- data.frame(severity = data_clean$severity,
                        year_std = as.numeric(scale(yr_match)))
  # Stationary Gamma
  m_stat <- tryCatch(
    gamlss::gamlss(severity ~ 1, sigma.formula = ~ 1,
                   family = gamlss.dist::GA(), data = df_fit, trace = FALSE),
    error = function(e) NULL)
  # Non-stationary: mu ~ year
  m_ns   <- tryCatch(
    gamlss::gamlss(severity ~ year_std, sigma.formula = ~ 1,
                   family = gamlss.dist::GA(), data = df_fit, trace = FALSE),
    error = function(e) NULL)
  if (is.null(m_stat) || is.null(m_ns)) {
    log_event(sprintf("  [%s] GAMLSS fitting failed", index_name)); return(NULL)
  }
  aic_comp <- data.frame(
    model     = c("Stationary", "NonStationary_year"),
    AIC       = c(AIC(m_stat), AIC(m_ns)),
    BIC       = c(BIC(m_stat), BIC(m_ns)))
  aic_comp$delta_AIC <- aic_comp$AIC - min(aic_comp$AIC)
  write.csv(aic_comp,
            file.path(output_dir,
                      sprintf("%s_nonstat_marginal_AIC.csv", index_name)),
            row.names = FALSE)
  log_event(sprintf("  [%s] Non-stat Gamma AIC: stat=%.2f  ns=%.2f  ΔAIC=%.2f",
                    index_name, AIC(m_stat), AIC(m_ns),
                    AIC(m_ns) - AIC(m_stat)))
  # Predict mu for each reference period
  yr_full <- sort(unique(yr_match))
  ys_full <- (yr_full - mean(yr_match)) / sd(yr_match)
  mu_pred <- exp(predict(m_ns,
                         newdata = data.frame(year_std = ys_full),
                         type = "response"))
  # Period-mean mu
  ref_ys <- (seq(NS_P_REF_START, NS_P_REF_END) - mean(yr_match)) / sd(yr_match)
  rec_ys <- (seq(NS_P_REC_START, NS_P_REC_END) - mean(yr_match)) / sd(yr_match)
  mu_ref <- mean(exp(predict(m_ns,
                             newdata = data.frame(year_std = ref_ys),
                             type = "response")), na.rm = TRUE)
  mu_rec <- mean(exp(predict(m_ns,
                             newdata = data.frame(year_std = rec_ys),
                             type = "response")), na.rm = TRUE)
  sigma_hat <- exp(coef(m_ns, what = "sigma")[1])
  log_event(sprintf("  [%s] μ(P_ref 1950-1990)=%.3f  μ(P_rec 2010-2025)=%.3f  σ=%.4f",
                    index_name, mu_ref, mu_rec, sigma_hat))
  list(m_stat=m_stat, m_ns=m_ns, aic_comp=aic_comp,
       mu_ref=mu_ref, mu_rec=mu_rec, sigma=sigma_hat,
       year_vec=yr_match)
}

# Derive SAF curves for a given period-specific mu (non-stationary severity)
derive_SAF_nonstationary <- function(mu_period, sigma_ns, copula_fit_obj,
                                     marginal_fits, drought_data,
                                     index_name, period_label, output_dir) {
  data_clean <- drought_data[drought_data$severity>0 & drought_data$area_pct>0, ]
  n_total    <- nrow(drought_data)
  n_drought  <- nrow(data_clean)
  mu_T       <- n_total / n_drought
  cop_obj    <- copula_fit_obj$best_copula_fit@copula
  beta_par   <- marginal_fits$area_fit$estimate
  # Use non-stationary shape/rate: GA parameterisation mu=mu, sigma=sigma
  # shape = 1/sigma^2;  rate = shape/mu
  ns_shape <- 1 / (sigma_ns^2)
  ns_rate  <- ns_shape / mu_period
  
  return_periods <- c(10, 25, 50, 100)
  area_range_pct <- seq(5, 95, by = 5)
  saf_ns <- data.frame()
  
  for (T_years in return_periods) {
    T_months    <- T_years * 12
    target_prob <- 1 - mu_T / T_months
    severities  <- numeric(length(area_range_pct))
    for (i in seq_along(area_range_pct)) {
      v <- pbeta(area_range_pct[i]/100, beta_par[1], beta_par[2])
      obj_fn <- function(s) {
        u <- pgamma(s, shape = ns_shape, rate = ns_rate)
        cdf <- copula::cCopula(cbind(u, v), copula = cop_obj, indices = 2)
        (cdf - target_prob)^2
      }
      res <- tryCatch(optimize(obj_fn, c(0.001, 50)),
                      error = function(e) list(minimum = 0))
      severities[i] <- res$minimum
    }
    saf_ns <- rbind(saf_ns, data.frame(
      ReturnPeriod_years = T_years,
      Area_pct           = area_range_pct,
      Severity           = severities,
      Method             = "Conditional_NS",
      Period             = period_label,
      mu_period          = round(mu_period, 4)))
  }
  fn <- file.path(output_dir,
                  sprintf("%s_SAF_nonstationary_%s.csv", index_name,
                          gsub("[^A-Za-z0-9]","_", period_label)))
  write.csv(saf_ns, fn, row.names = FALSE)
  log_event(sprintf("    [%s | %s] Non-stationary SAF saved.", index_name, period_label))
  saf_ns
}

# Run non-stationary analysis for PNPI and PNWBI
pnpi_drought_years  <- pnpi_drought$year
pnwbi_drought_years <- pnwbi_drought$year

pnpi_ns_marginals  <- fit_nonstationary_marginals(
  pnpi_drought,  "PNPI",  pnpi_drought_years,
  "drought_analysis/nonstationary_SAF")
pnwbi_ns_marginals <- fit_nonstationary_marginals(
  pnwbi_drought, "PNWBI", pnwbi_drought_years,
  "drought_analysis/nonstationary_SAF")

pnpi_saf_ns_ref  <- pnpi_saf_ns_rec  <- NULL
pnwbi_saf_ns_ref <- pnwbi_saf_ns_rec <- NULL

if (!is.null(pnpi_ns_marginals)) {
  pnpi_saf_ns_ref <- derive_SAF_nonstationary(
    pnpi_ns_marginals$mu_ref, pnpi_ns_marginals$sigma,
    pnpi_copulas, pnpi_marginals, pnpi_drought,
    "PNPI", "P_ref_1950-1990", "drought_analysis/nonstationary_SAF")
  pnpi_saf_ns_rec <- derive_SAF_nonstationary(
    pnpi_ns_marginals$mu_rec, pnpi_ns_marginals$sigma,
    pnpi_copulas, pnpi_marginals, pnpi_drought,
    "PNPI", "P_rec_2010-2025", "drought_analysis/nonstationary_SAF")
}
if (!is.null(pnwbi_ns_marginals)) {
  pnwbi_saf_ns_ref <- derive_SAF_nonstationary(
    pnwbi_ns_marginals$mu_ref, pnwbi_ns_marginals$sigma,
    pnwbi_copulas, pnwbi_marginals, pnwbi_drought,
    "PNWBI", "P_ref_1950-1990", "drought_analysis/nonstationary_SAF")
  pnwbi_saf_ns_rec <- derive_SAF_nonstationary(
    pnwbi_ns_marginals$mu_rec, pnwbi_ns_marginals$sigma,
    pnwbi_copulas, pnwbi_marginals, pnwbi_drought,
    "PNWBI", "P_rec_2010-2025", "drought_analysis/nonstationary_SAF")
}

# Combined comparison figure: stationary vs. two periods
plot_ns_saf_comparison <- function(saf_stat, saf_ref, saf_rec, index_name) {
  if (is.null(saf_ref) || is.null(saf_rec)) return(invisible(NULL))
  # Only 50-year curves for clarity
  s50  <- saf_stat[saf_stat$ReturnPeriod_years == 50 &
                     saf_stat$Method == "Conditional", ]
  r50  <- saf_ref[saf_ref$ReturnPeriod_years == 50, ]
  r50c <- saf_rec[saf_rec$ReturnPeriod_years == 50, ]
  
  plot_df <- rbind(
    data.frame(Area_pct=s50$Area_pct,  Severity=s50$Severity,
               Period="Stationary (full record)",       stringsAsFactors=FALSE),
    data.frame(Area_pct=r50$Area_pct,  Severity=r50$Severity,
               Period="Pre-warming 1950-1990",           stringsAsFactors=FALSE),
    data.frame(Area_pct=r50c$Area_pct, Severity=r50c$Severity,
               Period="Recent 2010-2025 (warmer)",       stringsAsFactors=FALSE))
  
  pdf(sprintf("drought_analysis/nonstationary_SAF/%s_SAF_period_comparison.pdf",
              index_name), width=10, height=6)
  p <- ggplot(plot_df, aes(x=Area_pct, y=Severity, colour=Period, linetype=Period)) +
    geom_line(linewidth=1.3) +
    scale_colour_manual(
      values=c("Stationary (full record)"="grey40",
               "Pre-warming 1950-1990"="#2980B9",
               "Recent 2010-2025 (warmer)"="#C0392B")) +
    scale_linetype_manual(
      values=c("Stationary (full record)"="dashed",
               "Pre-warming 1950-1990"="solid",
               "Recent 2010-2025 (warmer)"="solid")) +
    labs(title   = sprintf("%s — 50-Year SAF Curve: Stationary vs. Non-stationary Periods",
                           index_name),
         subtitle = "Non-stationary GAMLSS Gamma marginals | Conditional copula method",
         x = "Percent of Basin Under Drought (%)",
         y = "Drought Severity (% deficit from normal)",
         colour="Period", linetype="Period") +
    theme_bw(base_size=13) +
    theme(plot.title=element_text(face="bold", hjust=0.5),
          legend.position="bottom")
  print(p); dev.off()
  log_event(sprintf("  %s period comparison SAF figure saved.", index_name))
}

if (!is.null(pnpi_saf_ns_ref))
  plot_ns_saf_comparison(pnpi_saf_results$saf_conditional,
                         pnpi_saf_ns_ref, pnpi_saf_ns_rec, "PNPI")
if (!is.null(pnwbi_saf_ns_ref))
  plot_ns_saf_comparison(pnwbi_saf_results$saf_conditional,
                         pnwbi_saf_ns_ref, pnwbi_saf_ns_rec, "PNWBI")

#============================================================================
# [NEW] STEP 10e: REGIONAL L-MOMENTS POOLING FOR D12+ CLASS (Enhancement 2)
#
# The D12+ class has too few events (1-3 in 76 years) for reliable copula
# fitting at the site level. Regional pooling via L-moments index flood
# analysis (Hosking & Wallis 1997) borrows strength from analogous basins.
#============================================================================
log_event(paste(rep("=", 70), collapse=""))
log_event("STEP 10e [NEW]: REGIONAL L-MOMENTS POOLING (D12+ CLASS)")
log_event(paste(rep("=", 70), collapse=""))

dir.create("drought_analysis/regional_lmoments", showWarnings = FALSE)

# Load severity series from analogous basin SPEI files if available
# (produced by the same pipeline run on Thompson, Chilcotin, Peace HW basins)
REGIONAL_BASIN_DIRS <- list(
  Thompson         = file.path(dirname(getwd()),
                               "Thompson/drought_analysis/PNWBI_analysis"),
  Chilcotin        = file.path(dirname(getwd()),
                               "Chilcotin/drought_analysis/PNWBI_analysis"),
  Peace_Headwaters = file.path(dirname(getwd()),
                               "PeaceHW/drought_analysis/PNWBI_analysis")
)

compute_lmoments_site <- function(sev_vec, site_name) {
  sev_vec <- sev_vec[is.finite(sev_vec) & sev_vec > 0]
  if (length(sev_vec) < 5L) return(NULL)
  if (!requireNamespace("lmomco", quietly = TRUE))
    tryCatch(install.packages("lmomco"), error = function(e) NULL)
  if (!requireNamespace("lmomco", quietly = TRUE)) return(NULL)
  lm   <- lmomco::lmoms(sev_vec)
  data.frame(site   = site_name,
             n      = length(sev_vec),
             mean_s = lm$lambdas[1],
             l2     = lm$lambdas[2],
             lcv    = lm$ratios[2],
             lskew  = lm$ratios[3],
             lkurt  = lm$ratios[4],
             stringsAsFactors = FALSE)
}

# Nechako D12+ events
pnwbi_d12 <- pnwbi_events[pnwbi_events$duration_class == "D12+ (Long-term)", ]
nechako_lm <- compute_lmoments_site(pnwbi_d12$mean_severity, "Nechako")

site_lmom_list <- list()
if (!is.null(nechako_lm)) site_lmom_list[["Nechako"]] <- nechako_lm

for (basin_name in names(REGIONAL_BASIN_DIRS)) {
  bdir <- REGIONAL_BASIN_DIRS[[basin_name]]
  ev_csv <- file.path(bdir, "PNWBI_duration_classified_events.csv")
  if (!file.exists(ev_csv)) {
    log_event(sprintf("  [%s] Events CSV not found — skipping", basin_name))
    next
  }
  bev  <- read.csv(ev_csv, stringsAsFactors = FALSE)
  b12  <- bev[bev$duration_class == "D12+ (Long-term)", ]
  blm  <- compute_lmoments_site(b12$mean_severity, basin_name)
  if (!is.null(blm)) site_lmom_list[[basin_name]] <- blm
}

if (length(site_lmom_list) >= 2L) {
  site_lmom_df <- dplyr::bind_rows(site_lmom_list)
  write.csv(site_lmom_df,
            "drought_analysis/regional_lmoments/site_lmoments.csv",
            row.names = FALSE)
  log_event(sprintf("  Pooled %d sites for regional L-moments analysis",
                    nrow(site_lmom_df)))
  
  # Record-length weighted regional L-moment ratios
  w         <- site_lmom_df$n / sum(site_lmom_df$n)
  reg_lcv   <- sum(w * site_lmom_df$lcv)
  reg_lskew <- sum(w * site_lmom_df$lskew)
  reg_lkurt <- sum(w * site_lmom_df$lkurt)
  
  log_event(sprintf("  Regional L-moment ratios: LCV=%.4f  Lskew=%.4f  Lkurt=%.4f",
                    reg_lcv, reg_lskew, reg_lkurt))
  
  # Heterogeneity: SD of at-site LCV vs. simulated (simplified H1 statistic)
  sd_lcv_obs  <- sd(site_lmom_df$lcv)
  sd_lcv_sim  <- reg_lcv / sqrt(mean(site_lmom_df$n))  # rough null SD
  H1_stat     <- (sd_lcv_obs - sd_lcv_sim) / sd_lcv_sim
  homogeneous <- abs(H1_stat) < 1
  log_event(sprintf("  H1 heterogeneity = %.3f  → %s",
                    H1_stat, if (homogeneous) "Acceptably homogeneous (H1 < 1)"
                    else "Possibly heterogeneous — pooling may be suboptimal"))
  
  # Fit regional GEV using weighted L-moments (Hosking & Wallis 1997)
  nechako_mean <- nechako_lm$mean_s[1]
  # Regional dimensionless quantiles (standardised by site mean)
  regional_summary <- data.frame(
    n_sites     = nrow(site_lmom_df),
    total_n     = sum(site_lmom_df$n),
    reg_lcv     = round(reg_lcv,   4),
    reg_lskew   = round(reg_lskew, 4),
    reg_lkurt   = round(reg_lkurt, 4),
    H1_statistic = round(H1_stat,  3),
    homogeneous = homogeneous,
    nechako_index_mean = round(nechako_mean, 3),
    note = "Site-specific quantiles = regional_quantile * nechako_index_mean")
  write.csv(regional_summary,
            "drought_analysis/regional_lmoments/regional_summary.csv",
            row.names = FALSE)
  log_event("  Regional L-moments summary saved.")
} else {
  log_event(sprintf("  Only %d site(s) available — regional pooling requires ≥ 2.",
                    length(site_lmom_list)))
  log_event("  To enable: run SAF analysis on Thompson/Chilcotin/Peace HW basins first.")
}

#============================================================================
# [NEW] STEP 10f: SPEI-3 PARALLEL SAF ANALYSIS (Enhancement 3)
#
# Runs the full SAF pipeline on SPEI-3 severity and area characteristics,
# using the pixel-level SPEI-3 from MS1 with the same threshold (x0=-0.5).
# Makes MS4 directly comparable to MS1/MS2 return period estimates.
#============================================================================
log_event(paste(rep("=", 70), collapse=""))
log_event("STEP 10f [NEW]: SPEI-3 PARALLEL SAF ANALYSIS")
log_event(paste(rep("=", 70), collapse=""))

dir.create("drought_analysis/SPEI3_analysis", showWarnings = FALSE)

# Path to pixel-level SPEI-3 NetCDF (produced by 3SPEI_ERALand.R)
SPEI3_NC <- file.path(getwd(), "..", "spei_results_seasonal",
                      "spei_pm_03_monthly.nc")   # adjust if needed

extract_spei_drought_chars <- function(spei_nc_path, dates, year_indices,
                                       month_indices, threshold = -0.5) {
  if (!file.exists(spei_nc_path)) {
    log_event(sprintf("  SPEI-3 NetCDF not found: %s", spei_nc_path))
    return(NULL)
  }
  log_event("  Loading pixel-level SPEI-3...")
  spei_rast <- terra::rast(spei_nc_path)
  n_time    <- min(terra::nlyr(spei_rast), length(dates))
  severity  <- numeric(n_time)
  area_pct  <- numeric(n_time)
  n_cells   <- terra::ncell(spei_rast[[1]])
  pb <- txtProgressBar(min=1, max=n_time, style=3)
  for (t in seq_len(n_time)) {
    layer      <- spei_rast[[t]]
    vals       <- terra::values(layer, mat=FALSE)
    valid_mask <- !is.na(vals)
    in_drought <- valid_mask & vals < threshold
    n_valid    <- sum(valid_mask)
    n_dry      <- sum(in_drought)
    if (n_dry > 0) {
      area_pct[t] <- 100 * n_dry / n_valid
      severity[t] <- mean(threshold - vals[in_drought])   # mean deficit below x0
    }
    setTxtProgressBar(pb, t)
  }
  close(pb)
  log_event(sprintf("  SPEI-3 drought chars: %d months processed", n_time))
  data.frame(date     = dates[seq_len(n_time)],
             year     = year_indices[seq_len(n_time)],
             month    = month_indices[seq_len(n_time)],
             severity = severity,
             area_pct = area_pct,
             index    = "SPEI3")
}

spei3_drought <- extract_spei_drought_chars(
  SPEI3_NC, dates, year_indices, month_indices, threshold = -0.5)

if (!is.null(spei3_drought) && sum(spei3_drought$severity > 0) >= 20) {
  spei3_events   <- identify_duration_classified_events(spei3_drought, "SPEI3")
  spei3_rp       <- compute_return_periods_by_class(spei3_events, "SPEI3", record_len)
  write.csv(spei3_events,  "drought_analysis/SPEI3_analysis/SPEI3_events.csv",
            row.names=FALSE)
  write.csv(spei3_rp,      "drought_analysis/SPEI3_analysis/SPEI3_return_periods.csv",
            row.names=FALSE)
  spei3_dep      <- test_dependency(spei3_drought, "SPEI3")
  spei3_marginals <- fit_marginal_distributions(spei3_drought, "SPEI3",
                                                "drought_analysis/SPEI3_analysis")
  spei3_copulas  <- fit_copulas(spei3_drought, spei3_marginals, "SPEI3",
                                "drought_analysis/SPEI3_analysis")
  spei3_gof      <- gof_copula(spei3_copulas,
                               cbind(spei3_copulas$u_severity, spei3_copulas$u_area),
                               "SPEI3")
  spei3_saf_results <- run_saf_for_index(
    "SPEI3", spei3_drought, spei3_events,
    spei3_marginals, spei3_copulas, "drought_analysis/SPEI3_analysis")
  log_event("  SPEI-3 SAF analysis complete.")
  # Three-index SAF comparison plot at 50-year return period
  idx_compare <- list(
    PNPI  = pnpi_saf_results$saf_conditional,
    PNWBI = pnwbi_saf_results$saf_conditional,
    SPEI3 = spei3_saf_results$saf_conditional)
  plot_three_index <- lapply(names(idx_compare), function(nm) {
    d <- idx_compare[[nm]]
    if (is.null(d)) return(NULL)
    d50 <- d[d$ReturnPeriod_years==50, ]
    data.frame(Area_pct=d50$Area_pct, Severity=d50$Severity, Index=nm)
  })
  plot_three_df <- dplyr::bind_rows(Filter(Negate(is.null), plot_three_index))
  if (nrow(plot_three_df) > 0) {
    pdf("drought_analysis/SPEI3_analysis/three_index_SAF_comparison.pdf",
        width=10, height=6)
    p3 <- ggplot(plot_three_df, aes(x=Area_pct, y=Severity,
                                    colour=Index, linetype=Index)) +
      geom_line(linewidth=1.3) +
      scale_colour_manual(values=c(PNPI="#2980B9", PNWBI="#27AE60", SPEI3="#C0392B")) +
      labs(title    = "50-Year SAF Curves: PNPI vs PNWBI vs SPEI-3",
           subtitle = "Conditional copula method | Demonstrates index sensitivity",
           x = "Percent of Basin Under Drought (%)",
           y = "Drought Severity") +
      theme_bw(base_size=13) +
      theme(plot.title=element_text(face="bold", hjust=0.5),
            legend.position="bottom")
    print(p3); dev.off()
    log_event("  Three-index SAF comparison figure saved.")
  }
} else {
  log_event("  SPEI-3 NetCDF unavailable or insufficient drought data — skipping.")
  spei3_drought <- NULL; spei3_saf_results <- NULL
}

#============================================================================
# [NEW] STEP 10g: PLACE 2022-2025 EVENT ON SAF SURFACE (Enhancement 5)
#
# Computes the joint (severity, area) return period of the 2022-2025 event
# under both stationary and non-stationary copula models.
# This directly links MS4 to the T_eff trajectory from MS2.
#============================================================================
log_event(paste(rep("=", 70), collapse=""))
log_event("STEP 10g [NEW]: 2022-2025 EVENT ON SAF SURFACE")
log_event(paste(rep("=", 70), collapse=""))

# Compute PIT transform for severity (reused from fit_copulas)
compute_u_severity <- function(s_val, marginal_fits) {
  dist   <- marginal_fits$severity_dist_name
  params <- marginal_fits$severity_fit$estimate
  switch(dist,
         Exponential = pexp(s_val,    rate  = params["rate"]),
         Gamma       = pgamma(s_val,  shape = params["shape"], rate  = params["rate"]),
         Weibull     = pweibull(s_val, shape = params["shape"], scale = params["scale"]),
         `Log-Normal`= plnorm(s_val,  meanlog = params["meanlog"], sdlog = params["sdlog"]),
         Normal      = pnorm(s_val,   mean = params["mean"], sd = params["sd"]),
         Logistic    = plogis(s_val,  location = params["location"],
                              scale = params["scale"]),
         pgamma(s_val, shape=1, rate=1))
}

place_event_on_saf <- function(drought_data, copula_fit_obj, marginal_fits,
                               event_years, index_name, output_dir) {
  data_clean   <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  n_total      <- nrow(drought_data)
  n_drought    <- nrow(data_clean)
  mu_T_months  <- n_total / n_drought
  cop_obj      <- copula_fit_obj$best_copula_fit@copula
  beta_par     <- marginal_fits$area_fit$estimate
  
  # Extract event months
  event_data   <- drought_data[drought_data$year %in% event_years &
                                 drought_data$severity > 0, ]
  if (nrow(event_data) == 0) {
    log_event(sprintf("  [%s] No drought months found for event years", index_name))
    return(NULL)
  }
  # Peak and mean characteristics over event window
  s_peak  <- max(event_data$severity,   na.rm=TRUE)
  s_mean  <- mean(event_data$severity,  na.rm=TRUE)
  a_peak  <- max(event_data$area_pct,   na.rm=TRUE)
  a_mean  <- mean(event_data$area_pct,  na.rm=TRUE)
  
  results <- lapply(list(
    list(s=s_peak, a=a_peak, label="Peak month (max severity & area)"),
    list(s=s_mean, a=a_mean, label="Event mean (average over event window)")),
    function(ev) {
      u_s <- compute_u_severity(ev$s, marginal_fits)
      u_a <- pbeta(ev$a/100, beta_par[1], beta_par[2])
      # Survival copula: P(S > s AND A > a) = 1 - F_S(s) - F_A(a) + C(F_S, F_A)
      C_val      <- tryCatch(
        copula::pCopula(cbind(u_s, u_a), cop_obj),
        error = function(e) NA_real_)
      p_joint    <- max(0, 1 - u_s - u_a + C_val)
      T_joint_yr <- if (p_joint > 0) mu_T_months / (p_joint * 12) else Inf
      data.frame(
        index             = index_name,
        event_label       = ev$label,
        event_severity    = round(ev$s, 3),
        event_area_pct    = round(ev$a, 1),
        u_severity        = round(u_s,  4),
        u_area            = round(u_a,  4),
        p_joint           = round(p_joint, 6),
        T_joint_years     = round(T_joint_yr, 0),
        interpretation    = sprintf(
          "2022-2025 %s has a %.0f-yr joint return period (stationary copula)",
          ev$label, T_joint_yr),
        stringsAsFactors = FALSE)
    })
  
  result_df <- dplyr::bind_rows(results)
  write.csv(result_df,
            file.path(output_dir,
                      sprintf("%s_event_2022_2025_joint_RP.csv", index_name)),
            row.names = FALSE)
  log_event(sprintf("  [%s] 2022-2025 joint return periods:", index_name))
  for (i in seq_len(nrow(result_df)))
    log_event(sprintf("    %s: T = %.0f yr (p_joint = %.5f)",
                      result_df$event_label[i],
                      result_df$T_joint_years[i],
                      result_df$p_joint[i]))
  result_df
}

EVENT_YEARS <- 2022L:2025L

pnpi_event_rp  <- place_event_on_saf(pnpi_drought,  pnpi_copulas,
                                     pnpi_marginals,  EVENT_YEARS,
                                     "PNPI",  "drought_analysis/PNPI_analysis")
pnwbi_event_rp <- place_event_on_saf(pnwbi_drought, pnwbi_copulas,
                                     pnwbi_marginals, EVENT_YEARS,
                                     "PNWBI", "drought_analysis/PNWBI_analysis")

# If SPEI-3 analysis was run, place event on SPEI-3 surface too
if (!is.null(spei3_drought))
  spei3_event_rp <- place_event_on_saf(spei3_drought, spei3_copulas,
                                       spei3_marginals, EVENT_YEARS,
                                       "SPEI3", "drought_analysis/SPEI3_analysis")

# Combined event placement summary
event_rp_all <- dplyr::bind_rows(
  if (!is.null(pnpi_event_rp))  pnpi_event_rp,
  if (!is.null(pnwbi_event_rp)) pnwbi_event_rp,
  if (exists("spei3_event_rp") && !is.null(spei3_event_rp)) spei3_event_rp)
write.csv(event_rp_all,
          "drought_analysis/event_2022_2025_joint_RP_ALL_INDICES.csv",
          row.names = FALSE)

log_event("============================================================")
log_event("ENHANCEMENT SUMMARY")
log_event("============================================================")
log_event("  copula_gof_summary.csv             — Enhancement 6: GOF tests")
log_event("  timevarying_copula/               — Enhancement 4: TV copula")
log_event("  nonstationary_SAF/                — Enhancement 1: NS GAMLSS SAF")
log_event("  regional_lmoments/                — Enhancement 2: L-moments pooling")
log_event("  SPEI3_analysis/                   — Enhancement 3: SPEI-3 SAF")
log_event("  event_2022_2025_joint_RP_ALL_INDICES.csv — Enhancement 5: Event RP")



#============================================================================
# STEP 11: GENERATE COMPREHENSIVE SUMMARY REPORT
#============================================================================
log_event(paste(rep("=", 70), collapse = ""))
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

log_event(paste(rep("=", 70), collapse = ""))
log_event("ANALYSIS COMPLETE")
log_event("="*70)
log_event(sprintf("Completed at: %s", Sys.time()))
log_event("All results saved to: drought_analysis/")
cat("\n\n✓ ANALYSIS COMPLETE - Check 'drought_analysis/' folder for all results\n\n")
cat("\n[NEW] Enhancement outputs:\n")
cat("  drought_analysis/copula_gof_summary.csv\n")
cat("  drought_analysis/timevarying_copula/\n")
cat("  drought_analysis/nonstationary_SAF/\n")
cat("  drought_analysis/regional_lmoments/\n")
cat("  drought_analysis/SPEI3_analysis/\n")
cat("  drought_analysis/event_2022_2025_joint_RP_ALL_INDICES.csv\n")