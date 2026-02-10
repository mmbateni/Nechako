
# spi_gamma.R
# Compute Standardized Precipitation Index (SPI) at user-defined timescales
# using a Gamma distribution (maximum-likelihood fit), per McKee/Edwards.

# ---- Dependencies ----
# install.packages(c("zoo", "MASS"))
library(zoo)   # rollapply for moving sums
library(MASS)  # fitdistr for Gamma MLE

# ---- Helper: rolling accumulation (end-aligned) ----
accumulate_series <- function(x, scale) {
  rollapply(x, width = scale, by = 1, align = "right",
            FUN = sum, na.rm = FALSE, partial = FALSE)
}

# ---- Helper: month factor from Date vector ----
month_factor <- function(dates) as.integer(format(dates, "%m"))

# ---- Fit Gamma per calendar month and transform to SPI ----
spi_gamma <- function(prcp_monthly, dates,
                      scales = c(1, 2, 3, 6, 9, 12, 24, 36, 48, 60, 72),
                      baseline_start = NULL, baseline_end = NULL,
                      min_samples_per_month = 30) {

  stopifnot(length(prcp_monthly) == length(dates))
  z <- zoo(prcp_monthly, order.by = as.Date(dates))

  # Baseline mask
  if (!is.null(baseline_start) && !is.null(baseline_end)) {
    baseline_idx <- index(z) >= as.Date(baseline_start) & index(z) <= as.Date(baseline_end)
  } else {
    baseline_idx <- rep(TRUE, length(z))
  }

  out <- list()

  for (k in scales) {
    ak <- accumulate_series(z, k)         # k-month totals
    idx <- index(ak)
    mon <- month_factor(idx)

    spi_vals <- rep(NA_real_, length(ak))

    for (m in 1:12) {
      sel_all  <- which(mon == m)
      # Baseline subset within month m
      if (!is.null(baseline_start) && !is.null(baseline_end)) {
        sel_base <- sel_all[idx[sel_all] >= as.Date(baseline_start) & idx[sel_all] <= as.Date(baseline_end)]
      } else {
        sel_base <- sel_all
      }

      x_base <- coredata(ak)[sel_base]
      x_obs  <- coredata(ak)[sel_all]

      # Ensure enough baseline samples per calendar month
      if (sum(!is.na(x_base)) < min_samples_per_month) next

      # Fraction of zeros in baseline (q)
      n_base <- sum(!is.na(x_base))
      m_zero <- sum(x_base == 0, na.rm = TRUE)
      q_zero <- m_zero / n_base

      # Fit Gamma on strictly positive baseline values
      x_pos <- x_base[x_base > 0 & !is.na(x_base)]
      if (length(x_pos) < 10) next
      fit <- try(fitdistr(x_pos, densfun = "gamma"), silent = TRUE)
      if (inherits(fit, "try-error")) next

      shape <- fit$estimate[["shape"]]
      rate  <- fit$estimate[["rate"]]

      # CDF mixing with zeros: H(x) = q + (1-q)*G(x)
      Gx <- rep(NA_real_, length(x_obs))
      pos_mask <- !is.na(x_obs)
      Gx[pos_mask] <- pgamma(pmax(x_obs[pos_mask], 0), shape = shape, rate = rate)
      Hx <- q_zero + (1 - q_zero) * Gx

      # SPI = inverse-normal
      spi_vals[sel_all] <- qnorm(Hx)
    }

    out[[paste0("SPI_", k, "mo_gamma")]] <- zoo(spi_vals, order.by = idx)
  }

  return(out)
}

# ---- Example ----
# df <- data.frame(date = seq(as.Date("1980-01-01"), as.Date("2024-12-01"), by = "month"),
#                  prcp = pmax(rnorm(541, mean = 60, sd = 25), 0))
# res_gamma <- spi_gamma(prcp_monthly = df$prcp, dates = df$date,
#                        baseline_start = "1991-01-01", baseline_end = "2020-12-01")
# plot(res_gamma$SPI_12mo_gamma, main = "SPI-12 (Gamma)", ylab = "SPI", xlab = "Time")
``
