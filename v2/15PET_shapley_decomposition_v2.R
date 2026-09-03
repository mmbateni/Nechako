################################################################################
# 15PET_shapley_decomposition_v2.R
# SHAPLEY DECOMPOSITION OF PENMAN-MONTEITH PET ANOMALIES
#
# Purpose
#   Decompose the monthly basin-mean Penman-Monteith PET departure from its
#   calendar-month climatological reference state into exact Shapley allocations
#   across the meteorological inputs used by the PM equation.
#
# IMPORTANT INTERPRETATION
#   This is an additive decomposition of the PM PET model output anomaly,
#   NOT a decomposition of SPEI_PM, SPEI_Thw, or delta_SPEI_PET.
#   delta_SPEI_PET remains a PET-formulation sensitivity residual.
#
# Shapley inputs
#   Temperature, net radiation, VPD, 2-m wind speed, surface pressure.
#   Net radiation is treated as an already-computed PM input; therefore the
#   temperature allocation is conditional on the diagnosed net-radiation field
#   rather than an attribution of every upstream process that generated Rn.
#
# Reference state
#   Calendar-month climatology over 1950-2025 for each input.
#   This removes the seasonal cycle and asks how the observed monthly PM PET
#   differs from the expected PM PET under a typical state for that month.
#
# Method
#   Exact permutation Shapley values. Five inputs require 5! = 120 permutations
#   per month, so no Monte-Carlo SHAP approximation is required.
#
# Inputs
#   monthly_data_direct/meteorological_controls_monthly_summary.csv
#   monthly_data_direct/ERA5Land_Nechako_PET_monthly_summary.csv
#   pet_sensitivity/PET_sensitivity_monthly.csv
#
# Outputs -> pet_sensitivity/shapley/
#   PET_PM_Shapley_monthly.csv
#   PET_PM_Shapley_summary.csv
#   PET_PM_Shapley_2022_2025.csv
#   Fig_PET_PM_Shapley_2022_2025.pdf/.png
#   Fig_PET_PM_Shapley_annual.pdf/.png
################################################################################

setwd(Sys.getenv("NECHAKO_WD", "D:/Nechako_Drought/Nechako/"))
source("DROUGHT_ANALYSIS_utils_v2.R")
utils_load_packages(c("dplyr", "tidyr", "ggplot2", "data.table", "lubridate"))

if (!dir.exists(WD_PATH)) stop("Working directory not found: ", WD_PATH)
setwd(WD_PATH)

OUT_DIR <- file.path(WD_PATH, "pet_sensitivity", "shapley")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

ANALYSIS_START <- as.Date("1950-01-01")
ANALYSIS_END   <- as.Date("2025-12-01")
RECENT_START   <- as.Date("2022-01-01")
RECENT_END     <- as.Date("2025-12-01")

MET_FILE <- file.path(WD_PATH, "monthly_data_direct",
                      "meteorological_controls_monthly_summary.csv")
PET_FILE <- file.path(WD_PATH, "monthly_data_direct",
                      "ERA5Land_Nechako_PET_monthly_summary.csv")
SENS_FILE <- file.path(WD_PATH, "pet_sensitivity",
                       "PET_sensitivity_monthly.csv")

for (f in c(MET_FILE, PET_FILE, SENS_FILE))
  if (!file.exists(f)) stop("Required input not found: ", f)

met <- data.table::fread(MET_FILE)
pet <- data.table::fread(PET_FILE)

required_met <- c("date", "year", "month", "t2m_degC", "rn_MJ_m2_d",
                  "vpd_kPa", "wind2m_ms", "pressure_kPa")
miss <- setdiff(required_met, names(met))
if (length(miss)) {
  stop("Meteorological summary is missing pressure_kPa or another required field: ",
       paste(miss, collapse = ", "),
       "\nRun 2b_PET_ERALand_v3.R once to regenerate the summary CSV with pressure_kPa.")
}

met$date <- as.Date(paste0(as.character(met$date), "-01"))
pet$date <- as.Date(paste0(as.character(pet$date), "-01"))

validate_monthly_record(met$date, ANALYSIS_START, ANALYSIS_END,
                        label = "meteorological-controls record",
                        require_complete = TRUE)
validate_monthly_record(pet$date, ANALYSIS_START, ANALYSIS_END,
                        label = "PM PET summary record",
                        require_complete = TRUE)

met <- met |>
  dplyr::select(date, year, month, t2m_degC, rn_MJ_m2_d, vpd_kPa,
                wind2m_ms, pressure_kPa) |>
  dplyr::arrange(date)

pet <- pet |>
  dplyr::select(date, mean_pet) |>
  dplyr::rename(pm_pet_summary_mm_d = mean_pet)

base <- met |>
  dplyr::left_join(pet, by = "date") |>
  dplyr::arrange(date) |>
  dplyr::mutate(
    year = as.integer(format(date, "%Y")),
    month = as.integer(format(date, "%m")),
    recent = date >= RECENT_START & date <= RECENT_END,
    warm_season = month %in% 5:9
  )

if (anyNA(base$pm_pet_summary_mm_d))
  stop("PM PET summary failed to match all meteorological-control dates.")

# ------------------------------------------------------------------------------
# PM MODEL FUNCTION
#
# This reproduces the FAO-56 monthly PM algebra used in 2b at the basin-input
# level. The current-month temperature and the previous-month temperature enter
# the soil-heat-flux term G. When temperature is switched between observed and
# reference states, both are switched consistently.
# ------------------------------------------------------------------------------

pm_pet_from_state <- function(T, Rn, VPD, wind, pressure, T_prev) {
  if (any(!is.finite(c(T, Rn, VPD, wind, pressure, T_prev)))) return(NA_real_)
  es <- 0.6108 * exp((17.27 * T) / (T + 237.3))
  delta_slope <- (4098 * es) / (T + 237.3)^2
  gamma <- 0.000665 * pressure
  
  G <- 0.14 * (T - T_prev)
  if (T < 0) {
    G_cap <- 0.3 * abs(Rn)
    G <- max(-G_cap, min(G, G_cap))
  }
  
  num <- 0.408 * delta_slope * (Rn - G) +
    gamma * (900 / (T + 273)) * wind * VPD
  den <- delta_slope + gamma * (1 + 0.34 * wind)
  et0 <- num / den
  
  max(0, as.numeric(et0))
}

# ------------------------------------------------------------------------------
# CALENDAR-MONTH REFERENCE STATES
# ------------------------------------------------------------------------------

clim <- met |>
  dplyr::group_by(month) |>
  dplyr::summarise(
    T_ref = mean(t2m_degC, na.rm = TRUE),
    Rn_ref = mean(rn_MJ_m2_d, na.rm = TRUE),
    VPD_ref = mean(vpd_kPa, na.rm = TRUE),
    wind_ref = mean(wind2m_ms, na.rm = TRUE),
    pressure_ref = mean(pressure_kPa, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# EXACT PERMUTATION SHAPLEY
# ------------------------------------------------------------------------------

players <- c("T", "Rn", "VPD", "wind", "pressure")

# Generate all n! permutations recursively.
# Five players -> exactly 120 permutations.
permute <- function(x) {
  if (length(x) <= 1L) {
    return(list(x))
  }
  
  out <- vector("list", factorial(length(x)))
  k <- 0L
  
  for (i in seq_along(x)) {
    rest <- x[-i]
    sub <- permute(rest)
    
    for (z in sub) {
      k <- k + 1L
      out[[k]] <- c(x[i], z)
    }
  }
  
  out
}

# The recursive function returns a list of character vectors.  Bind that list
# directly into the permutation matrix.  The previous implementation applied
# do.call(rbind, ...) twice: the inner call already returned a matrix, so the
# outer call failed with "second argument must be a list".
perm_list <- permute(players)
perm_matrix <- do.call(rbind, perm_list)
storage.mode(perm_matrix) <- "character"

if (!is.matrix(perm_matrix) ||
    nrow(perm_matrix) != factorial(length(players)) ||
    ncol(perm_matrix) != length(players)) {
  stop(
    "Permutation construction failed: expected ",
    factorial(length(players)), " x ", length(players),
    " permutation matrix; got ",
    paste(dim(perm_matrix), collapse = " x ")
  )
}

# Every row must contain every player exactly once.
if (any(apply(perm_matrix, 1L, function(z) !setequal(z, players)))) {
  stop("Permutation matrix contains an invalid permutation.")
}

state_value <- function(state, obs, ref, prev_obs_T, prev_ref_T) {
  vals <- ref
  if ("T" %in% state) vals$T <- obs$T
  if ("Rn" %in% state) vals$Rn <- obs$Rn
  if ("VPD" %in% state) vals$VPD <- obs$VPD
  if ("wind" %in% state) vals$wind <- obs$wind
  if ("pressure" %in% state) vals$pressure <- obs$pressure
  tprev <- if ("T" %in% state) prev_obs_T else prev_ref_T
  pm_pet_from_state(vals$T, vals$Rn, vals$VPD, vals$wind,
                    vals$pressure, tprev)
}

shapley_one <- function(i, d, clim_tbl) {
  m <- d$month[i]
  cr <- clim_tbl[clim_tbl$month == m, , drop = FALSE]
  if (!nrow(cr)) stop("No climatological reference for month ", m)
  
  # Previous-month reference temperature is the climatological mean of the
  # previous calendar month; December follows November within the monthly
  # climatology, so no year-boundary special case is needed.
  prev_m <- if (m == 1L) 12L else m - 1L
  cr_prev <- clim_tbl$T_ref[match(prev_m, clim_tbl$month)]
  
  if (i == 1L) {
    prev_obs_T <- d$t2m_degC[i]
  } else {
    prev_obs_T <- d$t2m_degC[i - 1L]
  }
  
  obs <- list(T=d$t2m_degC[i], Rn=d$rn_MJ_m2_d[i], VPD=d$vpd_kPa[i],
              wind=d$wind2m_ms[i], pressure=d$pressure_kPa[i])
  ref <- list(T=cr$T_ref, Rn=cr$Rn_ref, VPD=cr$VPD_ref,
              wind=cr$wind_ref, pressure=cr$pressure_ref)
  
  if (any(!is.finite(unlist(obs))) || any(!is.finite(unlist(ref))) ||
      !is.finite(prev_obs_T) || !is.finite(cr_prev)) {
    return(rep(NA_real_, length(players)))
  }
  
  phi <- setNames(rep(0, length(players)), players)
  for (p in players) {
    vals <- numeric(nrow(perm_matrix))
    for (r in seq_len(nrow(perm_matrix))) {
      perm <- perm_matrix[r, ]
      pos <- match(p, perm)
      S_before <- if (pos > 1L) perm[seq_len(pos - 1L)] else character(0)
      S_after <- perm[seq_len(pos)]
      v_before <- state_value(S_before, obs, ref, prev_obs_T, cr_prev)
      v_after  <- state_value(S_after,  obs, ref, prev_obs_T, cr_prev)
      vals[r] <- v_after - v_before
    }
    phi[p] <- mean(vals, na.rm = TRUE)
  }
  unname(phi)
}

message("Computing exact five-player Shapley allocations for ", nrow(base), " months...")
phi_mat <- t(vapply(seq_len(nrow(base)), shapley_one,
                    numeric(length(players)), d=base, clim_tbl=clim))
colnames(phi_mat) <- paste0("phi_", players)

shapley <- cbind(base, as.data.frame(phi_mat, check.names = FALSE))
shapley <- shapley |>
  dplyr::mutate(
    pm_pet_reconstructed_mm_d = mapply(
      pm_pet_from_state, t2m_degC, rn_MJ_m2_d, vpd_kPa,
      wind2m_ms, pressure_kPa,
      dplyr::lag(t2m_degC, default = t2m_degC[1])
    ),
    pm_pet_reference_mm_d = vapply(seq_len(nrow(shapley)), function(i) {
      cr <- clim[clim$month == shapley$month[i], , drop=FALSE]
      prev_m <- if (shapley$month[i] == 1L) 12L else shapley$month[i] - 1L
      prev_ref <- clim$T_ref[match(prev_m, clim$month)]
      pm_pet_from_state(cr$T_ref, cr$Rn_ref, cr$VPD_ref,
                        cr$wind_ref, cr$pressure_ref, prev_ref)
    }, numeric(1)),
    pm_pet_anomaly_mm_d = pm_pet_reconstructed_mm_d - pm_pet_reference_mm_d,
    shapley_sum_mm_d = rowSums(as.data.frame(phi_mat), na.rm = TRUE),
    shapley_closure_error_mm_d = shapley_sum_mm_d - pm_pet_anomaly_mm_d,
    wb_phi_T = -phi_T,
    wb_phi_Rn = -phi_Rn,
    wb_phi_VPD = -phi_VPD,
    wb_phi_wind = -phi_wind,
    wb_phi_pressure = -phi_pressure
  )

closure_tol <- 1e-8
bad <- is.finite(shapley$shapley_closure_error_mm_d) &
  abs(shapley$shapley_closure_error_mm_d) > closure_tol
if (any(bad))
  stop("Shapley closure failed for ", sum(bad), " month(s); max error = ",
       max(abs(shapley$shapley_closure_error_mm_d[bad]), na.rm=TRUE), " mm/day")

write.csv(shapley, file.path(OUT_DIR, "PET_PM_Shapley_monthly.csv"), row.names=FALSE)

# ------------------------------------------------------------------------------
# SUMMARY TABLES
# ------------------------------------------------------------------------------

summary_long <- shapley |>
  dplyr::select(date, year, month, recent, warm_season,
                pm_pet_anomaly_mm_d, starts_with("phi_")) |>
  tidyr::pivot_longer(starts_with("phi_"), names_to="component", values_to="phi_mm_d") |>
  dplyr::mutate(component=sub("^phi_", "", component),
                component= factor(component, levels=players,
                                  labels=c("Temperature","Net radiation","VPD",
                                           "Wind speed","Surface pressure")))

summary_periods <- list(
  "Full record" = rep(TRUE, nrow(shapley)),
  "2022-2025" = shapley$recent,
  "2022-2025 warm season" = shapley$recent & shapley$warm_season,
  "Historical warm season" = !shapley$recent & shapley$warm_season
)

summary_tbl <- dplyr::bind_rows(lapply(names(summary_periods), function(nm) {
  ii <- summary_periods[[nm]] & is.finite(shapley$pm_pet_anomaly_mm_d)
  z <- summary_long[summary_long$date %in% shapley$date[ii], , drop=FALSE]
  z |>
    dplyr::group_by(component) |>
    dplyr::summarise(
      period = nm,
      n_months = sum(is.finite(phi_mm_d)),
      mean_phi_mm_d = mean(phi_mm_d, na.rm=TRUE),
      mean_abs_phi_mm_d = mean(abs(phi_mm_d), na.rm=TRUE),
      sd_phi_mm_d = sd(phi_mm_d, na.rm=TRUE),
      .groups="drop"
    )
}))

write.csv(summary_tbl, file.path(OUT_DIR, "PET_PM_Shapley_summary.csv"), row.names=FALSE)

recent_tbl <- shapley |>
  dplyr::filter(recent) |>
  dplyr::select(date, year, month, pm_pet_summary_mm_d,
                pm_pet_reconstructed_mm_d, pm_pet_reference_mm_d,
                pm_pet_anomaly_mm_d, phi_T, phi_Rn, phi_VPD,
                phi_wind, phi_pressure, shapley_closure_error_mm_d)
write.csv(recent_tbl, file.path(OUT_DIR, "PET_PM_Shapley_2022_2025.csv"), row.names=FALSE)

# Link the PET-level Shapley allocations to the existing formulation-sensitivity
# residual without treating that residual as decomposable. This table asks only
# whether months with large PM PET allocations also tend to have large |Delta SPEI|.
sens_link <- data.table::fread(SENS_FILE) |>
  dplyr::mutate(date = as.Date(date)) |>
  dplyr::select(date, scale, abs_delta_spei_pet) |>
  dplyr::left_join(
    shapley |>
      dplyr::select(date, phi_T, phi_Rn, phi_VPD, phi_wind, phi_pressure,
                    pm_pet_anomaly_mm_d),
    by = "date"
  )

shapley_sensitivity_link <- dplyr::bind_rows(lapply(split(sens_link, sens_link$scale), function(df) {
  df <- df[is.finite(df$abs_delta_spei_pet), , drop=FALSE]
  comps <- c("phi_T","phi_Rn","phi_VPD","phi_wind","phi_pressure")
  dplyr::bind_rows(lapply(comps, function(cc) {
    ok <- is.finite(df[[cc]]) & is.finite(df$abs_delta_spei_pet)
    data.frame(
      scale = unique(df$scale), component = cc, n = sum(ok),
      pearson_r = if(sum(ok)>=3) cor(df[[cc]][ok], df$abs_delta_spei_pet[ok]) else NA_real_,
      spearman_rho = if(sum(ok)>=3) suppressWarnings(cor(df[[cc]][ok], df$abs_delta_spei_pet[ok], method="spearman")) else NA_real_,
      mean_abs_phi_mm_d = if(sum(ok)) mean(abs(df[[cc]][ok])) else NA_real_
    )
  }))
}))
write.csv(shapley_sensitivity_link,
          file.path(OUT_DIR, "PET_PM_Shapley_vs_PET_sensitivity.csv"),
          row.names=FALSE)

# ------------------------------------------------------------------------------
# FIGURES
# ------------------------------------------------------------------------------

plot_recent <- summary_long |>
  dplyr::filter(date >= RECENT_START & date <= RECENT_END)

p_recent <- ggplot2::ggplot(plot_recent,
                            ggplot2::aes(x=date, y=phi_mm_d, linetype=component)) +
  ggplot2::geom_hline(yintercept=0, linewidth=0.35) +
  ggplot2::geom_line(linewidth=0.45) +
  ggplot2::labs(
    x=NULL, y="Shapley allocation to PM PET anomaly (mm/day)",
    linetype="PM input",
    title="Meteorological Shapley decomposition of PM PET anomaly",
    subtitle="2022–2025; reference = 1950–2025 calendar-month climatology",
    caption="Positive values increase PM PET relative to the monthly climatological reference; negative values decrease it."
  ) +
  ggplot2::theme_bw(base_size=10) +
  ggplot2::theme(legend.position="bottom")

ggplot2::ggsave(file.path(OUT_DIR,"Fig_PET_PM_Shapley_2022_2025.png"),p_recent,
                width=11,height=6.5,dpi=FIG_DPI)
ggplot2::ggsave(file.path(OUT_DIR,"Fig_PET_PM_Shapley_2022_2025.pdf"),p_recent,
                width=11,height=6.5)

annual <- shapley |>
  dplyr::filter(year >= 1950, year <= 2025) |>
  dplyr::group_by(year) |>
  dplyr::summarise(
    phi_T=mean(phi_T,na.rm=TRUE), phi_Rn=mean(phi_Rn,na.rm=TRUE),
    phi_VPD=mean(phi_VPD,na.rm=TRUE), phi_wind=mean(phi_wind,na.rm=TRUE),
    phi_pressure=mean(phi_pressure,na.rm=TRUE),
    .groups="drop"
  ) |>
  tidyr::pivot_longer(-year,names_to="component",values_to="mean_phi_mm_d") |>
  dplyr::mutate(component=sub("^phi_","",component))

p_annual <- ggplot2::ggplot(annual,
                            ggplot2::aes(year,mean_phi_mm_d,linetype=component)) +
  ggplot2::geom_hline(yintercept=0,linewidth=0.35) +
  ggplot2::geom_line(linewidth=0.45) +
  ggplot2::labs(x="Year",y="Mean annual Shapley allocation (mm/day)",
                linetype="PM input",
                title="Annual mean Shapley allocations to PM PET anomaly") +
  ggplot2::theme_bw(base_size=10) +
  ggplot2::theme(legend.position="bottom")

ggplot2::ggsave(file.path(OUT_DIR,"Fig_PET_PM_Shapley_annual.png"),p_annual,
                width=11,height=6.5,dpi=FIG_DPI)
ggplot2::ggsave(file.path(OUT_DIR,"Fig_PET_PM_Shapley_annual.pdf"),p_annual,
                width=11,height=6.5)

cat("\n15PET_shapley_decomposition_v2.R completed.\n")
cat("Outputs: ", OUT_DIR, "\n", sep="")
cat("Shapley closure max abs error: ",
    sprintf("%.3e", max(abs(shapley$shapley_closure_error_mm_d), na.rm=TRUE)),
    " mm/day\n", sep="")
