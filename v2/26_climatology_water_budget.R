
# 16_climatology_water_budget.R
# REGIONAL CLIMATOLOGY AND WATER-BUDGET ANALYSIS
#
# Regions: whole Nechako basin (dissolved union of 7 supplied subbasins) +
#          the 11 individual subbasins.
#
# Spatial means are authoritative outputs from 4pr_pet_trends.r and use
# fractional polygon coverage on the native ERA5-Land analysis grid.
#
# AET is now a required component of the water-budget framework. Script 4 must
# supply aet_mm_month in regional_avg_monthly, derived from ERA5-Land
# total_evapotranspiration_monthly.nc (total evaporation / AET). Script 16 consumes
# that regional AET field; it does not independently read or spatially
# aggregate the raw NetCDF. PET is never substituted for AET.
################################################################################

setwd(Sys.getenv("NECHAKO_WD", "D:/Nechako_Drought/Nechako/"))
source("DROUGHT_ANALYSIS_utils.R")
utils_load_packages(c("data.table", "dplyr", "lubridate", "Kendall", "trend", "ggplot2"))

if (!dir.exists(WD_PATH)) stop("Working directory not found: ", WD_PATH)
setwd(WD_PATH)

IN_DIR <- file.path(WD_PATH, "trend_analysis_pr_pet")
OUT_DIR <- file.path(WD_PATH, "temporal_drought", "climatology_water_budget")
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

META_FILE <- file.path(IN_DIR, "analysis_metadata.rds")
if (!file.exists(META_FILE)) stop("Missing ", META_FILE, ". Run 4pr_pet_trends.r first.")
metadata <- readRDS(META_FILE)
if (is.null(metadata$regional_avg_monthly))
  stop("analysis_metadata.rds does not contain regional_avg_monthly. Re-run 4pr_pet_trends.r.")

d <- data.table::as.data.table(metadata$regional_avg_monthly)
setnames(d, trimws(names(d)))
d[, date := as.Date(date)]
setorder(d, region_id, date)
d[, year := as.integer(format(date,"%Y"))]
d[, month := as.integer(format(date,"%m"))]
d[, season := fifelse(month %in% c(12L,1L,2L),"DJF",
                      fifelse(month %in% 3:5,"MAM",fifelse(month %in% 6:8,"JJA","SON")))]

req <- c(
  "region_id", "region_name", "date", "precip_mm_month",
  "pet_mm_month", "pet_era5_mm_month", "tair_degC_month",
  "aet_mm_month"
)
miss <- setdiff(req, names(d))
if(length(miss)) {
  stop(
    "Regional climate data missing required field(s): ",
    paste(miss, collapse = ", "),
    "\nRe-run 4pr_pet_trends.r after the AET update. Script 4 must provide ",
    "aet_mm_month from ERA5-Land total evaporation (total_evapotranspiration_monthly.nc)."
  )
}

# AET is authoritative here and must be finite for the complete regional
# monthly record. Do not silently replace missing AET with PET or NA values.
aet_missing <- !is.finite(d$aet_mm_month)
if(any(aet_missing)) {
  bad_regions <- unique(d$region_name[aet_missing])
  stop(
    "AET (aet_mm_month) contains ", sum(aet_missing),
    " non-finite regional month(s). ",
    "The revised water-budget framework requires complete AET coverage ",
    "before climatology/anomaly/trend products are generated. Affected region(s): ",
    paste(bad_regions, collapse = ", "),
    ". Check 4pr_pet_trends.r and total_evapotranspiration_monthly.nc."
  )
}

has_aet <- TRUE


# Monthly water-budget diagnostics. Ratios are defined only when denominators
# are positive and finite.
d[, `:=`(
  AET_over_PET = fifelse(is.finite(aet_mm_month)&is.finite(pet_mm_month)&pet_mm_month>0,aet_mm_month/pet_mm_month,NA_real_),
  AET_over_P = fifelse(is.finite(aet_mm_month)&is.finite(precip_mm_month)&precip_mm_month>0,aet_mm_month/precip_mm_month,NA_real_),
  PET_over_P = fifelse(is.finite(pet_mm_month)&is.finite(precip_mm_month)&precip_mm_month>0,pet_mm_month/precip_mm_month,NA_real_),
  residual_runoff_over_P = fifelse(is.finite(aet_mm_month)&is.finite(precip_mm_month)&precip_mm_month>0,(precip_mm_month-aet_mm_month)/precip_mm_month,NA_real_),
  residual_water_balance_mm = ifelse(is.finite(aet_mm_month),precip_mm_month-aet_mm_month,NA_real_),
  atmospheric_demand_deficit_mm = ifelse(is.finite(aet_mm_month),pet_mm_month-aet_mm_month,NA_real_)
)]
d[, PET_ERA5_over_P := fifelse(is.finite(pet_era5_mm_month)&is.finite(precip_mm_month)&precip_mm_month>0,pet_era5_mm_month/precip_mm_month,NA_real_)]

period_summary <- function(x) {
  P <- sum(x$precip_mm_month,na.rm=TRUE)
  PM <- sum(x$pet_mm_month,na.rm=TRUE)
  E5 <- if(any(is.finite(x$pet_era5_mm_month))) sum(x$pet_era5_mm_month,na.rm=TRUE) else NA_real_
  A <- if(any(is.finite(x$aet_mm_month))) sum(x$aet_mm_month,na.rm=TRUE) else NA_real_
  out <- data.table(
    precipitation_mm=P,AET_mm=A,PET_PM_mm=PM,
    PET_ERA5_mm=if(is.finite(E5)) E5 else NA_real_,
    temperature_degC=if(any(is.finite(x$tair_degC_month))) mean(x$tair_degC_month,na.rm=TRUE) else NA_real_,
    P_minus_AET_mm=if(is.finite(A)) P-A else NA_real_,
    P_minus_PET_PM_mm=P-PM,
    AET_over_P=if(is.finite(A)&&P>0) A/P else NA_real_,
    PET_PM_over_P=if(P>0) PM/P else NA_real_,
    AET_over_PET_PM=if(is.finite(A)&&PM>0) A/PM else NA_real_,
    PET_ERA5_over_P=if(is.finite(E5)&&P>0) E5/P else NA_real_,
    residual_runoff_over_P=if(is.finite(A)&&P>0) (P-A)/P else NA_real_
  )
  out
}

# 1. Monthly climatology
monthly_clim <- d[, {
  P <- mean(precip_mm_month,na.rm=TRUE); PM <- mean(pet_mm_month,na.rm=TRUE); E5 <- mean(pet_era5_mm_month,na.rm=TRUE); A <- mean(aet_mm_month,na.rm=TRUE)
  list(month_name=month.name[unique(month)],precipitation_mm=P,AET_mm=if(is.finite(A))A else NA_real_,PET_PM_mm=PM,PET_ERA5_mm=if(is.finite(E5))E5 else NA_real_,temperature_degC=mean(tair_degC_month,na.rm=TRUE),
       P_minus_AET_mm=if(is.finite(A))P-A else NA_real_,P_minus_PET_PM_mm=P-PM,AET_over_P=if(is.finite(A)&&P>0)A/P else NA_real_,PET_PM_over_P=if(P>0)PM/P else NA_real_,AET_over_PET_PM=if(is.finite(A)&&PM>0)A/PM else NA_real_,PET_ERA5_over_P=if(is.finite(E5)&&P>0)E5/P else NA_real_,residual_runoff_over_P=if(is.finite(A)&&P>0)(P-A)/P else NA_real_,n_years=data.table::uniqueN(year))
},by=.(region_id,region_name,month)]
setorder(monthly_clim,region_id,month)
fwrite(monthly_clim,file.path(OUT_DIR,"climatology_monthly_1950_2025.csv"))

# 2. Seasonal climatology: seasonal totals are 3 × the mean monthly depth.
season_order <- c("DJF","MAM","JJA","SON")
season_clim <- d[, {
  P<-mean(precip_mm_month,na.rm=TRUE)*3; PM<-mean(pet_mm_month,na.rm=TRUE)*3; E5<-mean(pet_era5_mm_month,na.rm=TRUE)*3; A<-mean(aet_mm_month,na.rm=TRUE)*3
  list(precipitation_mm=P,AET_mm=if(is.finite(A))A else NA_real_,PET_PM_mm=PM,PET_ERA5_mm=if(is.finite(E5))E5 else NA_real_,temperature_degC=mean(tair_degC_month,na.rm=TRUE),P_minus_AET_mm=if(is.finite(A))P-A else NA_real_,P_minus_PET_PM_mm=P-PM,AET_over_P=if(is.finite(A)&&P>0)A/P else NA_real_,PET_PM_over_P=if(P>0)PM/P else NA_real_,AET_over_PET_PM=if(is.finite(A)&&PM>0)A/PM else NA_real_,PET_ERA5_over_P=if(is.finite(E5)&&P>0)E5/P else NA_real_,residual_runoff_over_P=if(is.finite(A)&&P>0)(P-A)/P else NA_real_,n_months=.N)
},by=.(region_id,region_name,season)]
season_clim[,season:=factor(season,levels=season_order)]; setorder(season_clim,region_id,season); season_clim[,season:=as.character(season)]
fwrite(season_clim,file.path(OUT_DIR,"climatology_seasonal_1950_2025.csv"))

# 3. Annual full-year water budget
annual <- d[, {
  P<-sum(precip_mm_month,na.rm=TRUE); PM<-sum(pet_mm_month,na.rm=TRUE); E5<-if(any(is.finite(pet_era5_mm_month)))sum(pet_era5_mm_month,na.rm=TRUE) else NA_real_; A<-if(any(is.finite(aet_mm_month)))sum(aet_mm_month,na.rm=TRUE) else NA_real_
  list(precipitation_mm=P,AET_mm=A,PET_PM_mm=PM,PET_ERA5_mm=if(is.finite(E5))E5 else NA_real_,temperature_degC=mean(tair_degC_month,na.rm=TRUE),P_minus_AET_mm=if(is.finite(A))P-A else NA_real_,P_minus_PET_PM_mm=P-PM,AET_over_P=if(is.finite(A)&&P>0)A/P else NA_real_,PET_PM_over_P=if(P>0)PM/P else NA_real_,AET_over_PET_PM=if(is.finite(A)&&PM>0)A/PM else NA_real_,PET_ERA5_over_P=if(is.finite(E5)&&P>0)E5/P else NA_real_,residual_runoff_over_P=if(is.finite(A)&&P>0)(P-A)/P else NA_real_,n_months=sum(is.finite(precip_mm_month)))
},by=.(region_id,region_name,year)]
setorder(annual,region_id,year)
fwrite(annual,file.path(OUT_DIR,"annual_water_budget_1950_2025.csv"))
annual_clim <- annual[,lapply(.SD,mean,na.rm=TRUE),by=.(region_id,region_name),.SDcols=setdiff(names(annual),c("region_id","region_name","year"))]
fwrite(annual_clim,file.path(OUT_DIR,"climatology_annual_1950_2025.csv"))

# 4. Recent 2022–2025 diagnostics and climatological anomalies.
recent <- d[year>=2022 & year<=2025]
fwrite(recent,file.path(OUT_DIR,"recent_2022_2025_monthly_climate_water_budget.csv"))
recent_seasonal <- recent[,period_summary(.SD),by=.(region_id,region_name,season)]
fwrite(recent_seasonal,file.path(OUT_DIR,"recent_2022_2025_seasonal_climate_water_budget.csv"))
recent_annual <- recent[,period_summary(.SD),by=.(region_id,region_name,year)]
fwrite(recent_annual,file.path(OUT_DIR,"recent_2022_2025_annual_climate_water_budget.csv"))

anom <- d[,.(
  region_id, region_name, date, year, month,
  precip_mm_month, pet_mm_month, pet_era5_mm_month,
  aet_mm_month, tair_degC_month
)]

clim <- monthly_clim[,.(
  region_id, month,
  precip_clim = precipitation_mm,
  pet_clim = PET_PM_mm,
  pet_era5_clim = PET_ERA5_mm,
  aet_clim = AET_mm,
  tair_clim = temperature_degC
)]

anom <- merge(
  anom[year >= 2022 & year <= 2025],
  clim,
  by = c("region_id", "month"),
  all.x = TRUE,
  sort = TRUE
)

# One anomaly is retained per regional calendar month. AET is a direct
# difference from the corresponding 1950-2025 calendar-month climatology,
# expressed in mm/month.
anom[,`:=`(
  precip_anomaly_mm = precip_mm_month - precip_clim,
  pet_PM_anomaly_mm = pet_mm_month - pet_clim,
  pet_ERA5_anomaly_mm = pet_era5_mm_month - pet_era5_clim,
  AET_anomaly_mm = aet_mm_month - aet_clim,
  tair_anomaly_degC = tair_degC_month - tair_clim
)]

fwrite(
  anom,
  file.path(OUT_DIR, "recent_2022_2025_monthly_climatological_anomalies.csv")
)

# 5. Annual trend tests for water-budget ratios. Mann–Kendall and Sen slope
# are calculated independently for every region and variable.
trend_vars <- c("PET_PM_over_P","PET_ERA5_over_P","AET_over_P","AET_over_PET_PM","residual_runoff_over_P")
trend_tbl <- rbindlist(lapply(trend_vars,function(v){
  rbindlist(lapply(split(annual,annual$region_id),function(z){
    x<-z[[v]]; ok<-is.finite(x); x<-x[ok]
    mk<-if(length(x)>=10)tryCatch(Kendall::MannKendall(x),error=function(e)NULL) else NULL
    ss<-if(length(x)>=10)tryCatch(trend::sens.slope(x),error=function(e)NULL) else NULL
    data.table(region_id=unique(z$region_id),region_name=unique(z$region_name),variable=v,n_years=length(x),tau=if(!is.null(mk))as.numeric(mk$tau)else NA_real_,p_value=if(!is.null(mk))as.numeric(mk$sl)else NA_real_,sen_slope_per_year=if(!is.null(ss))as.numeric(ss$estimates)else NA_real_)
  }),fill=TRUE)
}),fill=TRUE)
trend_tbl[,p_fdr:=if(sum(is.finite(p_value))>1)p.adjust(p_value,"BH")else p_value,by=.(variable)]
fwrite(trend_tbl,file.path(OUT_DIR,"annual_water_budget_ratio_trends_1950_2025.csv"))

# 6. Practitioner-facing table: annual P, PM PET and diagnostic residuals.
practitioner <- annual[,.(region_id,region_name,year,precipitation_mm,PET_PM_mm,PET_ERA5_mm,AET_mm,P_minus_AET_mm,P_minus_PET_PM_mm,AET_over_P,PET_PM_over_P,AET_over_PET_PM,residual_runoff_over_P)]
fwrite(practitioner,file.path(OUT_DIR,"practitioner_annual_water_budget_1950_2025.csv"))

# 7. Regional figures.
p1 <- ggplot2::ggplot(annual,ggplot2::aes(year,PET_PM_mm)) + ggplot2::geom_line(linewidth=.35) +
  ggplot2::facet_wrap(~region_name,ncol=3,scales="free_y") + ggplot2::labs(title="Annual full-year PM PET",x="Year",y="PET-PM (mm/year)") + ggplot2::theme_bw(base_size=9)
ggplot2::ggsave(file.path(OUT_DIR,"regional_annual_PM_PET.png"),p1,width=14,height=11,dpi=300)
p2 <- ggplot2::ggplot(anom,ggplot2::aes(date,precip_anomaly_mm)) + ggplot2::geom_hline(yintercept=0,linetype="dashed",linewidth=.3) + ggplot2::geom_line(linewidth=.3) + ggplot2::facet_wrap(~region_name,ncol=3) + ggplot2::labs(title="2022–2025 precipitation anomalies by Nechako region",x=NULL,y="P anomaly (mm/month)") + ggplot2::theme_bw(base_size=9)
ggplot2::ggsave(file.path(OUT_DIR,"regional_recent_precipitation_anomalies.png"),p2,width=14,height=11,dpi=300)

cat("==============================================================\n")
cat("Regional climatology / water-budget analysis complete.\n")
cat("Regions: ",data.table::uniqueN(d$region_id),"\n",sep="")
cat("AET available: ",has_aet,"\n",sep="")
cat("AET source: ERA5-Land total evaporation (total_evapotranspiration_monthly.nc), aggregated by Script 4\n")
cat("Output directory: ",OUT_DIR,"\n",sep="")
cat("==============================================================\n")
