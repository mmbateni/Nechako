###############################################################################
# run_all.R
#
# MASTER RUNNER — Nechako Drought Analysis Framework
#
# Runs the complete analysis pipeline in dependency order.
#
# IMPORTANT:
#   Each analysis script is executed in a NEW R process using Rscript.
#   This isolates memory, variables, parallel workers, and package state
#   between stages.
###############################################################################

# =============================================================================
# CONFIGURATION
# =============================================================================

PROJECT_DIR <- "D:/Nechako_Drought/Nechako"

if (!dir.exists(PROJECT_DIR)) {
  stop("Project directory does not exist: ", PROJECT_DIR)
}

setwd(PROJECT_DIR)

# Stop the complete framework immediately if any stage fails.
STOP_ON_ERROR <- TRUE

# Optional stages can be disabled here if needed.
RUN_PR_PET_TRENDS         <- TRUE
RUN_PR_PET_VISUALIZATION  <- TRUE
RUN_DROUGHT_TRENDS        <- TRUE
RUN_BASIN_TIMESERIES      <- TRUE
RUN_TREND_VISUALIZATION   <- TRUE
RUN_EVENT_RANKING         <- TRUE
RUN_DECOMPOSITION         <- TRUE
RUN_PET_BIAS              <- TRUE
RUN_MANUSCRIPT            <- TRUE


# =============================================================================
# PIPELINE
# =============================================================================

pipeline <- data.frame(
  stage = c(
    "01",
    "02a",
    "02b",
    "03",
    "04a",
    "04b",
    "06",
    "07",
    "08",
    "09",
    "10a",
    "10b",
    "11"
  ),
  
  description = c(
    "SPI calculation",
    "Temperature detrending",
    "PET calculation",
    "SPEI calculation",
    "Precipitation/PET/temperature trend analysis",
    "Precipitation/PET trend visualization",
    "Spatial drought-index trend analysis",
    "Basin-averaged time-series analysis",
    "Spatial drought trend visualization",
    "Drought event ranking",
    "Dynamic/thermodynamic decomposition",
    "PET bias non-stationarity analysis",
    "Manuscript tables, figures and statistics"
  ),
  
  script = c(
    "1SPI_ERALand.R",
    "2a_detrend_temperature.R",
    "2b_PET_ERALand.R",
    "3SPEI_ERALand.R",
    "4pr_pet_trends.r",
    "4pr_pet_trends_visualization.r",
    "6trend_test_ALL.R",
    "7basin_timeseries.R",
    "8trends_visualization.R",
    "9event_ranking.R",
    "10a_dynamic_thermodynamic_decomp.R",
    "10b_PET_bias_nonstationarity.R",
    "11nechako_drought_manuscript1.R"
  ),
  
  run = c(
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    RUN_PR_PET_TRENDS,
    RUN_PR_PET_VISUALIZATION,
    RUN_DROUGHT_TRENDS,
    RUN_BASIN_TIMESERIES,
    RUN_TREND_VISUALIZATION,
    RUN_EVENT_RANKING,
    RUN_DECOMPOSITION,
    RUN_PET_BIAS,
    RUN_MANUSCRIPT
  ),
  
  stringsAsFactors = FALSE
)


# =============================================================================
# PRE-FLIGHT CHECKS
# =============================================================================

cat("\n")
cat("============================================================\n")
cat(" NECHAKO DROUGHT ANALYSIS — MASTER FRAMEWORK\n")
cat("============================================================\n")
cat("Project directory:\n  ", normalizePath(PROJECT_DIR), "\n\n", sep = "")

# Find Rscript executable from the currently running R installation.
rscript <- file.path(R.home("bin"), "Rscript.exe")

if (!file.exists(rscript)) {
  # Linux/macOS fallback
  rscript <- file.path(R.home("bin"), "Rscript")
}

if (!file.exists(rscript)) {
  stop("Could not locate Rscript executable.")
}

cat("Rscript:\n  ", rscript, "\n\n", sep = "")


# Check required scripts.
required_scripts <- pipeline$script[pipeline$run]

missing_scripts <- required_scripts[
  !file.exists(file.path(PROJECT_DIR, required_scripts))
]

if (length(missing_scripts) > 0L) {
  stop(
    "The following required scripts are missing:\n  ",
    paste(missing_scripts, collapse = "\n  ")
  )
}


# Several downstream scripts require the shared utility file.
utils_file <- file.path(PROJECT_DIR, "DROUGHT_ANALYSIS_utils.R")

downstream_enabled <- any(
  pipeline$run[
    pipeline$script %in% c(
      "4pr_pet_trends.r",
      "4pr_pet_trends_visualization.r",
      "6trend_test_ALL.R",
      "7basin_timeseries.R",
      "8trends_visualization.R",
      "9event_ranking.R",
      "10a_dynamic_thermodynamic_decomp.R"
    )
  ]
)

if (downstream_enabled && !file.exists(utils_file)) {
  stop(
    "Required shared utility file is missing:\n  ",
    utils_file
  )
}


# =============================================================================
# LOGGING
# =============================================================================

LOG_DIR <- file.path(PROJECT_DIR, "framework_logs")

if (!dir.exists(LOG_DIR)) {
  dir.create(LOG_DIR, recursive = TRUE)
}

run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")

MASTER_LOG <- file.path(
  LOG_DIR,
  paste0("framework_", run_id, ".log")
)

cat(
  "Nechako Drought Analysis Framework\n",
  "Run ID: ", run_id, "\n",
  "Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
  "Project: ", PROJECT_DIR, "\n",
  "R version: ", R.version.string, "\n",
  "============================================================\n\n",
  file = MASTER_LOG,
  sep = ""
)


log_master <- function(...) {
  
  txt <- paste0(...)
  
  cat(txt, "\n")
  
  cat(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " | ",
    txt,
    "\n",
    file = MASTER_LOG,
    append = TRUE,
    sep = ""
  )
  
  invisible(NULL)
}


# =============================================================================
# SCRIPT EXECUTION FUNCTION
# =============================================================================

run_script <- function(stage, description, script) {
  
  script_path <- file.path(PROJECT_DIR, script)
  
  stage_log <- file.path(
    LOG_DIR,
    paste0(
      run_id,
      "_",
      stage,
      "_",
      tools::file_path_sans_ext(basename(script)),
      ".log"
    )
  )
  
  cat("\n")
  cat("============================================================\n")
  cat(" STAGE ", stage, ": ", description, "\n", sep = "")
  cat(" Script: ", script, "\n", sep = "")
  cat("============================================================\n\n")
  
  start_time <- Sys.time()
  
  log_master(
    "START | Stage ", stage,
    " | ", description,
    " | ", script
  )
  
  # Run each analysis script in an independent R session.
  #
  # stdout and stderr are written to the same per-stage log.
  # This avoids filling the RStudio console with potentially millions
  # of lines while preserving the complete diagnostic output.
  
  exit_status <- tryCatch(
    {
      system2(
        command = rscript,
        args = c(
          "--vanilla",
          shQuote(script_path)
        ),
        stdout = stage_log,
        stderr = stage_log,
        wait = TRUE
      )
    },
    error = function(e) {
      attr(e, "runner_error") <- TRUE
      e
    }
  )
  
  elapsed <- as.numeric(
    difftime(Sys.time(), start_time, units = "mins")
  )
  
  # system2 itself failed.
  if (inherits(exit_status, "error")) {
    
    log_master(
      "FAILED | Stage ", stage,
      " | Runner error: ",
      conditionMessage(exit_status)
    )
    
    cat("\nStage log:\n", stage_log, "\n", sep = "")
    
    return(FALSE)
  }
  
  # Rscript returned a non-zero status.
  if (!identical(as.integer(exit_status), 0L)) {
    
    log_master(
      "FAILED | Stage ", stage,
      " | Exit status ", exit_status,
      " | Elapsed ",
      sprintf("%.2f", elapsed),
      " min"
    )
    
    cat("\nStage failed.\n")
    cat("Diagnostic log:\n  ", stage_log, "\n", sep = "")
    
    # Print the last part of the failed log directly to the console.
    if (file.exists(stage_log)) {
      
      log_lines <- readLines(
        stage_log,
        warn = FALSE,
        encoding = "UTF-8"
      )
      
      n <- length(log_lines)
      
      if (n > 0L) {
        from <- max(1L, n - 39L)
        
        cat("\n")
        cat("--------------- LAST LOG LINES ----------------\n")
        cat(paste(log_lines[from:n], collapse = "\n"))
        cat("\n-------------------------------------------------\n")
      }
    }
    
    return(FALSE)
  }
  
  log_master(
    "SUCCESS | Stage ", stage,
    " | Elapsed ",
    sprintf("%.2f", elapsed),
    " min"
  )
  
  cat(
    "\nCompleted in ",
    sprintf("%.2f", elapsed),
    " minutes.\n",
    sep = ""
  )
  
  cat(
    "Stage log: ",
    stage_log,
    "\n",
    sep = ""
  )
  
  TRUE
}


# =============================================================================
# RUN COMPLETE FRAMEWORK
# =============================================================================

framework_start <- Sys.time()

results <- data.frame(
  stage       = character(),
  script      = character(),
  status      = character(),
  stringsAsFactors = FALSE
)


for (i in seq_len(nrow(pipeline))) {
  
  p <- pipeline[i, ]
  
  if (!isTRUE(p$run)) {
    
    log_master(
      "SKIPPED | Stage ",
      p$stage,
      " | ",
      p$script
    )
    
    results <- rbind(
      results,
      data.frame(
        stage  = p$stage,
        script = p$script,
        status = "SKIPPED",
        stringsAsFactors = FALSE
      )
    )
    
    next
  }
  
  ok <- run_script(
    stage       = p$stage,
    description = p$description,
    script      = p$script
  )
  
  results <- rbind(
    results,
    data.frame(
      stage  = p$stage,
      script = p$script,
      status = if (ok) "SUCCESS" else "FAILED",
      stringsAsFactors = FALSE
    )
  )
  
  # Release memory held by the master process itself.
  gc(verbose = FALSE)
  
  if (!ok && STOP_ON_ERROR) {
    
    log_master(
      "FRAMEWORK ABORTED after failure in ",
      p$script
    )
    
    break
  }
}


# =============================================================================
# FINAL SUMMARY
# =============================================================================

framework_end <- Sys.time()

elapsed_hours <- as.numeric(
  difftime(framework_end, framework_start, units = "hours")
)

summary_file <- file.path(
  LOG_DIR,
  paste0("framework_", run_id, "_summary.csv")
)

write.csv(
  results,
  summary_file,
  row.names = FALSE
)

cat("\n\n")
cat("============================================================\n")
cat(" FRAMEWORK EXECUTION SUMMARY\n")
cat("============================================================\n\n")

print(results, row.names = FALSE)

cat("\nTotal elapsed time: ",
    sprintf("%.2f", elapsed_hours),
    " hours\n",
    sep = "")

cat("Master log:\n  ", MASTER_LOG, "\n", sep = "")
cat("Summary:\n  ", summary_file, "\n", sep = "")

if (any(results$status == "FAILED")) {
  
  failed <- results$script[results$status == "FAILED"]
  
  cat("\nFRAMEWORK FAILED\n")
  cat(
    "Failed stage(s): ",
    paste(failed, collapse = ", "),
    "\n",
    sep = ""
  )
  
  stop(
    "Framework terminated because one or more stages failed.",
    call. = FALSE
  )
  
} else {
  
  log_master(
    "FRAMEWORK COMPLETED SUCCESSFULLY | Total elapsed ",
    sprintf("%.2f", elapsed_hours),
    " hours"
  )
  
  cat("\n============================================================\n")
  cat(" FRAMEWORK COMPLETED SUCCESSFULLY\n")
  cat("============================================================\n")
}