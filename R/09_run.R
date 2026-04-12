#' @title Run the Complete TwoTankModel Pipeline
#' @description Single entry-point function that runs the whole workflow
#'   from a config file.


#' Run the Full Model Pipeline from a Config File
#'
#' Reads a config.R file, loads data, calibrates the model, runs
#' uncertainty + sensitivity analyses, generates plots, and saves
#' everything to a dedicated results folder.
#'
#' This is the main user-facing function for running the complete
#' pipeline after installing the package.
#'
#' @param config_file Character. Path to config.R file. Default "config.R".
#'
#' @return Invisibly returns a list with all results:
#'   \code{cal}, \code{unc}, \code{sens}, \code{metrics}, \code{run_dir}.
#'
#' @examples
#' \dontrun{
#' # From within R
#' library(TwoTankModel)
#' run_model("config.R")
#'
#' # From the terminal
#' # Rscript -e 'TwoTankModel::run_model("config.R")'
#' }
#'
#' @export
run_model <- function(config_file = "config.R") {

  # ── Check config exists ──
  if (!file.exists(config_file)) {
    stop("\n  Config file not found: ", config_file,
         "\n  Create one using: create_config_template()\n")
  }

  # ── Source config into a new environment ──
  cfg <- new.env()
  source(config_file, local = cfg)

  # ── Validate required settings ──
  required <- c("RAINFALL_FILE", "DISCHARGE_FILE", "DISCHARGE_UNITS",
                "N_MC_SAMPLES", "K1_BOUNDS", "K2_BOUNDS", "K3_BOUNDS",
                "NSE_THRESHOLD", "OBJECTIVE", "RUN_NAME", "RESULTS_DIR")
  missing_settings <- setdiff(required, ls(cfg))
  if (length(missing_settings) > 0)
    stop("Missing config values: ", paste(missing_settings, collapse = ", "))

  # ── Defaults for optional settings ──
  if (!exists("CATCHMENT_AREA_KM2", envir = cfg))  cfg$CATCHMENT_AREA_KM2  <- NULL
  if (!exists("RUN_SENSITIVITY",    envir = cfg))  cfg$RUN_SENSITIVITY    <- TRUE
  if (!exists("SENSITIVITY_PERTURB",envir = cfg))  cfg$SENSITIVITY_PERTURB <- 0.20
  if (!exists("SAVE_PLOTS",         envir = cfg))  cfg$SAVE_PLOTS         <- TRUE
  if (!exists("SAVE_CSV_RESULTS",   envir = cfg))  cfg$SAVE_CSV_RESULTS   <- TRUE
  if (!exists("SAVE_PARAMETERS",    envir = cfg))  cfg$SAVE_PARAMETERS    <- TRUE
  if (!exists("SAVE_CONFIG_COPY",   envir = cfg))  cfg$SAVE_CONFIG_COPY   <- TRUE
  if (!exists("DATE_FORMAT",        envir = cfg))  cfg$DATE_FORMAT        <- "auto"
  if (!exists("RANDOM_SEED",        envir = cfg))  cfg$RANDOM_SEED        <- 42
  if (!exists("VERBOSE",            envir = cfg))  cfg$VERBOSE            <- TRUE

  # ── Banner ──
  cat("\n")
  cat("  ██████████████████████████████████████████████████████████████\n")
  cat(sprintf("  █    TWO-TANK MODEL — RUN: %-32s █\n", cfg$RUN_NAME))
  cat(sprintf("  █    Config file: %-42s █\n", config_file))
  cat("  ██████████████████████████████████████████████████████████████\n\n")

  # ── Step 1: Results folder ──
  cat("━━━ STEP 1: Creating results folder ━━━\n")
  run_dir  <- create_run_folder(cfg$RESULTS_DIR, cfg$RUN_NAME)
  plot_dir <- file.path(run_dir, "plots")
  csv_dir  <- file.path(run_dir, "csv")
  if (cfg$SAVE_CONFIG_COPY) {
    file.copy(config_file, file.path(run_dir, "config_used.R"), overwrite = TRUE)
    cat("  ✓ config_used.R saved\n\n")
  }

  # ── Step 2: Load data ──
  cat("━━━ STEP 2: Loading data ━━━\n")
  data <- load_data(
    rainfall_file   = cfg$RAINFALL_FILE,
    discharge_file  = cfg$DISCHARGE_FILE,
    discharge_units = cfg$DISCHARGE_UNITS,
    area_km2        = cfg$CATCHMENT_AREA_KM2,
    date_format     = cfg$DATE_FORMAT
  )

  # ── Step 3: Calibration ──
  cat("━━━ STEP 3: Calibration ━━━\n")
  cal <- calibrate_montecarlo(
    times      = data$times,
    precip_vec = data$precip,
    Q_obs      = data$Q_obs,
    n_samples  = cfg$N_MC_SAMPLES,
    k1_range   = cfg$K1_BOUNDS,
    k2_range   = cfg$K2_BOUNDS,
    k3_range   = cfg$K3_BOUNDS,
    area_km2   = cfg$CATCHMENT_AREA_KM2,
    objective  = cfg$OBJECTIVE,
    seed       = cfg$RANDOM_SEED,
    verbose    = cfg$VERBOSE
  )

  # ── Step 4: Uncertainty ──
  cat("━━━ STEP 4: Uncertainty analysis ━━━\n")
  unc <- extract_uncertainty(
    cal_result    = cal,
    times         = data$times,
    precip_vec    = data$precip,
    Q_obs         = data$Q_obs,
    nse_threshold = cfg$NSE_THRESHOLD,
    verbose       = cfg$VERBOSE
  )

  # ── Step 5: Sensitivity ──
  sens <- NULL
  if (cfg$RUN_SENSITIVITY) {
    cat("━━━ STEP 5: Sensitivity analysis ━━━\n")
    sens <- sensitivity_analysis(
      base_params = cal$best_params,
      times       = data$times,
      precip_vec  = data$precip,
      Q_obs       = data$Q_obs,
      area_km2    = cfg$CATCHMENT_AREA_KM2,
      perturb     = cfg$SENSITIVITY_PERTURB,
      verbose     = cfg$VERBOSE
    )
  }

  # ── Step 6: Performance ──
  cat("━━━ STEP 6: Performance metrics ━━━\n")
  Q_sim_cmp <- if (!is.null(cfg$CATCHMENT_AREA_KM2)) cal$best_sim$Q_total_m3s
               else cal$best_sim$Q_total
  metrics <- calc_all_metrics(data$Q_obs, Q_sim_cmp)
  units_str <- if (!is.null(cfg$CATCHMENT_AREA_KM2)) "m3s" else "mm"
  print_model_summary(cal$best_sim, label = "Calibrated Model", units = units_str)

  # ── Step 7: Hydrograph metrics & monthly balance ──
  cat("━━━ STEP 7: Hydrograph & monthly analysis ━━━\n")
  hydro_met <- extract_metrics(cal$best_sim, data$dates,
                               label = "Calibrated Model", verbose = TRUE)
  mon <- monthly_summary(cal$best_sim, data$dates)
  if (cfg$SAVE_CSV_RESULTS)
    write.csv(mon, file.path(csv_dir, "monthly_summary.csv"), row.names = FALSE)

  # ── Step 8: Mass balance ──
  cat("━━━ STEP 8: Mass balance check ━━━\n")
  mb <- check_mass_balance(cal$best_sim)

  # ── Step 9: Plots ──
  if (cfg$SAVE_PLOTS) {
    cat("━━━ STEP 9: Generating plots ━━━\n")
    plots <- plot_all(data$dates, data$precip, data$Q_obs, cal, unc,
                      output_dir = plot_dir, prefix = cfg$RUN_NAME)
  }

  # ── Step 10: Save ──
  cat("━━━ STEP 10: Saving results ━━━\n")
  if (cfg$SAVE_CSV_RESULTS) {
    export_results(cal, unc, sens,
                   dates = data$dates, Q_obs = data$Q_obs,
                   precip = data$precip,
                   output_dir = csv_dir, prefix = cfg$RUN_NAME)
  }
  if (cfg$SAVE_PARAMETERS)
    save_parameters(cal, file.path(run_dir, "parameters.txt"))

  # ── Done ──
  cat("\n")
  cat("  ██████████████████████████████████████████████████████████████\n")
  cat("  █                                                            █\n")
  cat("  █    ✓ RUN COMPLETE                                          █\n")
  cat(sprintf("  █    Results: %-46s █\n", run_dir))
  cat("  █                                                            █\n")
  cat("  ██████████████████████████████████████████████████████████████\n\n")

  invisible(list(cal = cal, unc = unc, sens = sens,
                 metrics = metrics, run_dir = run_dir))
}


#' Create a Config Template in the Current Directory
#'
#' Copies a default \code{config.R} template to the working directory
#' so you can edit it and run the model.
#'
#' @param path Character. Where to save the template. Default "config.R".
#' @param overwrite Logical. Default FALSE.
#'
#' @export
create_config_template <- function(path = "config.R", overwrite = FALSE) {
  if (file.exists(path) && !overwrite)
    stop(path, " already exists. Use overwrite = TRUE to replace it.")

  # Try to find shipped template first
  shipped <- system.file("config_template.R", package = "TwoTankModel")
  if (nzchar(shipped) && file.exists(shipped)) {
    file.copy(shipped, path, overwrite = overwrite)
  } else {
    # Fall back to embedded template
    writeLines(config_template_text(), path)
  }
  cat("  ✓ Config template created: ", path, "\n", sep = "")
  cat("  Edit it, then run: TwoTankModel::run_model('", path, "')\n", sep = "")
}


# ── Embedded config template (used if installed package copy not found) ─────
config_template_text <- function() {
'################################################################################
#  CONFIG.R — TwoTankModel Configuration
#  Edit this file, then run:  Rscript -e \'TwoTankModel::run_model("config.R")\'
################################################################################

# ── Input data ──
RAINFALL_FILE      <- "data/rainfall.csv"      # columns: date, P
DISCHARGE_FILE     <- "data/discharge.csv"     # columns: date, Q
DISCHARGE_UNITS    <- "m3s"                    # "m3s" or "mm_day"
DATE_FORMAT        <- "auto"                   # example "%m/%d/%Y" 
CATCHMENT_AREA_KM2 <- 150                      # required if m3s

# ── Calibration ──
N_MC_SAMPLES       <- 5000                     # 1000 quick, 5000 standard
K1_BOUNDS          <- c(0.01, 0.80)            # surface runoff coef
K2_BOUNDS          <- c(0.01, 0.50)            # percolation coef
K3_BOUNDS          <- c(0.001, 0.15)           # baseflow coef
NSE_THRESHOLD      <- 0.5                      # for behavioural sets
OBJECTIVE          <- "NSE"                    # "NSE", "KGE", "LogNSE"

# ── Sensitivity ──
RUN_SENSITIVITY    <- TRUE
SENSITIVITY_PERTURB <- 0.20                    # ±20%

# ── Output ──
RUN_NAME           <- "run_01"
RESULTS_DIR        <- "results"
SAVE_PLOTS         <- TRUE
SAVE_CSV_RESULTS   <- TRUE
SAVE_PARAMETERS    <- TRUE
SAVE_CONFIG_COPY   <- TRUE

# ── Advanced ──
RANDOM_SEED        <- 42
VERBOSE            <- TRUE
'
}
