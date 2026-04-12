################################################################################
#  MAIN.R — MAIN ENTRY POINT
#  ─────────────────────────
#  This script reads config.R and runs the complete analysis:
#    1. Loads data (rainfall + observed discharge) from files in config
#    2. Calibrates the two-tank model
#    3. Runs uncertainty analysis
#    4. Runs sensitivity analysis
#    5. Saves all results to a dedicated folder
#
#  USAGE:
#    1. Edit config.R with your data paths and settings
#    2. From terminal:  Rscript main.R
#       Or in R:        source("main.R")
################################################################################

# ── Set CRAN mirror for terminal use ─────────────────────────────────────────
options(repos = c(CRAN = "https://cloud.r-project.org"))

# ── Auto-install required packages ───────────────────────────────────────────
required_pkgs <- c("deSolve", "lhs", "ggplot2", "gridExtra")
for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ── Load configuration ───────────────────────────────────────────────────────
if (!file.exists("config.R")) {
  stop("config.R not found. Copy config.R from the package and edit it.")
}
source("config.R")

# ── Load package functions ───────────────────────────────────────────────────
# (If installed as a package, replace this loop with: library(TwoTankModel))
for (f in list.files("R", full.names = TRUE, pattern = "\\.R$")) source(f)


# ═══════════════════════════════════════════════════════════════════════════════

cat("\n")
cat("  ██████████████████████████████████████████████████████████████\n")
cat("  █                                                            █\n")
cat("  █    TWO-TANK MODEL — RUN:  ", sprintf("%-25s", RUN_NAME), "█\n")
cat("  █                                                            █\n")
cat("  ██████████████████████████████████████████████████████████████\n\n")


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1: CREATE RESULTS FOLDER
# ═══════════════════════════════════════════════════════════════════════════════

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  STEP 1: Setting up results folder\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

run_dir  <- create_run_folder(RESULTS_DIR, RUN_NAME)
plot_dir <- file.path(run_dir, "plots")
csv_dir  <- file.path(run_dir, "csv")

# Save a copy of config.R for reproducibility
if (SAVE_CONFIG_COPY) {
  file.copy("config.R", file.path(run_dir, "config_used.R"), overwrite = TRUE)
  cat("  ✓ config_used.R saved\n\n")
}


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2: LOAD DATA
# ═══════════════════════════════════════════════════════════════════════════════

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  STEP 2: Loading data\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

data <- load_data(
  rainfall_file   = RAINFALL_FILE,
  discharge_file  = DISCHARGE_FILE,
  discharge_units = DISCHARGE_UNITS,
  area_km2        = CATCHMENT_AREA_KM2
)


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3: CALIBRATE
# ═══════════════════════════════════════════════════════════════════════════════

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  STEP 3: Calibration\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

cal <- calibrate_montecarlo(
  times      = data$times,
  precip_vec = data$precip,
  Q_obs      = data$Q_obs,
  n_samples  = N_MC_SAMPLES,
  k1_range   = K1_BOUNDS,
  k2_range   = K2_BOUNDS,
  k3_range   = K3_BOUNDS,
  area_km2   = CATCHMENT_AREA_KM2,
  objective  = OBJECTIVE,
  seed       = RANDOM_SEED,
  verbose    = VERBOSE
)


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4: UNCERTAINTY
# ═══════════════════════════════════════════════════════════════════════════════

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  STEP 4: Uncertainty analysis\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

unc <- extract_uncertainty(
  cal_result    = cal,
  times         = data$times,
  precip_vec    = data$precip,
  Q_obs         = data$Q_obs,
  nse_threshold = NSE_THRESHOLD,
  verbose       = VERBOSE
)


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 5: SENSITIVITY ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════════

sens <- NULL
if (RUN_SENSITIVITY) {
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("  STEP 5: Sensitivity analysis\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

  sens <- sensitivity_analysis(
    base_params = cal$best_params,
    times       = data$times,
    precip_vec  = data$precip,
    Q_obs       = data$Q_obs,
    area_km2    = CATCHMENT_AREA_KM2,
    perturb     = SENSITIVITY_PERTURB,
    verbose     = VERBOSE
  )
}


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 6: PERFORMANCE & HYDROGRAPH METRICS
# ═══════════════════════════════════════════════════════════════════════════════

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  STEP 6: Performance & hydrograph metrics\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

# Compare in the appropriate units
Q_sim_cmp <- if (!is.null(CATCHMENT_AREA_KM2)) cal$best_sim$Q_total_m3s
             else cal$best_sim$Q_total
metrics <- calc_all_metrics(data$Q_obs, Q_sim_cmp)

# Hydrograph metrics in m³/s (or mm) as appropriate
units_str <- if (!is.null(CATCHMENT_AREA_KM2)) "m3s" else "mm"
print_model_summary(cal$best_sim, label = "Calibrated Model", units = units_str)
hydro_met <- extract_metrics(cal$best_sim, data$dates, label = "Calibrated Model",
                             verbose = TRUE)


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 7: MONTHLY SUMMARY
# ═══════════════════════════════════════════════════════════════════════════════

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  STEP 7: Monthly water balance\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

mon <- monthly_summary(cal$best_sim, data$dates)
write.csv(mon, file.path(csv_dir, "monthly_summary.csv"), row.names = FALSE)


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 8: MASS BALANCE
# ═══════════════════════════════════════════════════════════════════════════════

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  STEP 8: Mass balance check\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

mb <- check_mass_balance(cal$best_sim)


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 9: GENERATE PLOTS
# ═══════════════════════════════════════════════════════════════════════════════

if (SAVE_PLOTS) {
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("  STEP 9: Generating plots\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

  plots <- plot_all(data$dates, data$precip, data$Q_obs, cal, unc,
                    output_dir = plot_dir, prefix = RUN_NAME)
}


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 10: SAVE CSVs AND PARAMETERS
# ═══════════════════════════════════════════════════════════════════════════════

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  STEP 10: Saving results\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

if (SAVE_CSV_RESULTS) {
  export_results(cal, unc, sens,
                 dates = data$dates, Q_obs = data$Q_obs, precip = data$precip,
                 output_dir = csv_dir, prefix = RUN_NAME)
}

if (SAVE_PARAMETERS) {
  save_parameters(cal, file.path(run_dir, "parameters.txt"))
}


# ═══════════════════════════════════════════════════════════════════════════════
# DONE
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n")
cat("  ██████████████████████████████████████████████████████████████\n")
cat("  █                                                            █\n")
cat("  █    ✓ RUN COMPLETE                                          █\n")
cat("  █                                                            █\n")
cat(sprintf("  █    All results saved to: %-32s █\n", run_dir))
cat("  █                                                            █\n")
cat("  █    Folder contents:                                        █\n")
cat("  █      parameters.txt        (best k1, k2, k3 + metrics)    █\n")
cat("  █      config_used.R         (copy of your config)           █\n")
cat("  █      plots/                (PDF + individual PNGs)         █\n")
cat("  █      csv/                  (samples, behavioural, etc.)    █\n")
cat("  █                                                            █\n")
cat("  ██████████████████████████████████████████████████████████████\n\n")
