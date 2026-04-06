################################################################################
#
#    ████████╗██╗    ██╗ ██████╗    ████████╗ █████╗ ███╗   ██╗██╗  ██╗
#    ╚══██╔══╝██║    ██║██╔═══██╗   ╚══██╔══╝██╔══██╗████╗  ██║██║ ██╔╝
#       ██║   ██║ █╗ ██║██║   ██║      ██║   ███████║██╔██╗ ██║█████╔╝
#       ██║   ██║███╗██║██║   ██║      ██║   ██╔══██║██║╚██╗██║██╔═██╗
#       ██║   ╚███╔███╔╝╚██████╔╝      ██║   ██║  ██║██║ ╚████║██║  ██╗
#       ╚═╝    ╚══╝╚══╝  ╚═════╝       ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═══╝╚═╝  ╚═╝
#
#    Two-Tank Conceptual Rainfall-Runoff Model
#    Monte Carlo Calibration with Latin Hypercube Sampling
#
#    COMPLETE WORKFLOW — Run this script from top to bottom.
#    Each step is clearly marked. Modify only the USER SETTINGS below.
#
################################################################################


# ═══════════════════════════════════════════════════════════════════════════════
# USER SETTINGS — CHANGE THESE TO MATCH YOUR DATA
# ═══════════════════════════════════════════════════════════════════════════════

# OPTION 1: Use your own data files
#   Set use_own_data <- TRUE and fill in the file paths below.

use_own_data <- FALSE    # ← Change to TRUE when you have real data

rainfall_file  <- "daily_rainfall.csv"    # Columns: date, P (mm/day)
discharge_file <- "daily_discharge.csv"   # Columns: date, Q (mm/day)
catchment_area <- NULL                    # Set to area in km² if Q is in m³/s

# OPTION 2: Use synthetic demo data (default)
#   Leave use_own_data <- FALSE to run with generated data.

# Calibration settings
n_mc_samples   <- 5000     # Monte Carlo samples (1000=quick, 5000=standard, 10000=fine)
nse_threshold  <- 0.5      # Behavioural threshold

# Parameter search bounds [1/day]
k1_bounds <- c(0.01, 0.80)    # Surface runoff coefficient
k2_bounds <- c(0.01, 0.50)    # Percolation coefficient
k3_bounds <- c(0.001, 0.15)   # Baseflow coefficient

# Output
output_directory <- "."       # Where to save plots and CSVs
output_prefix    <- "twotank" # Filename prefix


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 0: Load the package
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n")
cat("  ╔══════════════════════════════════════════════════════════╗\n")
cat("  ║    TWO-TANK MODEL — STARTING WORKFLOW                   ║\n")
cat("  ╚══════════════════════════════════════════════════════════╝\n\n")

# Load all functions (use library(TwoTankModel) if installed as package)
for (f in list.files("R", full.names = TRUE, pattern = "\\.R$")) source(f)
cat("  ✓ All functions loaded.\n\n")


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1: Load or generate data
# ═══════════════════════════════════════════════════════════════════════════════

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  STEP 1: Loading data\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

if (use_own_data) {
  # ── Load your CSV files ──
  data <- load_data(rainfall_file, discharge_file, area_km2 = catchment_area)
  dates  <- data$dates
  times  <- data$times
  precip <- data$precip
  Q_obs  <- data$Q_obs
  n_days <- data$n_days

} else {
  # ── Generate synthetic demo data ──
  cat("  Using synthetic demo data (set use_own_data <- TRUE for real data).\n\n")

  dates  <- seq(as.Date("2024-01-01"), as.Date("2024-12-31"), by = "day")
  n_days <- length(dates)
  times  <- 0:(n_days - 1)
  precip <- generate_daily_rainfall(n_days, dates)

  # Create synthetic observations with known true parameters
  true_params <- c(k1 = 0.30, k2 = 0.12, k3 = 0.025)
  sim_true    <- run_two_tank(true_params[1], true_params[2], true_params[3],
                              times, precip)
  set.seed(123)
  Q_obs <- pmax(0, sim_true$Q_total + rnorm(n_days, 0, 0.10 * max(sim_true$Q_total)))

  cat(sprintf("  Synthetic rainfall: %d days | Total = %.0f mm\n", n_days, sum(precip)))
  cat(sprintf("  True parameters: k1=%.3f  k2=%.3f  k3=%.3f\n\n",
              true_params[1], true_params[2], true_params[3]))
}


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2: Calibrate the model
# ═══════════════════════════════════════════════════════════════════════════════

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  STEP 2: Calibrating model (Monte Carlo with LHS)\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

cal <- calibrate_montecarlo(
  times      = times,
  precip_vec = precip,
  Q_obs      = Q_obs,
  n_samples  = n_mc_samples,
  k1_range   = k1_bounds,
  k2_range   = k2_bounds,
  k3_range   = k3_bounds
)


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3: Analyse uncertainty
# ═══════════════════════════════════════════════════════════════════════════════

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  STEP 3: Uncertainty analysis\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

unc <- extract_uncertainty(
  cal_result    = cal,
  times         = times,
  precip_vec    = precip,
  Q_obs         = Q_obs,
  nse_threshold = nse_threshold
)


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4: Performance report
# ═══════════════════════════════════════════════════════════════════════════════

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  STEP 4: Performance evaluation\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

metrics <- calc_all_metrics(Q_obs, cal$best_sim$Q_total)


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 5: Hydrograph analysis
# ═══════════════════════════════════════════════════════════════════════════════

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  STEP 5: Hydrograph metrics\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

hydro_met <- extract_metrics(cal$best_sim, dates, label = "Calibrated Model")
print_model_summary(cal$best_sim, label = "Best Calibration")


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 6: Monthly water balance
# ═══════════════════════════════════════════════════════════════════════════════

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  STEP 6: Monthly water balance\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

mon <- monthly_summary(cal$best_sim, dates)


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 7: Mass balance check
# ═══════════════════════════════════════════════════════════════════════════════

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  STEP 7: Mass balance verification\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

mb <- check_mass_balance(cal$best_sim)


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 8: Generate all plots
# ═══════════════════════════════════════════════════════════════════════════════

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  STEP 8: Generating diagnostic plots\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

plots <- plot_all(dates, precip, Q_obs, cal, unc,
                  output_dir = output_directory, prefix = output_prefix)


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 9: Export results
# ═══════════════════════════════════════════════════════════════════════════════

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  STEP 9: Exporting results to CSV\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

export_results(cal, unc, dates, Q_obs,
               output_dir = output_directory, prefix = output_prefix)


# ═══════════════════════════════════════════════════════════════════════════════
# DONE!
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n")
cat("  ╔══════════════════════════════════════════════════════════╗\n")
cat("  ║                                                        ║\n")
cat("  ║    ✓ WORKFLOW COMPLETE                                  ║\n")
cat("  ║                                                        ║\n")
cat("  ║    Results saved in: ", sprintf("%-33s", output_directory),   "║\n")
cat("  ║                                                        ║\n")
cat("  ║    Files created:                                      ║\n")
cat(sprintf("  ║      %s_results.pdf         (all plots)   ║\n", output_prefix))
cat(sprintf("  ║      %s_mc_samples.csv      (all MC data) ║\n", output_prefix))
cat(sprintf("  ║      %s_behavioural.csv     (good sets)   ║\n", output_prefix))
cat(sprintf("  ║      %s_best_simulation.csv (daily Q)     ║\n", output_prefix))
cat("  ║      + 12 individual PNG plots                        ║\n")
cat("  ║                                                        ║\n")
cat("  ║    Next steps:                                         ║\n")
cat("  ║      • Set use_own_data <- TRUE with your CSV files    ║\n")
cat("  ║      • Adjust k1/k2/k3 bounds if needed               ║\n")
cat("  ║      • Increase n_mc_samples for publication           ║\n")
cat("  ║      • Use compare_watersheds() for two catchments     ║\n")
cat("  ║                                                        ║\n")
cat("  ╚══════════════════════════════════════════════════════════╝\n\n")
