################################################################################
#  QUICK START — Run this to see the model in action in under 1 minute.
#  No data files needed. Just source and go!
################################################################################

# Load functions
for (f in list.files("R", full.names = TRUE, pattern = "\\.R$")) source(f)

# Generate demo data
dates  <- seq(as.Date("2024-01-01"), as.Date("2024-12-31"), by = "day")
times  <- 0:(length(dates) - 1)
precip <- generate_daily_rainfall(length(dates), dates)

# Create fake observed discharge
sim_true <- run_two_tank(0.30, 0.12, 0.025, times, precip)
Q_obs    <- pmax(0, sim_true$Q_total + rnorm(length(dates), 0, 0.3))

# Calibrate (1000 samples for speed — use 5000+ for real work)
cal <- calibrate_montecarlo(times, precip, Q_obs, n_samples = 1000)

# See how well it did
calc_all_metrics(Q_obs, cal$best_sim$Q_total)

# Check mass balance
check_mass_balance(cal$best_sim)

# Plot results (saves to current directory)
unc <- extract_uncertainty(cal, times, precip, Q_obs)
plots <- plot_all(dates, precip, Q_obs, cal, unc)

cat("\n  Done! Check your working directory for plots.\n")
