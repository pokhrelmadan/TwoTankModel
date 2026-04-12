################################################################################
#  GENERATE SAMPLE DATA
#  ────────────────────
#  Creates sample rainfall and discharge CSV files so you can test the model
#  immediately. Run this script once to create data/rainfall.csv and
#  data/discharge.csv, then run main.R.
#
#  Usage:  Rscript generate_sample_data.R
################################################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!require("deSolve", quietly = TRUE)) install.packages("deSolve")
library(deSolve)

# Source the model functions
for (f in list.files("R", full.names = TRUE, pattern = "\\.R$")) source(f)

# Generate one year of data
dates  <- seq(as.Date("2024-01-01"), as.Date("2024-12-31"), by = "day")
n_days <- length(dates)
times  <- 0:(n_days - 1)

# Synthetic rainfall with seasonal pattern
precip <- generate_daily_rainfall(n_days, dates)

# "True" parameters (what we're trying to recover during calibration)
true_params <- c(k1 = 0.30, k2 = 0.12, k3 = 0.025)
area_km2    <- 150   # example catchment area

# Run model to get discharge, then add some noise to simulate observations
sim <- run_two_tank(true_params[1], true_params[2], true_params[3],
                    times, precip, area_km2 = area_km2)
set.seed(123)
Q_m3s <- pmax(0, sim$Q_total_m3s +
              rnorm(n_days, 0, 0.08 * max(sim$Q_total_m3s)))

# Create data folder and write CSVs
dir.create("data", showWarnings = FALSE)

write.csv(data.frame(date = dates, P = precip),
          "data/rainfall.csv", row.names = FALSE)

write.csv(data.frame(date = dates, Q = round(Q_m3s, 3)),
          "data/discharge.csv", row.names = FALSE)

cat("\n  Sample data created:\n")
cat("    data/rainfall.csv    (", n_days, "days,",
    sum(precip), "mm total)\n")
cat("    data/discharge.csv   (", n_days, "days, in m³/s)\n")
cat("\n  True parameters used:\n")
cat(sprintf("    k1 = %.3f\n    k2 = %.3f\n    k3 = %.3f\n",
            true_params[1], true_params[2], true_params[3]))
cat(sprintf("    Catchment area = %.1f km²\n", area_km2))
cat("\n  Now run:  Rscript main.R\n\n")
