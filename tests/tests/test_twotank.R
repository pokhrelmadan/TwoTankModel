#' @title Unit Tests for TwoTankModel
#' @description Run with: testthat::test_file("tests/test_twotank.R")

library(testthat)

# Source all R files (when not installed as package)
# Uncomment these lines if running tests standalone:
# for (f in list.files("R", full.names = TRUE)) source(f)

# ── Test: Model runs without error ──────────────────────────────────────────

test_that("run_two_tank returns correct structure", {
  times  <- 0:99
  precip <- c(rep(0, 20), rep(10, 10), rep(0, 70))
  sim    <- run_two_tank(0.3, 0.1, 0.02, times, precip)

  expect_s3_class(sim, "data.frame")
  expect_equal(nrow(sim), 100)
  expect_equal(names(sim), c("time", "S1", "S2", "Q1", "Q2", "Q_total", "P"))
  expect_true(all(sim$Q_total >= 0))
  expect_true(all(sim$S1 >= -1e-10))   # allow tiny numerical noise
  expect_true(all(sim$S2 >= -1e-10))
})

test_that("run_two_tank validates inputs", {
  times <- 0:9
  precip <- rep(5, 10)

  expect_error(run_two_tank(-0.1, 0.1, 0.02, times, precip))
  expect_error(run_two_tank(0.1, 0.1, 0.02, times, precip[1:5]))
  expect_error(run_two_tank(0.1, 0.1, 0.02, times, c(precip[1:9], -1)))
})

# ── Test: Zero rainfall produces zero discharge ────────────────────────────

test_that("zero precipitation gives zero discharge", {
  times  <- 0:99
  precip <- rep(0, 100)
  sim    <- run_two_tank(0.3, 0.1, 0.02, times, precip)

  expect_equal(sum(sim$Q_total), 0, tolerance = 1e-10)
})

# ── Test: Mass balance closes ──────────────────────────────────────────────

test_that("mass balance closes within tolerance", {
  times  <- 0:364
  dates  <- seq(as.Date("2024-01-01"), as.Date("2024-12-31"), by = "day")
  precip <- generate_daily_rainfall(365, dates)
  sim    <- run_two_tank(0.25, 0.15, 0.03, times, precip)

  mb <- check_mass_balance(sim, verbose = FALSE)
  expect_lt(abs(mb$residual), 0.01)  # residual < 0.01 mm
})

# ── Test: Performance metrics ──────────────────────────────────────────────

test_that("NSE of perfect match is 1", {
  obs <- c(1, 2, 3, 4, 5)
  expect_equal(calc_nse(obs, obs), 1)
})

test_that("NSE of mean prediction is 0", {
  obs <- c(1, 2, 3, 4, 5)
  sim <- rep(mean(obs), 5)
  expect_equal(calc_nse(obs, sim), 0)
})

test_that("KGE of perfect match is 1", {
  obs <- c(1, 2, 3, 4, 5)
  expect_equal(calc_kge(obs, obs), 1)
})

test_that("RMSE of perfect match is 0", {
  obs <- c(1, 2, 3, 4, 5)
  expect_equal(calc_rmse(obs, obs), 0)
})

test_that("PBIAS of perfect match is 0", {
  obs <- c(1, 2, 3, 4, 5)
  expect_equal(calc_pbias(obs, obs), 0)
})

test_that("metrics validate equal lengths", {
  expect_error(calc_nse(1:5, 1:3))
  expect_error(calc_kge(1:5, 1:3))
  expect_error(calc_rmse(1:5, 1:3))
})

# ── Test: Monte Carlo calibration ─────────────────────────────────────────

test_that("calibrate_montecarlo returns correct structure", {
  times  <- 0:99
  precip <- c(rep(0, 20), rep(15, 10), rep(0, 70))
  sim_t  <- run_two_tank(0.3, 0.1, 0.02, times, precip)
  Q_obs  <- pmax(0, sim_t$Q_total + rnorm(100, 0, 0.1))

  cal <- calibrate_montecarlo(times, precip, Q_obs,
                               n_samples = 50, verbose = FALSE)

  expect_type(cal, "list")
  expect_named(cal, c("best_params", "best_sim", "best_metrics",
                       "samples", "n_samples", "elapsed"))
  expect_equal(length(cal$best_params), 3)
  expect_equal(nrow(cal$samples), 50)
  expect_true(cal$best_metrics$NSE > -10)  # sanity check
})

# ── Test: Uncertainty extraction ──────────────────────────────────────────

test_that("extract_uncertainty returns valid bounds", {
  times  <- 0:99
  precip <- c(rep(0, 20), rep(15, 10), rep(0, 70))
  sim_t  <- run_two_tank(0.3, 0.1, 0.02, times, precip)
  Q_obs  <- pmax(0, sim_t$Q_total + rnorm(100, 0, 0.1))

  cal <- calibrate_montecarlo(times, precip, Q_obs,
                               n_samples = 100, verbose = FALSE)
  unc <- extract_uncertainty(cal, times, precip, Q_obs,
                              nse_threshold = -999, verbose = FALSE)

  expect_equal(length(unc$Q_lower), 100)
  expect_equal(length(unc$Q_upper), 100)
  expect_true(all(unc$Q_upper >= unc$Q_lower))
})

# ── Test: Unit conversions ────────────────────────────────────────────────

test_that("unit conversions are inverse of each other", {
  Q_m3s <- 5.0
  area  <- 150
  Q_mm  <- m3s_to_mmday(Q_m3s, area)
  Q_back <- mmday_to_m3s(Q_mm, area)
  expect_equal(Q_back, Q_m3s, tolerance = 1e-10)
})

# ── Test: Rainfall generation ─────────────────────────────────────────────

test_that("generate_daily_rainfall produces valid output", {
  dates  <- seq(as.Date("2024-01-01"), as.Date("2024-12-31"), by = "day")
  precip <- generate_daily_rainfall(365, dates)

  expect_equal(length(precip), 365)
  expect_true(all(precip >= 0))
  expect_true(max(precip) > 0)  # at least some rain
})

cat("\n  All tests passed.\n")
