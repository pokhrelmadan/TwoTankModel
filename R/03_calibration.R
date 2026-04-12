#' @title Monte Carlo Calibration with Latin Hypercube Sampling

#' Calibrate the Two-Tank Model
#'
#' Generates \code{n_samples} parameter sets using Latin Hypercube Sampling,
#' runs the model for each, and compares to observed discharge.
#'
#' \strong{Units note:} If \code{area_km2} is provided, the model's simulated
#' Q is converted to m³/s before comparison with \code{Q_obs}. This is the
#' recommended approach when gauge data is in m³/s — you keep the original
#' units for reporting. If \code{area_km2} is NULL, both Q_obs and the
#' simulated Q must already be in mm/day.
#'
#' @param times Numeric vector. Time steps.
#' @param precip_vec Numeric vector. Daily precipitation \[mm/day\].
#' @param Q_obs Numeric vector. Observed daily discharge
#'   (m³/s if area_km2 is given, else mm/day).
#' @param n_samples Integer. Monte Carlo samples. Default 5000.
#' @param k1_range,k2_range,k3_range Numeric(2). Parameter bounds.
#' @param area_km2 Numeric. Catchment area \[km²\]. Default NULL.
#' @param objective Character. "NSE", "KGE", or "LogNSE". Default "NSE".
#' @param seed Integer. Random seed.
#' @param verbose Logical. Print progress.
#'
#' @return A list with best_params, best_sim, best_metrics, samples,
#'   n_samples, elapsed, area_km2, obs_units.
#'
#' @importFrom lhs randomLHS
#' @export
calibrate_montecarlo <- function(times, precip_vec, Q_obs,
                                 n_samples = 5000,
                                 k1_range = c(0.01, 0.80),
                                 k2_range = c(0.01, 0.50),
                                 k3_range = c(0.001, 0.15),
                                 area_km2 = NULL,
                                 objective = "NSE",
                                 seed = 42,
                                 verbose = TRUE) {

  if (length(times) != length(precip_vec))
    stop("times and precip_vec must have equal length.")
  if (length(times) != length(Q_obs))
    stop("times and Q_obs must have equal length.")
  if (n_samples < 10) stop("n_samples must be >= 10")
  if (any(Q_obs < 0)) { Q_obs[Q_obs < 0] <- 0; warning("Set negative Q_obs to 0.") }
  if (!objective %in% c("NSE", "KGE", "LogNSE"))
    stop("objective must be 'NSE', 'KGE', or 'LogNSE'.")

  obs_units <- if (!is.null(area_km2)) "m3s" else "mm_day"

  par_lower <- c(k1_range[1], k2_range[1], k3_range[1])
  par_upper <- c(k1_range[2], k2_range[2], k3_range[2])

  if (verbose) {
    cat("\n  ╔══════════════════════════════════════════════════════╗\n")
    cat("  ║          MONTE CARLO CALIBRATION (LHS)              ║\n")
    cat("  ╠══════════════════════════════════════════════════════╣\n")
    cat(sprintf("  ║  Samples        : %d\n", n_samples))
    cat(sprintf("  ║  Objective      : %s\n", objective))
    cat(sprintf("  ║  Observed units : %s\n", obs_units))
    if (!is.null(area_km2))
      cat(sprintf("  ║  Catchment area : %.1f km²\n", area_km2))
    cat(sprintf("  ║  k1 range       : [%.3f – %.3f]\n", k1_range[1], k1_range[2]))
    cat(sprintf("  ║  k2 range       : [%.3f – %.3f]\n", k2_range[1], k2_range[2]))
    cat(sprintf("  ║  k3 range       : [%.3f – %.3f]\n", k3_range[1], k3_range[2]))
    cat("  ╚══════════════════════════════════════════════════════╝\n\n")
  }

  set.seed(seed)
  lhs_mat <- lhs::randomLHS(n_samples, 3)

  samples <- data.frame(
    k1 = par_lower[1] + lhs_mat[,1] * (par_upper[1] - par_lower[1]),
    k2 = par_lower[2] + lhs_mat[,2] * (par_upper[2] - par_lower[2]),
    k3 = par_lower[3] + lhs_mat[,3] * (par_upper[3] - par_lower[3]),
    NSE = NA_real_, KGE = NA_real_, LogNSE = NA_real_,
    RMSE = NA_real_, PBIAS = NA_real_,
    Peak_Q = NA_real_, Vol_Q = NA_real_
  )

  t0 <- proc.time()
  update_at <- floor(seq(1, n_samples, length.out = 21))[-1]
  best_obj <- -Inf
  best_idx <- 1

  if (verbose) {
    cat("  Running simulations:\n")
    cat("  [", rep(" ", 40), "] 0%", sep = "")
  }

  for (i in 1:n_samples) {
    sim_i <- run_two_tank(samples$k1[i], samples$k2[i], samples$k3[i],
                          times, precip_vec, area_km2 = area_km2)

    # Compare in the user's units
    Q_sim <- if (!is.null(area_km2)) sim_i$Q_total_m3s else sim_i$Q_total

    samples$NSE[i]    <- calc_nse(Q_obs, Q_sim)
    samples$KGE[i]    <- calc_kge(Q_obs, Q_sim)
    samples$LogNSE[i] <- calc_lognse(Q_obs, Q_sim)
    samples$RMSE[i]   <- calc_rmse(Q_obs, Q_sim)
    samples$PBIAS[i]  <- calc_pbias(Q_obs, Q_sim)
    samples$Peak_Q[i] <- max(Q_sim)
    samples$Vol_Q[i]  <- sum(Q_sim)

    current_obj <- switch(objective,
      "NSE"    = samples$NSE[i],
      "KGE"    = samples$KGE[i],
      "LogNSE" = samples$LogNSE[i])

    if (current_obj > best_obj) {
      best_obj <- current_obj
      best_idx <- i
    }

    if (verbose && i %in% update_at) {
      pct <- round(i / n_samples * 100)
      bar <- round(i / n_samples * 40)
      cat("\r  [", rep("█", bar), rep(" ", 40 - bar), "] ",
          pct, "%  (best ", objective, ": ",
          sprintf("%.4f", best_obj), ")", sep = "")
    }
  }

  elapsed <- (proc.time() - t0)[3]
  if (verbose) cat("\n\n")

  best_params <- unlist(samples[best_idx, 1:3])
  names(best_params) <- c("k1", "k2", "k3")
  best_sim <- run_two_tank(best_params[1], best_params[2], best_params[3],
                           times, precip_vec, area_km2 = area_km2)

  best_metrics <- list(
    NSE    = samples$NSE[best_idx],
    KGE    = samples$KGE[best_idx],
    LogNSE = samples$LogNSE[best_idx],
    RMSE   = samples$RMSE[best_idx],
    PBIAS  = samples$PBIAS[best_idx]
  )

  if (verbose) {
    u <- if (!is.null(area_km2)) "m³/s" else "mm/day"
    cat("  ┌─────────────────────────────────────────────────────┐\n")
    cat("  │               CALIBRATION RESULTS                   │\n")
    cat("  ├─────────────────────────────────────────────────────┤\n")
    cat(sprintf("  │  k1 = %.4f  (surface runoff)                       │\n", best_params[1]))
    cat(sprintf("  │  k2 = %.4f  (percolation)                          │\n", best_params[2]))
    cat(sprintf("  │  k3 = %.4f  (baseflow)                            │\n", best_params[3]))
    cat("  ├─────────────────────────────────────────────────────┤\n")
    cat(sprintf("  │  NSE    = %.4f                                    │\n", best_metrics$NSE))
    cat(sprintf("  │  KGE    = %.4f                                    │\n", best_metrics$KGE))
    cat(sprintf("  │  LogNSE = %.4f                                    │\n", best_metrics$LogNSE))
    cat(sprintf("  │  RMSE   = %.4f %s                            │\n", best_metrics$RMSE, u))
    cat(sprintf("  │  PBIAS  = %.2f %%                                 │\n", best_metrics$PBIAS))
    cat("  ├─────────────────────────────────────────────────────┤\n")
    cat(sprintf("  │  %d simulations in %.1f sec (%.0f sims/s)         │\n",
                n_samples, elapsed, n_samples/elapsed))
    cat("  └─────────────────────────────────────────────────────┘\n\n")
  }

  return(list(
    best_params  = best_params,
    best_sim     = best_sim,
    best_metrics = best_metrics,
    samples      = samples,
    n_samples    = n_samples,
    elapsed      = elapsed,
    area_km2     = area_km2,
    obs_units    = obs_units,
    objective    = objective
  ))
}
