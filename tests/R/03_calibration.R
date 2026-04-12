#' @title Monte Carlo Calibration with Latin Hypercube Sampling
#' @description Calibrates the two-tank model by sampling the parameter
#'   space using Latin Hypercube Sampling (LHS) and evaluating each
#'   sample against observed discharge.


#' Calibrate the Two-Tank Model
#'
#' Generates \code{n_samples} parameter sets using Latin Hypercube
#' Sampling, runs the model for each, and finds the best-performing
#' set. Also records all results for uncertainty analysis.
#'
#' \strong{Quick start:}
#' \preformatted{
#' cal <- calibrate_montecarlo(times, precip, Q_obs)
#' print(cal$best_params)   # best k1, k2, k3
#' }
#'
#' @param times Numeric vector. Time steps (0, 1, 2, ...).
#' @param precip_vec Numeric vector. Daily precipitation \[mm/day\].
#' @param Q_obs Numeric vector. Observed daily discharge \[mm/day\].
#' @param n_samples Integer. Number of Monte Carlo samples.
#'   Default 5000. Use 1000 for quick tests, 10000+ for publication.
#' @param k1_range Numeric vector of length 2. Min/max for k1.
#'   Default c(0.01, 0.80).
#' @param k2_range Numeric vector of length 2. Min/max for k2.
#'   Default c(0.01, 0.50).
#' @param k3_range Numeric vector of length 2. Min/max for k3.
#'   Default c(0.001, 0.15).
#' @param seed Integer. Random seed for reproducibility. Default 42.
#' @param verbose Logical. Print progress. Default TRUE.
#'
#' @return A list with:
#'   \describe{
#'     \item{best_params}{Named vector: best k1, k2, k3}
#'     \item{best_sim}{Data.frame: simulation with best parameters}
#'     \item{best_metrics}{List: NSE, KGE, LogNSE, RMSE, PBIAS of best}
#'     \item{samples}{Data.frame: all samples with metrics (for uncertainty)}
#'     \item{n_samples}{Number of samples run}
#'     \item{elapsed}{Runtime in seconds}
#'   }
#'
#' @examples
#' \dontrun{
#' cal <- calibrate_montecarlo(times, precip, Q_obs, n_samples = 2000)
#' print(cal$best_params)
#' }
#'
#' @importFrom lhs randomLHS
#' @export
calibrate_montecarlo <- function(times, precip_vec, Q_obs,
                                 n_samples = 5000,
                                 k1_range = c(0.01, 0.80),
                                 k2_range = c(0.01, 0.50),
                                 k3_range = c(0.001, 0.15),
                                 seed = 42,
                                 verbose = TRUE) {

  # ── Input validation with helpful messages ──
  if (length(times) != length(precip_vec)) {
    stop("'times' and 'precip_vec' must have the same length.\n",
         "  Got: times=", length(times), ", precip_vec=", length(precip_vec))
  }
  if (length(times) != length(Q_obs)) {
    stop("'times' and 'Q_obs' must have the same length.\n",
         "  Got: times=", length(times), ", Q_obs=", length(Q_obs))
  }
  if (n_samples < 10) {
    stop("n_samples must be at least 10. Got: ", n_samples, "\n",
         "  Recommended: 1000 (quick), 5000 (standard), 10000+ (publication)")
  }
  if (any(Q_obs < 0)) {
    warning("Q_obs contains ", sum(Q_obs < 0), " negative values. ",
            "Setting them to 0.")
    Q_obs[Q_obs < 0] <- 0
  }

  par_lower <- c(k1_range[1], k2_range[1], k3_range[1])
  par_upper <- c(k1_range[2], k2_range[2], k3_range[2])

  # ── Header ──
  if (verbose) {
    cat("\n")
    cat("  ╔══════════════════════════════════════════════════════╗\n")
    cat("  ║          MONTE CARLO CALIBRATION (LHS)              ║\n")
    cat("  ╠══════════════════════════════════════════════════════╣\n")
    cat(sprintf("  ║  Samples        : %d\n", n_samples))
    cat(sprintf("  ║  k1 range       : [%.3f – %.3f]  surface runoff\n",
                k1_range[1], k1_range[2]))
    cat(sprintf("  ║  k2 range       : [%.3f – %.3f]  percolation\n",
                k2_range[1], k2_range[2]))
    cat(sprintf("  ║  k3 range       : [%.3f – %.3f]  baseflow\n",
                k3_range[1], k3_range[2]))
    cat(sprintf("  ║  Data length    : %d days\n", length(times)))
    cat("  ╚══════════════════════════════════════════════════════╝\n\n")
  }

  # ── Generate LHS samples ──
  set.seed(seed)
  lhs_mat <- lhs::randomLHS(n_samples, 3)

  samples <- data.frame(
    k1     = par_lower[1] + lhs_mat[, 1] * (par_upper[1] - par_lower[1]),
    k2     = par_lower[2] + lhs_mat[, 2] * (par_upper[2] - par_lower[2]),
    k3     = par_lower[3] + lhs_mat[, 3] * (par_upper[3] - par_lower[3]),
    NSE    = NA_real_,
    KGE    = NA_real_,
    LogNSE = NA_real_,
    RMSE   = NA_real_,
    PBIAS  = NA_real_,
    Peak_Q = NA_real_,
    Vol_Q  = NA_real_
  )

  # ── Run all samples with progress bar ──
  t0 <- proc.time()
  n_updates <- 20  # number of progress bar updates
  update_at <- floor(seq(1, n_samples, length.out = n_updates + 1))[-1]

  if (verbose) {
    cat("  Running simulations:\n")
    cat("  [", rep(" ", 40), "] 0%", sep = "")
  }

  best_nse <- -Inf
  best_idx <- 1

  for (i in 1:n_samples) {
    sim_i <- run_two_tank(samples$k1[i], samples$k2[i], samples$k3[i],
                          times, precip_vec)
    Q_sim <- sim_i$Q_total

    samples$NSE[i]    <- calc_nse(Q_obs, Q_sim)
    samples$KGE[i]    <- calc_kge(Q_obs, Q_sim)
    samples$LogNSE[i] <- calc_lognse(Q_obs, Q_sim)
    samples$RMSE[i]   <- calc_rmse(Q_obs, Q_sim)
    samples$PBIAS[i]  <- calc_pbias(Q_obs, Q_sim)
    samples$Peak_Q[i] <- max(Q_sim)
    samples$Vol_Q[i]  <- sum(Q_sim)

    if (samples$NSE[i] > best_nse) {
      best_nse <- samples$NSE[i]
      best_idx <- i
    }

    # Update progress bar
    if (verbose && i %in% update_at) {
      pct  <- round(i / n_samples * 100)
      bar  <- round(i / n_samples * 40)
      cat("\r  [", rep("█", bar), rep(" ", 40 - bar), "] ",
          pct, "%  (best NSE so far: ", sprintf("%.4f", best_nse), ")",
          sep = "")
    }
  }

  elapsed <- (proc.time() - t0)[3]
  if (verbose) cat("\n\n")

  # ── Extract best result ──
  best_params <- unlist(samples[best_idx, 1:3])
  names(best_params) <- c("k1", "k2", "k3")
  best_sim <- run_two_tank(best_params[1], best_params[2], best_params[3],
                           times, precip_vec)

  best_metrics <- list(
    NSE    = samples$NSE[best_idx],
    KGE    = samples$KGE[best_idx],
    LogNSE = samples$LogNSE[best_idx],
    RMSE   = samples$RMSE[best_idx],
    PBIAS  = samples$PBIAS[best_idx]
  )

  # ── Print results ──
  if (verbose) {
    cat("  ┌─────────────────────────────────────────────────────┐\n")
    cat("  │               CALIBRATION RESULTS                   │\n")
    cat("  ├─────────────────────────────────────────────────────┤\n")
    cat(sprintf("  │  Best k1 = %.4f  (surface runoff)                 │\n",
                best_params[1]))
    cat(sprintf("  │  Best k2 = %.4f  (percolation)                    │\n",
                best_params[2]))
    cat(sprintf("  │  Best k3 = %.4f  (baseflow)                      │\n",
                best_params[3]))
    cat("  ├─────────────────────────────────────────────────────┤\n")
    cat(sprintf("  │  Upper tank residence ≈ %.1f days                  │\n",
                1 / (best_params[1] + best_params[2])))
    cat(sprintf("  │  Lower tank residence ≈ %.1f days                  │\n",
                1 / best_params[3]))
    cat("  ├─────────────────────────────────────────────────────┤\n")
    cat(sprintf("  │  NSE    = %.4f                                    │\n", best_metrics$NSE))
    cat(sprintf("  │  KGE    = %.4f                                    │\n", best_metrics$KGE))
    cat(sprintf("  │  LogNSE = %.4f                                    │\n", best_metrics$LogNSE))
    cat(sprintf("  │  RMSE   = %.4f mm/day                             │\n", best_metrics$RMSE))
    cat(sprintf("  │  PBIAS  = %.2f %%                                  │\n", best_metrics$PBIAS))
    cat("  ├─────────────────────────────────────────────────────┤\n")
    cat(sprintf("  │  Completed %d samples in %.1f seconds              │\n",
                n_samples, elapsed))
    cat(sprintf("  │  Speed: %.0f simulations/second                    │\n",
                n_samples / elapsed))
    cat("  └─────────────────────────────────────────────────────┘\n\n")
  }

  return(list(
    best_params  = best_params,
    best_sim     = best_sim,
    best_metrics = best_metrics,
    samples      = samples,
    n_samples    = n_samples,
    elapsed      = elapsed
  ))
}
