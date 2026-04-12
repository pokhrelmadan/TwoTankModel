#' @title Uncertainty Analysis
#' @description Extract behavioural parameter sets and compute prediction
#'   uncertainty envelopes from Monte Carlo calibration results.


#' Analyse Parameter Uncertainty
#'
#' Identifies 'behavioural' parameter sets (those with NSE above a
#' threshold), computes prediction uncertainty envelopes (90\% bounds),
#' and reports how well the bounds capture the observations.
#'
#' \strong{What are behavioural sets?}
#' In hydrology, a parameter set is called 'behavioural' if it produces
#' an acceptable simulation (traditionally NSE > 0.5). The range of
#' behavioural parameters tells you how uncertain each parameter is.
#'
#' @param cal_result List. Output from \code{\link{calibrate_montecarlo}}.
#' @param times Numeric vector. Time steps.
#' @param precip_vec Numeric vector. Daily precipitation \[mm/day\].
#' @param Q_obs Numeric vector. Observed daily discharge \[mm/day\].
#' @param nse_threshold Numeric. Minimum NSE for behavioural.
#'   Default 0.5 (standard in hydrology).
#' @param max_behavioural Integer. Cap on sets to run for speed. Default 500.
#' @param verbose Logical. Print results. Default TRUE.
#'
#' @return A list with:
#'   \describe{
#'     \item{behavioural}{Data.frame of behavioural parameter sets}
#'     \item{n_behavioural}{Count of behavioural sets}
#'     \item{param_ranges}{Data.frame: min/max for k1, k2, k3}
#'     \item{Q_lower}{5th percentile envelope \[mm/day\]}
#'     \item{Q_upper}{95th percentile envelope \[mm/day\]}
#'     \item{Q_median}{Median ensemble prediction \[mm/day\]}
#'     \item{containment_pct}{\% of observations within 90\% bounds}
#'   }
#'
#' @export
extract_uncertainty <- function(cal_result, times, precip_vec, Q_obs,
                                nse_threshold = 0.5,
                                max_behavioural = 500,
                                verbose = TRUE) {

  samples <- cal_result$samples
  n_days  <- length(times)

  # ── Identify behavioural sets ──
  behav   <- samples[samples$NSE >= nse_threshold, ]
  behav   <- behav[order(-behav$NSE), ]
  n_behav <- nrow(behav)
  pct     <- round(n_behav / nrow(samples) * 100, 1)

  # Parameter ranges
  param_ranges <- data.frame(
    parameter = c("k1 (surface runoff)", "k2 (percolation)", "k3 (baseflow)"),
    min   = round(c(min(behav$k1), min(behav$k2), min(behav$k3)), 4),
    max   = round(c(max(behav$k1), max(behav$k2), max(behav$k3)), 4),
    range = round(c(max(behav$k1) - min(behav$k1),
                    max(behav$k2) - min(behav$k2),
                    max(behav$k3) - min(behav$k3)), 4),
    stringsAsFactors = FALSE
  )

  if (verbose) {
    cat("\n")
    cat("  ╔══════════════════════════════════════════════════════╗\n")
    cat("  ║           PARAMETER UNCERTAINTY ANALYSIS            ║\n")
    cat("  ╠══════════════════════════════════════════════════════╣\n")
    cat(sprintf("  ║  NSE threshold      : >= %.2f\n", nse_threshold))
    cat(sprintf("  ║  Behavioural sets    : %d of %d (%s%%)\n",
                n_behav, nrow(samples), pct))
    cat("  ║\n")
    cat("  ║  Parameter ranges (behavioural):\n")
    for (i in 1:3) {
      cat(sprintf("  ║    %-20s: [%.4f – %.4f]  (range %.4f)\n",
                  param_ranges$parameter[i], param_ranges$min[i],
                  param_ranges$max[i], param_ranges$range[i]))
    }
  }

  # ── Build uncertainty envelope ──
  n_run <- min(n_behav, max_behavioural)

  if (verbose) {
    cat("  ║\n")
    cat(sprintf("  ║  Computing 90%% envelope from %d simulations...\n", n_run))
  }

  Q_ensemble <- matrix(NA, nrow = n_days, ncol = n_run)
  for (b in 1:n_run) {
    sim_b <- run_two_tank(behav$k1[b], behav$k2[b], behav$k3[b],
                          times, precip_vec)
    Q_ensemble[, b] <- sim_b$Q_total
  }

  Q_lower  <- apply(Q_ensemble, 1, quantile, probs = 0.05)
  Q_upper  <- apply(Q_ensemble, 1, quantile, probs = 0.95)
  Q_median <- apply(Q_ensemble, 1, median)

  containment <- sum(Q_obs >= Q_lower & Q_obs <= Q_upper) / n_days * 100

  if (verbose) {
    cat(sprintf("  ║  Observations within 90%% bounds: %.1f%%\n", containment))
    cat("  ║\n")

    # Interpretation
    if (containment >= 80) {
      cat("  ║  ✓ Good: envelope captures most observations.\n")
    } else if (containment >= 50) {
      cat("  ║  ~ Fair: consider widening bounds or more samples.\n")
    } else {
      cat("  ║  ✗ Poor: model may have structural issues.\n")
    }

    # Identifiability
    cat("  ║\n")
    cat("  ║  Parameter identifiability:\n")
    for (i in 1:3) {
      total_range <- c(diff(c(0.01, 0.80)), diff(c(0.01, 0.50)),
                       diff(c(0.001, 0.15)))[i]
      ident_pct <- (1 - param_ranges$range[i] / total_range) * 100
      if (ident_pct > 80) tag <- "Well identified ✓"
      else if (ident_pct > 50) tag <- "Moderately identified ~"
      else tag <- "Poorly identified ✗ (equifinal)"
      cat(sprintf("  ║    %-20s: %s\n", param_ranges$parameter[i], tag))
    }

    cat("  ╚══════════════════════════════════════════════════════╝\n\n")
  }

  return(list(
    behavioural     = behav,
    n_behavioural   = n_behav,
    param_ranges    = param_ranges,
    Q_lower         = Q_lower,
    Q_upper         = Q_upper,
    Q_median        = Q_median,
    containment_pct = containment,
    nse_threshold   = nse_threshold
  ))
}
