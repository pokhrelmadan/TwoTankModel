#' @title Uncertainty Analysis

#' Analyse Parameter Uncertainty
#'
#' Identifies behavioural parameter sets (NSE >= threshold) and
#' builds 90% prediction uncertainty bounds.
#'
#' @param cal_result List. Output from \code{\link{calibrate_montecarlo}}.
#' @param times Numeric vector. Time steps.
#' @param precip_vec Numeric vector. Daily precipitation \[mm/day\].
#' @param Q_obs Numeric vector. Observed daily discharge.
#' @param nse_threshold Numeric. Default 0.5.
#' @param max_behavioural Integer. Max sets for envelope. Default 500.
#' @param verbose Logical. Default TRUE.
#'
#' @return List with behavioural, n_behavioural, param_ranges,
#'   Q_lower, Q_upper, Q_median, containment_pct, nse_threshold.
#' @export
extract_uncertainty <- function(cal_result, times, precip_vec, Q_obs,
                                nse_threshold = 0.5,
                                max_behavioural = 500,
                                verbose = TRUE) {

  samples  <- cal_result$samples
  area_km2 <- cal_result$area_km2
  n_days   <- length(times)

  behav   <- samples[samples$NSE >= nse_threshold, ]
  behav   <- behav[order(-behav$NSE), ]
  n_behav <- nrow(behav)
  pct     <- round(n_behav / nrow(samples) * 100, 1)

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
    cat("\n  ╔══════════════════════════════════════════════════════╗\n")
    cat("  ║          UNCERTAINTY ANALYSIS                       ║\n")
    cat("  ╠══════════════════════════════════════════════════════╣\n")
    cat(sprintf("  ║  Threshold      : NSE >= %.2f\n", nse_threshold))
    cat(sprintf("  ║  Behavioural    : %d of %d (%s%%)\n",
                n_behav, nrow(samples), pct))
    cat("  ║  Parameter uncertainty ranges:\n")
    for (i in 1:3)
      cat(sprintf("  ║    %-22s: [%.4f – %.4f]\n",
                  param_ranges$parameter[i], param_ranges$min[i],
                  param_ranges$max[i]))
  }

  n_run <- min(n_behav, max_behavioural)
  if (verbose) cat(sprintf("  ║  Computing envelope (%d sims)...\n", n_run))

  Q_ensemble <- matrix(NA, nrow = n_days, ncol = n_run)
  for (b in 1:n_run) {
    sim_b <- run_two_tank(behav$k1[b], behav$k2[b], behav$k3[b],
                          times, precip_vec, area_km2 = area_km2)
    Q_ensemble[, b] <- if (!is.null(area_km2)) sim_b$Q_total_m3s else sim_b$Q_total
  }

  Q_lower  <- apply(Q_ensemble, 1, quantile, probs = 0.05)
  Q_upper  <- apply(Q_ensemble, 1, quantile, probs = 0.95)
  Q_median <- apply(Q_ensemble, 1, median)

  containment <- sum(Q_obs >= Q_lower & Q_obs <= Q_upper) / n_days * 100

  if (verbose) {
    cat(sprintf("  ║  Containment (90%%): %.1f%%\n", containment))
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
