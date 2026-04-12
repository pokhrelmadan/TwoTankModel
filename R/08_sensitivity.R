#' @title Sensitivity Analysis
#' @description Local sensitivity analysis by perturbing each parameter
#'   around calibrated values.


#' Sensitivity Analysis
#'
#' Perturbs each parameter by ±\code{perturb} (as a fraction) and records
#' the resulting change in model outputs (peak Q, total volume, BFI, NSE).
#' This is a local one-at-a-time (OAT) sensitivity analysis.
#'
#' @param base_params Named numeric vector. c(k1, k2, k3).
#' @param times Numeric vector. Time steps.
#' @param precip_vec Numeric vector. Daily precipitation \[mm/day\].
#' @param Q_obs Numeric vector. Observed discharge. Optional — if provided,
#'   NSE change is also reported.
#' @param area_km2 Numeric. Catchment area. Default NULL.
#' @param perturb Numeric. Fractional perturbation. Default 0.20 (±20%).
#' @param verbose Logical. Print results. Default TRUE.
#'
#' @return A data.frame with Parameter, Change, Peak_Q_pct, Volume_pct,
#'   BFI_change, NSE_change.
#'
#' @export
sensitivity_analysis <- function(base_params, times, precip_vec,
                                 Q_obs = NULL, area_km2 = NULL,
                                 perturb = 0.20, verbose = TRUE) {

  # Baseline simulation
  base_sim <- run_two_tank(base_params[1], base_params[2], base_params[3],
                           times, precip_vec, area_km2 = area_km2)
  Q_base   <- if (!is.null(area_km2)) base_sim$Q_total_m3s else base_sim$Q_total
  base_peak <- max(Q_base)
  base_vol  <- sum(Q_base)
  base_bfi  <- sum(base_sim$Q2) / sum(base_sim$Q_total) * 100
  base_nse  <- if (!is.null(Q_obs)) calc_nse(Q_obs, Q_base) else NA

  results <- data.frame()
  param_names <- c("k1", "k2", "k3")

  for (i in 1:3) {
    for (direction in c(-1, 1)) {
      p <- base_params
      p[i] <- p[i] * (1 + direction * perturb)

      sim_p <- run_two_tank(p[1], p[2], p[3], times, precip_vec,
                            area_km2 = area_km2)
      Q_p <- if (!is.null(area_km2)) sim_p$Q_total_m3s else sim_p$Q_total

      delta_peak <- (max(Q_p) - base_peak) / base_peak * 100
      delta_vol  <- (sum(Q_p) - base_vol) / base_vol * 100
      delta_bfi  <- (sum(sim_p$Q2) / sum(sim_p$Q_total) * 100) - base_bfi
      delta_nse  <- if (!is.null(Q_obs)) calc_nse(Q_obs, Q_p) - base_nse else NA

      results <- rbind(results, data.frame(
        Parameter   = param_names[i],
        Change      = sprintf("%+.0f%%", direction * perturb * 100),
        Peak_Q_pct  = round(delta_peak, 2),
        Volume_pct  = round(delta_vol, 2),
        BFI_change  = round(delta_bfi, 2),
        NSE_change  = round(delta_nse, 4)
      ))
    }
  }

  if (verbose) {
    cat("\n  ╔══════════════════════════════════════════════════════╗\n")
    cat(sprintf("  ║   SENSITIVITY ANALYSIS (±%.0f%% perturbation)        ║\n",
                perturb * 100))
    cat("  ╠══════════════════════════════════════════════════════╣\n")
    cat("  Baseline: k1=", sprintf("%.4f", base_params[1]),
        " k2=", sprintf("%.4f", base_params[2]),
        " k3=", sprintf("%.4f", base_params[3]), "\n", sep = "")
    cat(sprintf("  Base peak=%.3f, volume=%.1f, BFI=%.1f%%",
                base_peak, base_vol, base_bfi))
    if (!is.null(Q_obs)) cat(sprintf(", NSE=%.4f", base_nse))
    cat("\n\n")
    print(results, row.names = FALSE)
    cat("\n")

    # Rank parameters by sensitivity
    avg_sens <- aggregate(abs(results$Peak_Q_pct),
                          list(Parameter = results$Parameter), mean)
    avg_sens <- avg_sens[order(-avg_sens$x), ]
    cat("  Ranking by peak-Q sensitivity (most → least):\n")
    for (i in 1:nrow(avg_sens))
      cat(sprintf("    %d. %s (avg |ΔPeak| = %.1f%%)\n",
                  i, avg_sens$Parameter[i], avg_sens$x[i]))
    cat("\n")
  }

  return(results)
}
