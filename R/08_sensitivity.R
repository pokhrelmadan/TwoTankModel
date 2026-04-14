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
 
  # Check if k4 was calibrated
  has_k4 <- length(base_params) >= 4 && !is.na(base_params[4])
  base_k4 <- if (has_k4) base_params[4] else 0
 
  # Baseline simulation
  base_sim <- run_two_tank(base_params[1], base_params[2], base_params[3],
                           times, precip_vec, area_km2 = area_km2,
                           k4 = base_k4)
  Q_base   <- if (!is.null(area_km2)) base_sim$Q_total_m3s else base_sim$Q_total
 
  if (any(is.na(Q_base))) {
    stop("Baseline simulation failed. The calibrated parameters may be ",
         "numerically unstable. Try different bounds in config.R.")
  }
 
  base_peak <- max(Q_base)
  base_vol  <- sum(Q_base)
  base_bfi  <- sum(base_sim$Q2) / sum(base_sim$Q_total) * 100
  base_nse  <- if (!is.null(Q_obs)) calc_nse(Q_obs, Q_base) else NA
 
  results <- data.frame()
  n_params <- if (has_k4) 4 else 3
  param_names <- c("k1", "k2", "k3", "k4")[1:n_params]
 
  for (i in 1:n_params) {
    for (direction in c(-1, 1)) {
      p <- base_params
      k4_p <- base_k4
      if (i <= 3) {
        p[i] <- p[i] * (1 + direction * perturb)
      } else {
        k4_p <- base_k4 * (1 + direction * perturb)
      }
 
      sim_p <- run_two_tank(p[1], p[2], p[3], times, precip_vec,
                            area_km2 = area_km2, k4 = k4_p)
      Q_p <- if (!is.null(area_km2)) sim_p$Q_total_m3s else sim_p$Q_total
 
      if (any(is.na(Q_p))) {
        results <- rbind(results, data.frame(
          Parameter   = param_names[i],
          Change      = sprintf("%+.0f%%", direction * perturb * 100),
          Peak_Q_pct  = NA, Volume_pct = NA,
          BFI_change  = NA, NSE_change = NA))
        next
      }
 
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
    base_str <- sprintf("  Baseline: k1=%.4f k2=%.4f k3=%.4f",
                        base_params[1], base_params[2], base_params[3])
    if (has_k4) base_str <- paste0(base_str, sprintf(" k4=%.4f", base_k4))
    cat(base_str, "\n", sep = "")
    cat(sprintf("  Base peak=%.3f, volume=%.1f, BFI=%.1f%%",
                base_peak, base_vol, base_bfi))
    if (!is.null(Q_obs)) cat(sprintf(", NSE=%.4f", base_nse))
    cat("\n\n")
    print(results, row.names = FALSE)
    cat("\n")
 
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
