#' @title Sensitivity Analysis
#' @description Local sensitivity analysis by perturbing each parameter
#'   around calibrated values. Supports k4 (ET) and b1 (nonlinear exponent).

#' Sensitivity Analysis
#'
#' @param base_params Named numeric vector from calibration (e.g., c(k1, k2, k3, k4, b1)).
#' @param times Numeric vector. Time steps.
#' @param precip_vec Numeric vector. Daily precipitation.
#' @param Q_obs Numeric vector. Observed discharge. Optional.
#' @param area_km2 Numeric. Catchment area. Default NULL.
#' @param perturb Numeric. Fractional perturbation. Default 0.20.
#' @param verbose Logical. Default TRUE.
#'
#' @return A data.frame with sensitivity results.
#' @export
sensitivity_analysis <- function(base_params, times, precip_vec,
                                 Q_obs = NULL, area_km2 = NULL,
                                 perturb = 0.20, verbose = TRUE) {

  # Detect which parameters were calibrated
  pnames <- names(base_params)
  has_k4 <- "k4" %in% pnames
  has_b1 <- "b1" %in% pnames

  base_k4 <- if (has_k4) base_params["k4"] else 0
  base_b1 <- if (has_b1) base_params["b1"] else 1

  # Baseline simulation
  base_sim <- run_two_tank(base_params["k1"], base_params["k2"], base_params["k3"],
                           times, precip_vec, area_km2 = area_km2,
                           k4 = base_k4, b1 = base_b1)
  Q_base <- if (!is.null(area_km2)) base_sim$Q_total_m3s else base_sim$Q_total

  if (any(is.na(Q_base))) {
    stop("Baseline simulation failed. The calibrated parameters may be ",
         "numerically unstable. Try different bounds in config.R.")
  }

  base_peak <- max(Q_base)
  base_vol  <- sum(Q_base)
  base_bfi  <- sum(base_sim$Q2) / sum(base_sim$Q_total) * 100
  base_nse  <- if (!is.null(Q_obs)) calc_nse(Q_obs, Q_base) else NA

  # Build list of parameters to perturb
  param_names <- c("k1", "k2", "k3")
  if (has_k4) param_names <- c(param_names, "k4")
  if (has_b1) param_names <- c(param_names, "b1")

  results <- data.frame()

  for (pname in param_names) {
    for (direction in c(-1, 1)) {
      k1_p <- base_params["k1"]
      k2_p <- base_params["k2"]
      k3_p <- base_params["k3"]
      k4_p <- base_k4
      b1_p <- base_b1

      if (pname == "k1") k1_p <- k1_p * (1 + direction * perturb)
      if (pname == "k2") k2_p <- k2_p * (1 + direction * perturb)
      if (pname == "k3") k3_p <- k3_p * (1 + direction * perturb)
      if (pname == "k4") k4_p <- k4_p * (1 + direction * perturb)
      if (pname == "b1") b1_p <- b1_p * (1 + direction * perturb)

      sim_p <- run_two_tank(k1_p, k2_p, k3_p, times, precip_vec,
                            area_km2 = area_km2, k4 = k4_p, b1 = b1_p)
      Q_p <- if (!is.null(area_km2)) sim_p$Q_total_m3s else sim_p$Q_total

      if (any(is.na(Q_p))) {
        results <- rbind(results, data.frame(
          Parameter = pname, Change = sprintf("%+.0f%%", direction * perturb * 100),
          Peak_Q_pct = NA, Volume_pct = NA, BFI_change = NA, NSE_change = NA))
        next
      }

      delta_peak <- (max(Q_p) - base_peak) / base_peak * 100
      delta_vol  <- (sum(Q_p) - base_vol) / base_vol * 100
      delta_bfi  <- (sum(sim_p$Q2) / sum(sim_p$Q_total) * 100) - base_bfi
      delta_nse  <- if (!is.null(Q_obs)) calc_nse(Q_obs, Q_p) - base_nse else NA

      results <- rbind(results, data.frame(
        Parameter   = pname,
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
                        base_params["k1"], base_params["k2"], base_params["k3"])
    if (has_k4) base_str <- paste0(base_str, sprintf(" k4=%.4f", base_k4))
    if (has_b1) base_str <- paste0(base_str, sprintf(" b1=%.4f", base_b1))
    cat(base_str, "\n", sep = "")
    cat(sprintf("  Base peak=%.3f, volume=%.1f, BFI=%.1f%%",
                base_peak, base_vol, base_bfi))
    if (!is.null(Q_obs)) cat(sprintf(", NSE=%.4f", base_nse))
    cat("\n\n")
    print(results, row.names = FALSE)
    cat("\n")

    avg_sens <- aggregate(abs(results$Peak_Q_pct),
                          list(Parameter = results$Parameter), mean, na.rm = TRUE)
    avg_sens <- avg_sens[order(-avg_sens$x), ]
    cat("  Ranking by peak-Q sensitivity (most → least):\n")
    for (i in 1:nrow(avg_sens))
      cat(sprintf("    %d. %s (avg |ΔPeak| = %.1f%%)\n",
                  i, avg_sens$Parameter[i], avg_sens$x[i]))
    cat("\n")
  }

  return(results)
}
