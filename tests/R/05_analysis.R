#' @title Hydrograph Analysis
#' @description Extract hydrograph metrics, monthly water balance,
#'   mass balance check, and watershed comparison.


#' Extract Hydrograph Metrics
#'
#' Computes key hydrological response characteristics from a simulation.
#'
#' @param sim Data.frame. Output from \code{\link{run_two_tank}}.
#' @param dates Date vector. Optional. If provided, peak is reported
#'   as a date instead of day index. Default NULL.
#' @param dt Numeric. Time step in days. Default 1.
#' @param label Character. Watershed name. Default "Watershed".
#' @param verbose Logical. Print results. Default TRUE.
#'
#' @return A data.frame with one row of metrics.
#' @export
extract_metrics <- function(sim, dates = NULL, dt = 1,
                            label = "Watershed", verbose = TRUE) {

  peak_Q     <- max(sim$Q_total)
  peak_idx   <- which.max(sim$Q_total)
  runoff_vol <- sum(sim$Q_total) * dt
  peak_Q1    <- max(sim$Q1)
  peak_Q2    <- max(sim$Q2)
  bfi        <- sum(sim$Q2) / sum(sim$Q_total) * 100
  rc         <- runoff_vol / sum(sim$P)

  if (!is.null(dates)) {
    peak_time <- as.character(dates[peak_idx])
  } else {
    peak_time <- as.character(sim$time[peak_idx])
  }

  if (verbose) {
    cat(sprintf("\n  ── %s ────────────────────────────────\n", label))
    cat(sprintf("    Peak discharge (Q)  : %.2f mm/day\n", peak_Q))
    cat(sprintf("    Peak date/time      : %s\n", peak_time))
    cat(sprintf("    Total runoff        : %.1f mm\n", runoff_vol))
    cat(sprintf("    Peak surface (Q1)   : %.2f mm/day\n", peak_Q1))
    cat(sprintf("    Peak baseflow (Q2)  : %.2f mm/day\n", peak_Q2))
    cat(sprintf("    Baseflow index      : %.1f %%\n", bfi))
    cat(sprintf("    Runoff coefficient  : %.3f\n\n", rc))
  }

  data.frame(
    Label = label, Peak_Q = round(peak_Q, 3),
    Peak_Time = peak_time, Runoff_Volume = round(runoff_vol, 1),
    Peak_Q1 = round(peak_Q1, 3), Peak_Q2 = round(peak_Q2, 3),
    BFI = round(bfi, 1), Runoff_Coeff = round(rc, 3),
    stringsAsFactors = FALSE
  )
}


#' Monthly Water Balance Summary
#'
#' @param sim Data.frame. Output from \code{\link{run_two_tank}}.
#' @param dates Date vector. Corresponding dates.
#' @param verbose Logical. Print table. Default TRUE.
#' @return A data.frame with monthly totals.
#' @export
monthly_summary <- function(sim, dates, verbose = TRUE) {
  months_lab <- format(dates, "%Y-%m")
  df  <- data.frame(month = months_lab, P = sim$P, Q = sim$Q_total,
                    Q1 = sim$Q1, Q2 = sim$Q2)
  agg <- aggregate(. ~ month, data = df, FUN = sum)
  agg$BFI <- round(agg$Q2 / agg$Q * 100, 1)
  agg$RC  <- round(agg$Q / agg$P, 3)
  agg[, 2:5] <- round(agg[, 2:5], 1)
  if (verbose) { cat("\n  Monthly water balance:\n"); print(agg, row.names = FALSE); cat("\n") }
  return(agg)
}


#' Check Mass Balance
#'
#' Verifies P = Q + delta_S. Residual should be near zero.
#'
#' @param sim Data.frame. Output from \code{\link{run_two_tank}}.
#' @param dt Numeric. Time step in days. Default 1.
#' @param verbose Logical. Print results. Default TRUE.
#' @return A named list with P_total, Q_total, delta_S, residual.
#' @export
check_mass_balance <- function(sim, dt = 1, verbose = TRUE) {
  P_total  <- sum(sim$P) * dt
  Q_total  <- sum(sim$Q_total) * dt
  delta_S  <- (tail(sim$S1, 1) + tail(sim$S2, 1)) - (sim$S1[1] + sim$S2[1])
  residual <- P_total - Q_total - delta_S

  if (verbose) {
    cat("\n  ┌──────────────────────────────────────┐\n")
    cat("  │        MASS BALANCE CHECK            │\n")
    cat("  ├──────────────────────────────────────┤\n")
    cat(sprintf("  │  Precipitation (P) : %8.1f mm     │\n", P_total))
    cat(sprintf("  │  Discharge (Q)     : %8.1f mm     │\n", Q_total))
    cat(sprintf("  │  Storage change    : %8.2f mm     │\n", delta_S))
    cat(sprintf("  │  Residual (P-Q-ΔS) : %8.6f mm   │\n", residual))
    if (abs(residual) < 0.01) {
      cat("  │  Status: ✓ Closed                   │\n")
    } else {
      cat("  │  Status: ✗ Imbalance detected        │\n")
    }
    cat("  └──────────────────────────────────────┘\n\n")
  }
  return(list(P_total = P_total, Q_total = Q_total,
              delta_S = delta_S, residual = residual))
}


#' Compare Two Watersheds
#'
#' Runs the model with two parameter sets under the same rainfall
#' and prints a side-by-side comparison.
#'
#' @param times Numeric vector. Time steps.
#' @param precip_vec Numeric vector. Daily precipitation.
#' @param params_A Named numeric vector c(k1, k2, k3) for watershed A.
#' @param params_B Named numeric vector c(k1, k2, k3) for watershed B.
#' @param label_A Character. Name for A. Default "Watershed A".
#' @param label_B Character. Name for B. Default "Watershed B".
#' @param dates Date vector. Optional. Default NULL.
#'
#' @return A list with sim_A, sim_B, metrics_A, metrics_B, comparison.
#' @export
compare_watersheds <- function(times, precip_vec,
                               params_A, params_B,
                               label_A = "Watershed A",
                               label_B = "Watershed B",
                               dates = NULL) {

  cat("\n")
  cat("  ╔══════════════════════════════════════════════════════╗\n")
  cat("  ║           WATERSHED COMPARISON                      ║\n")
  cat("  ╠══════════════════════════════════════════════════════╣\n")
  cat(sprintf("  ║  %-15s: k1=%.3f  k2=%.3f  k3=%.3f\n",
              label_A, params_A[1], params_A[2], params_A[3]))
  cat(sprintf("  ║  %-15s: k1=%.3f  k2=%.3f  k3=%.3f\n",
              label_B, params_B[1], params_B[2], params_B[3]))
  cat("  ║  Same precipitation forcing applied to both.\n")
  cat("  ╚══════════════════════════════════════════════════════╝\n")

  sim_A <- run_two_tank(params_A[1], params_A[2], params_A[3],
                        times, precip_vec)
  sim_B <- run_two_tank(params_B[1], params_B[2], params_B[3],
                        times, precip_vec)

  met_A <- extract_metrics(sim_A, dates, label = label_A)
  met_B <- extract_metrics(sim_B, dates, label = label_B)
  comparison <- rbind(met_A, met_B)

  cat("  Side-by-side comparison:\n")
  print(comparison, row.names = FALSE)
  cat("\n")

  return(list(sim_A = sim_A, sim_B = sim_B,
              metrics_A = met_A, metrics_B = met_B,
              comparison = comparison))
}
