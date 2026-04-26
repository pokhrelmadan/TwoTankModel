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

  # ── Adaptive threshold fallback ──
  # If no sets pass the threshold, use the top 5% of samples instead
  # so the user still gets uncertainty bounds (just wider ones).
  adaptive_threshold <- nse_threshold
  used_adaptive <- FALSE
  if (n_behav == 0) {
    used_adaptive <- TRUE
    cutoff_quantile <- quantile(samples$NSE, 0.95, na.rm = TRUE)
    behav <- samples[samples$NSE >= cutoff_quantile, ]
    behav <- behav[order(-behav$NSE), ]
    n_behav <- nrow(behav)
    adaptive_threshold <- cutoff_quantile

    if (verbose) {
      cat("\n  ⚠ WARNING: No samples passed NSE >= ",
          sprintf("%.2f", nse_threshold), "\n", sep = "")
      cat("    The calibrated model has poor fit (best NSE = ",
          sprintf("%.4f", max(samples$NSE, na.rm = TRUE)), ")\n", sep = "")
      cat("    Using top 5%% of samples for uncertainty bounds instead\n")
      cat("    (adaptive threshold: NSE >= ",
          sprintf("%.4f", cutoff_quantile), ")\n\n", sep = "")
      cat("    Possible issues to check:\n")
      cat("    - Catchment area (unit conversion depends on it)\n")
      cat("    - Date alignment between rainfall and discharge files\n")
      cat("    - Missing processes (evapotranspiration, snowmelt)\n")
      cat("    - Parameter bounds may be too restrictive\n\n")
    }
  }

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
    if (used_adaptive) {
      cat(sprintf("  ║  Threshold      : adaptive (top 5%%, NSE >= %.4f)\n",
                  adaptive_threshold))
    } else {
      cat(sprintf("  ║  Threshold      : NSE >= %.2f\n", nse_threshold))
    }
    cat(sprintf("  ║  Behavioural    : %d of %d (%s%%)\n",
                n_behav, nrow(samples), pct))
    cat("  ║  Parameter uncertainty ranges:\n")
    for (i in 1:3)
      cat(sprintf("  ║    %-22s: [%.4f – %.4f]\n",
                  param_ranges$parameter[i], param_ranges$min[i],
                  param_ranges$max[i]))
  }

  # Check if k4 (ET) and b1 (nonlinear) were calibrated
  has_k4 <- "k4" %in% names(behav)
  has_b1 <- "b1" %in% names(behav) && any(behav$b1 != 1)

  n_run <- min(n_behav, max_behavioural)
  if (verbose) cat(sprintf("  ║  Computing envelope (%d sims)...\n", n_run))

  Q_ensemble <- matrix(NA, nrow = n_days, ncol = n_run)

  # Use parallel if we have enough sims and multiple cores
  n_cores <- max(1, parallel::detectCores(logical = FALSE) - 1)
  if (n_cores > 1 && n_run >= 50) {
    param_list <- lapply(1:n_run, function(b) {
      k4_val <- if (has_k4) behav$k4[b] else 0
      b1_val <- if (has_b1) behav$b1[b] else 1
      c(behav$k1[b], behav$k2[b], behav$k3[b], k4_val, b1_val)
    })

    worker_fun <- function(p, times, precip_vec, area_km2) {
      sim_b <- run_two_tank(p[1], p[2], p[3], times, precip_vec,
                            area_km2 = area_km2, k4 = p[4], b1 = p[5])
      if (!is.null(area_km2)) sim_b$Q_total_m3s else sim_b$Q_total
    }

    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterExport(cl,
        c("run_two_tank", "two_tank_ode", "mmday_to_m3s"),
        envir = environment())
      parallel::clusterEvalQ(cl, library(deSolve))
      results <- parallel::parLapply(cl, param_list, worker_fun,
                                     times = times, precip_vec = precip_vec,
                                     area_km2 = area_km2)
    } else {
      results <- parallel::mclapply(param_list, worker_fun,
                                    times = times, precip_vec = precip_vec,
                                    area_km2 = area_km2,
                                    mc.cores = n_cores)
    }
    for (b in 1:n_run) Q_ensemble[, b] <- results[[b]]
  } else {
    for (b in 1:n_run) {
      k4_val <- if (has_k4) behav$k4[b] else 0
      b1_val <- if (has_b1) behav$b1[b] else 1
      sim_b <- run_two_tank(behav$k1[b], behav$k2[b], behav$k3[b],
                            times, precip_vec, area_km2 = area_km2,
                            k4 = k4_val, b1 = b1_val)
      Q_ensemble[, b] <- if (!is.null(area_km2)) sim_b$Q_total_m3s else sim_b$Q_total
    }
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
