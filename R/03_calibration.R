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
#' @param k1_range,k2_range,k3_range,k4_range Numeric(2). Parameter bounds.
#'   Set k4_range = c(0, 0) to disable ET. Default k4_range = c(0, 0.10).
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
                                 k4_range = c(0, 0.10),
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
 
  # Check if ET calibration is enabled
  use_k4 <- (k4_range[2] > 0)
  n_par  <- if (use_k4) 4 else 3
 
  par_lower <- c(k1_range[1], k2_range[1], k3_range[1])
  par_upper <- c(k1_range[2], k2_range[2], k3_range[2])
  if (use_k4) {
    par_lower <- c(par_lower, k4_range[1])
    par_upper <- c(par_upper, k4_range[2])
  }
 
  if (verbose) {
    cat("\n  ╔══════════════════════════════════════════════════════╗\n")
    cat("  ║          MONTE CARLO CALIBRATION (LHS)              ║\n")
    cat("  ╠══════════════════════════════════════════════════════╣\n")
    cat(sprintf("  ║  Samples        : %d\n", n_samples))
    cat(sprintf("  ║  Objective      : %s\n", objective))
    cat(sprintf("  ║  Observed units : %s\n", obs_units))
    if (!is.null(area_km2))
      cat(sprintf("  ║  Catchment area : %.1f km²\n", area_km2))
    cat(sprintf("  ║  k1 range       : [%.3f – %.3f]  (surface runoff)\n",
                k1_range[1], k1_range[2]))
    cat(sprintf("  ║  k2 range       : [%.3f – %.3f]  (percolation)\n",
                k2_range[1], k2_range[2]))
    cat(sprintf("  ║  k3 range       : [%.3f – %.3f]  (baseflow)\n",
                k3_range[1], k3_range[2]))
    if (use_k4) {
      cat(sprintf("  ║  k4 range       : [%.3f – %.3f]  (evapotranspiration)\n",
                  k4_range[1], k4_range[2]))
    } else {
      cat("  ║  k4 (ET)        : disabled (k4_range[2] = 0)\n")
    }
    cat("  ╚══════════════════════════════════════════════════════╝\n\n")
  }
 
  set.seed(seed)
  lhs_mat <- lhs::randomLHS(n_samples, n_par)
 
  samples <- data.frame(
    k1 = par_lower[1] + lhs_mat[,1] * (par_upper[1] - par_lower[1]),
    k2 = par_lower[2] + lhs_mat[,2] * (par_upper[2] - par_lower[2]),
    k3 = par_lower[3] + lhs_mat[,3] * (par_upper[3] - par_lower[3])
  )
  if (use_k4) {
    samples$k4 <- par_lower[4] + lhs_mat[,4] * (par_upper[4] - par_lower[4])
  } else {
    samples$k4 <- 0
  }
  samples$NSE    <- NA_real_
  samples$KGE    <- NA_real_
  samples$LogNSE <- NA_real_
  samples$RMSE   <- NA_real_
  samples$PBIAS  <- NA_real_
  samples$Peak_Q <- NA_real_
  samples$Vol_Q  <- NA_real_
 
  t0 <- proc.time()
  best_obj <- -Inf
  best_idx <- 1
 
  # ── Detect available cores for parallel processing ──
  n_cores <- max(1, parallel::detectCores(logical = FALSE) - 1)
  use_parallel <- n_cores > 1 && n_samples >= 100
 
  if (verbose) {
    if (use_parallel) {
      cat(sprintf("  Running %d simulations on %d cores (parallel)...\n",
                  n_samples, n_cores))
    } else {
      cat(sprintf("  Running %d simulations (sequential)...\n", n_samples))
    }
    flush.console()
  }
 
  if (use_parallel) {
    # ── PARALLEL execution ──
    # Build list of parameter sets
    param_list <- lapply(1:n_samples, function(i) {
      c(samples$k1[i], samples$k2[i], samples$k3[i], samples$k4[i])
    })
 
    # Worker function — must be self-contained for parallel workers
    worker_fun <- function(p, times, precip_vec, Q_obs, area_km2) {
      sim_i <- run_two_tank(p[1], p[2], p[3], times, precip_vec,
                            area_km2 = area_km2, k4 = p[4])
      Q_sim <- if (!is.null(area_km2)) sim_i$Q_total_m3s else sim_i$Q_total
 
      # Handle failed simulations
      if (any(is.na(Q_sim))) {
        return(c(NSE = -Inf, KGE = -Inf, LogNSE = -Inf,
                 RMSE = Inf, PBIAS = NA, Peak_Q = NA, Vol_Q = NA))
      }
 
      c(NSE    = calc_nse(Q_obs, Q_sim),
        KGE    = calc_kge(Q_obs, Q_sim),
        LogNSE = calc_lognse(Q_obs, Q_sim),
        RMSE   = calc_rmse(Q_obs, Q_sim),
        PBIAS  = calc_pbias(Q_obs, Q_sim),
        Peak_Q = max(Q_sim),
        Vol_Q  = sum(Q_sim))
    }
 
    # Run in chunks to show progress
    chunk_size <- max(1, floor(n_samples / 20))  # 20 progress updates
    chunks <- split(1:n_samples,
                    ceiling(seq_along(1:n_samples) / chunk_size))
 
    if (.Platform$OS.type == "windows") {
      # Windows: use PSOCK cluster
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
 
      # Export necessary objects and functions to workers
      parallel::clusterExport(cl,
        c("run_two_tank", "two_tank_ode", "calc_nse", "calc_kge",
          "calc_lognse", "calc_rmse", "calc_pbias", "mmday_to_m3s"),
        envir = environment())
      parallel::clusterEvalQ(cl, library(deSolve))
 
      for (ci in seq_along(chunks)) {
        chunk_idx <- chunks[[ci]]
        results <- parallel::parLapply(cl, param_list[chunk_idx],
                                       worker_fun,
                                       times = times,
                                       precip_vec = precip_vec,
                                       Q_obs = Q_obs,
                                       area_km2 = area_km2)
        # Store results
        for (j in seq_along(chunk_idx)) {
          i <- chunk_idx[j]
          samples[i, c("NSE", "KGE", "LogNSE", "RMSE", "PBIAS",
                       "Peak_Q", "Vol_Q")] <- results[[j]]
        }
 
        # Update progress + best
        completed <- tail(chunk_idx, 1)
        current_best_in_all <- switch(objective,
          "NSE"    = max(samples$NSE[1:completed], na.rm = TRUE),
          "KGE"    = max(samples$KGE[1:completed], na.rm = TRUE),
          "LogNSE" = max(samples$LogNSE[1:completed], na.rm = TRUE))
 
        if (verbose) {
          pct <- round(completed / n_samples * 100)
          bar <- round(completed / n_samples * 40)
          cat(sprintf("\r  [%s%s] %3d%%  (best %s: %.4f)",
                      paste(rep("█", bar), collapse = ""),
                      paste(rep(" ", 40 - bar), collapse = ""),
                      pct, objective, current_best_in_all))
          flush.console()
        }
      }
 
    } else {
      # Mac/Linux: use mclapply (faster, fork-based)
      for (ci in seq_along(chunks)) {
        chunk_idx <- chunks[[ci]]
        results <- parallel::mclapply(param_list[chunk_idx],
                                      worker_fun,
                                      times = times,
                                      precip_vec = precip_vec,
                                      Q_obs = Q_obs,
                                      area_km2 = area_km2,
                                      mc.cores = n_cores)
        for (j in seq_along(chunk_idx)) {
          i <- chunk_idx[j]
          samples[i, c("NSE", "KGE", "LogNSE", "RMSE", "PBIAS",
                       "Peak_Q", "Vol_Q")] <- results[[j]]
        }
 
        completed <- tail(chunk_idx, 1)
        current_best_in_all <- switch(objective,
          "NSE"    = max(samples$NSE[1:completed], na.rm = TRUE),
          "KGE"    = max(samples$KGE[1:completed], na.rm = TRUE),
          "LogNSE" = max(samples$LogNSE[1:completed], na.rm = TRUE))
 
        if (verbose) {
          pct <- round(completed / n_samples * 100)
          bar <- round(completed / n_samples * 40)
          cat(sprintf("\r  [%s%s] %3d%%  (best %s: %.4f)",
                      paste(rep("█", bar), collapse = ""),
                      paste(rep(" ", 40 - bar), collapse = ""),
                      pct, objective, current_best_in_all))
          flush.console()
        }
      }
    }
 
    # Find overall best
    best_idx <- switch(objective,
      "NSE"    = which.max(samples$NSE),
      "KGE"    = which.max(samples$KGE),
      "LogNSE" = which.max(samples$LogNSE))
 
  } else {
    # ── SEQUENTIAL execution (fallback) ──
    update_every <- max(1, floor(n_samples / 40))  # update every 2.5%
 
    for (i in 1:n_samples) {
      sim_i <- run_two_tank(samples$k1[i], samples$k2[i], samples$k3[i],
                            times, precip_vec, area_km2 = area_km2,
                            k4 = samples$k4[i])
      Q_sim <- if (!is.null(area_km2)) sim_i$Q_total_m3s else sim_i$Q_total
 
      # Skip failed simulations
      if (any(is.na(Q_sim))) {
        samples$NSE[i]    <- -Inf
        samples$KGE[i]    <- -Inf
        samples$LogNSE[i] <- -Inf
        samples$RMSE[i]   <- Inf
        samples$PBIAS[i]  <- NA
        samples$Peak_Q[i] <- NA
        samples$Vol_Q[i]  <- NA
        next
      }
 
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
 
      if (verbose && (i %% update_every == 0 || i == n_samples)) {
        pct <- round(i / n_samples * 100)
        bar <- round(i / n_samples * 40)
        cat(sprintf("\r  [%s%s] %3d%%  (best %s: %.4f)",
                    paste(rep("█", bar), collapse = ""),
                    paste(rep(" ", 40 - bar), collapse = ""),
                    pct, objective, best_obj))
        flush.console()
      }
    }
  }
 
  elapsed <- (proc.time() - t0)[3]
  if (verbose) cat("\n\n")
 
  if (use_k4) {
    best_params <- unlist(samples[best_idx, c("k1", "k2", "k3", "k4")])
    names(best_params) <- c("k1", "k2", "k3", "k4")
  } else {
    best_params <- unlist(samples[best_idx, c("k1", "k2", "k3")])
    names(best_params) <- c("k1", "k2", "k3")
  }
  best_k4 <- if (use_k4) best_params["k4"] else 0
  best_sim <- run_two_tank(best_params[1], best_params[2], best_params[3],
                           times, precip_vec, area_km2 = area_km2,
                           k4 = best_k4)
 
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
    if (use_k4) {
      cat(sprintf("  │  k4 = %.4f  (evapotranspiration)                  │\n", best_params[4]))
    }
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
