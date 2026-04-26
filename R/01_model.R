#' @title Two-Tank Model Core Engine
#' @description Two-tank model with optional NONLINEAR power-law surface runoff.
#'
#' Equations:
#'   Q1 = k1 * S1^b1        (surface runoff; b1=1 is linear, b1>1 is nonlinear)
#'   Q2 = k3 * S2            (linear baseflow)
#'   ET = k4 * S1            (linear evapotranspiration)
#'   dS1/dt = P - k1*S1^b1 - k2*S1 - k4*S1
#'   dS2/dt = k2*S1 - k3*S2
 
two_tank_ode <- function(t, state, pars, precip_fun) {
  S1 <- max(state["S1"], 0)
  S2 <- max(state["S2"], 0)
  P  <- precip_fun(t)
  Q1 <- pars["k1"] * S1
  Q2 <- pars["k3"] * S2
  dS1 <- P - Q1 - pars["k2"] * S1
  dS2 <- pars["k2"] * S1 - Q2
  list(c(dS1 = dS1, dS2 = dS2),
       Q1 = Q1, Q2 = Q2, Q_total = Q1 + Q2, P = P)
}
 
#' Run the Two-Tank Model
#'
#' Simulates daily discharge. Supports both linear (Q1 = k1*S1) and
#' nonlinear (Q1 = k1*S1^b1) surface runoff.
#'
#' @param k1 Numeric. Surface runoff coefficient. In linear mode (b1=1),
#'   Q1 = k1*S1 \[1/day\]. In nonlinear mode (b1>1), Q1 = k1*S1^b1.
#' @param k2 Numeric. Percolation coefficient \[1/day\].
#' @param k3 Numeric. Baseflow coefficient \[1/day\].
#' @param times Numeric vector. Time steps.
#' @param precip_vec Numeric vector. Daily precipitation \[mm/day\].
#' @param S1_0,S2_0 Numeric. Initial storage \[mm\]. Default 0.
#' @param area_km2 Numeric. Catchment area \[km²\]. Default NULL.
#' @param k4 Numeric. ET coefficient \[1/day\]. Default 0.
#' @param b1 Numeric. Power-law exponent for surface runoff. Default 1
#'   (linear). Set b1 > 1 for nonlinear (e.g., 1.5 – 3.0). When b1 > 1,
#'   the model produces sharper peaks because runoff increases
#'   disproportionately with storage.
#'
#' @return Data.frame with time, S1, S2, Q1, Q2, Q_total, P, ET.
#'   If area_km2 is provided, also Q1_m3s, Q2_m3s, Q_total_m3s.
#'
#' @importFrom deSolve ode
#' @export
run_two_tank <- function(k1, k2, k3, times, precip_vec,
                         S1_0 = 0, S2_0 = 0, area_km2 = NULL,
                         k4 = 0, b1 = 1) {
 
  if (!is.numeric(k1) || !is.numeric(k2) || !is.numeric(k3))
    stop("Parameters k1, k2, k3 must be numeric.")
  if (k1 <= 0 || k2 <= 0 || k3 <= 0)
    stop("Parameters must be positive. Got: k1=", k1, " k2=", k2, " k3=", k3)
  if (length(times) != length(precip_vec))
    stop("times and precip_vec must have equal length.")
  if (any(precip_vec < 0)) stop("Precipitation cannot be negative.")
  if (b1 < 0) stop("b1 must be >= 0. Got: ", b1)
 
  n <- length(times)
  S1 <- numeric(n)
  S2 <- numeric(n)
  Q1 <- numeric(n)
  Q2 <- numeric(n)
  ET <- numeric(n)
 
  S1[1] <- S1_0
  S2[1] <- S2_0
  Q1[1] <- k1 * (S1[1] ^ b1)
  Q2[1] <- k3 * S2[1]
  ET[1] <- k4 * S1[1]
 
  for (i in 2:n) {
    dt <- times[i] - times[i - 1]
    P  <- precip_vec[i - 1]
 
    if (b1 == 1) {
      # ── LINEAR MODE: analytical solution (always stable) ──
      k_out  <- k1 + k2 + k4
      S1_eq  <- P / k_out
      S1[i]  <- S1_eq + (S1[i-1] - S1_eq) * exp(-k_out * dt)
    } else {
      # ── NONLINEAR MODE: sub-step Euler integration ──
      # Use smaller sub-steps for stability with power-law
      n_sub <- max(10, ceiling(dt * 10))
      sub_dt <- dt / n_sub
      s1_temp <- S1[i-1]
      for (j in 1:n_sub) {
        q1_temp <- k1 * max(s1_temp, 0)^b1
        perc    <- k2 * s1_temp
        et_temp <- k4 * s1_temp
        ds1     <- P - q1_temp - perc - et_temp
        s1_temp <- max(s1_temp + ds1 * sub_dt, 0)
      }
      S1[i] <- s1_temp
    }
 
    # Surface runoff from end-of-step storage
    Q1[i] <- k1 * max(S1[i], 0)^b1
    ET[i] <- k4 * S1[i]
 
    # Average S1 for percolation
    S1_avg <- (S1[i-1] + S1[i]) / 2
 
    # Lower tank: always linear (physically appropriate for baseflow)
    S2_eq  <- (k2 * S1_avg) / k3
    S2[i]  <- S2_eq + (S2[i-1] - S2_eq) * exp(-k3 * dt)
    Q2[i]  <- k3 * S2[i]
  }
 
  df <- data.frame(
    time = times,
    S1 = S1, S2 = S2,
    Q1 = Q1, Q2 = Q2,
    Q_total = Q1 + Q2,
    P = precip_vec,
    ET = ET
  )
 
  if (!is.null(area_km2)) {
    df$Q1_m3s      <- mmday_to_m3s(df$Q1, area_km2)
    df$Q2_m3s      <- mmday_to_m3s(df$Q2, area_km2)
    df$Q_total_m3s <- mmday_to_m3s(df$Q_total, area_km2)
  }
  return(df)
}
 
 
#' Print Model Summary
#' @param sim Data.frame. Output from \code{\link{run_two_tank}}.
#' @param label Character. Name to display.
#' @param units Character. "mm" or "m3s". Default "mm".
#' @export
print_model_summary <- function(sim, label = "Simulation", units = "mm") {
  if (units == "m3s" && !"Q_total_m3s" %in% names(sim)) {
    warning("No m3s columns found. Using mm/day.")
    units <- "mm"
  }
  Q  <- if (units == "m3s") sim$Q_total_m3s else sim$Q_total
  Q1 <- if (units == "m3s") sim$Q1_m3s      else sim$Q1
  Q2 <- if (units == "m3s") sim$Q2_m3s      else sim$Q2
  u  <- if (units == "m3s") "m³/s" else "mm/day"
 
  cat("\n  ╔══════════════════════════════════════════════════════╗\n")
  cat(sprintf("  ║  %s — Summary\n", label))
  cat("  ╠══════════════════════════════════════════════════════╣\n")
  cat(sprintf("  ║  Duration           : %d days\n", nrow(sim)))
  cat(sprintf("  ║  Total rainfall     : %.1f mm\n", sum(sim$P)))
  cat(sprintf("  ║  Peak discharge     : %.3f %s (day %d)\n",
              max(Q), u, which.max(Q) - 1))
  cat(sprintf("  ║  Peak Q1 (surface)  : %.3f %s\n", max(Q1), u))
  cat(sprintf("  ║  Peak Q2 (baseflow) : %.3f %s\n", max(Q2), u))
  cat(sprintf("  ║  Baseflow index     : %.1f %%\n", sum(Q2)/sum(Q)*100))
  cat(sprintf("  ║  Runoff coefficient : %.3f\n", sum(sim$Q_total)/sum(sim$P)))
  if ("ET" %in% names(sim) && sum(sim$ET) > 0) {
    cat(sprintf("  ║  Total ET           : %.1f mm\n", sum(sim$ET)))
    cat(sprintf("  ║  ET fraction (ET/P) : %.1f %%\n",
                sum(sim$ET) / sum(sim$P) * 100))
  }
  cat("  ╚══════════════════════════════════════════════════════╝\n\n")
}
