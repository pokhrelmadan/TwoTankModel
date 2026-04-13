#' @title Two-Tank Model Core Engine
 
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
#' Simulates daily discharge using the two-tank linear reservoir model.
#' When \code{area_km2} is provided, the output includes m³/s columns.
#'
#' @param k1,k2,k3 Numeric. Parameters \[1/day\].
#' @param times Numeric vector. Time steps.
#' @param precip_vec Numeric vector. Daily precipitation \[mm/day\].
#' @param S1_0,S2_0 Numeric. Initial tank storage \[mm\]. Default 0.
#' @param area_km2 Numeric. Catchment area in km². Default NULL.
#'
#' @return Data.frame with simulation results in mm/day. If
#'   \code{area_km2} is provided, also includes Q1_m3s, Q2_m3s, Q_total_m3s.
#'
#' @importFrom deSolve ode
#' @export
run_two_tank <- function(k1, k2, k3, times, precip_vec,
                         S1_0 = 0, S2_0 = 0, area_km2 = NULL) {
 
  if (!is.numeric(k1) || !is.numeric(k2) || !is.numeric(k3))
    stop("Parameters k1, k2, k3 must be numeric.")
  if (k1 <= 0 || k2 <= 0 || k3 <= 0)
    stop("Parameters must be positive. Got: k1=", k1, " k2=", k2, " k3=", k3)
  if (length(times) != length(precip_vec))
    stop("times and precip_vec must have equal length.")
  if (any(precip_vec < 0)) stop("Precipitation cannot be negative.")
 
  precip_fun <- approxfun(times, precip_vec, method = "constant", rule = 2)
  pars  <- c(k1 = k1, k2 = k2, k3 = k3)
  state <- c(S1 = S1_0, S2 = S2_0)
 
  out <- tryCatch(
    deSolve::ode(y = state, times = times, func = two_tank_ode,
                 parms = pars, precip_fun = precip_fun,
                 method = "lsoda",
                 rtol = 1e-4, atol = 1e-6,
                 hmax = 1),
    error   = function(e) NULL,
    warning = function(w) {
      # Suppress DLSODA precision warnings — retry with looser tolerance
      if (grepl("accuracy|precision|early", w$message, ignore.case = TRUE)) {
        tryCatch(
          deSolve::ode(y = state, times = times, func = two_tank_ode,
                       parms = pars, precip_fun = precip_fun,
                       method = "euler", hini = 0.5),
          error = function(e) NULL)
      } else {
        deSolve::ode(y = state, times = times, func = two_tank_ode,
                     parms = pars, precip_fun = precip_fun,
                     method = "lsoda")
      }
    }
  )
 
  # If solver failed or returned incomplete output, return NA-filled data
  if (is.null(out) || nrow(out) != length(times)) {
    df <- data.frame(time = times, S1 = NA, S2 = NA,
                     Q1 = NA, Q2 = NA, Q_total = NA,
                     P = precip_vec)
    if (!is.null(area_km2)) {
      df$Q1_m3s <- NA; df$Q2_m3s <- NA; df$Q_total_m3s <- NA
    }
    return(df)
  }
 
  df <- as.data.frame(out)
  names(df) <- c("time", "S1", "S2", "Q1", "Q2", "Q_total", "P")
 
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
  cat("  ╚══════════════════════════════════════════════════════╝\n\n")
}
