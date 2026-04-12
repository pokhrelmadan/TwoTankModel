#' @title Two-Tank Model Core Engine
#' @description Core ODE system and simulation runner for the two-tank
#'   conceptual rainfall-runoff model.
#'
#' @details
#' The model consists of two vertically connected linear reservoirs:
#' \itemize{
#'   \item Upper tank: surface and interflow (drains via k1, percolates via k2)
#'   \item Lower tank: baseflow (drains via k3)
#' }
#'
#' Governing equations:
#' \deqn{dS1/dt = P(t) - k1 \cdot S1(t) - k2 \cdot S1(t)}
#' \deqn{dS2/dt = k2 \cdot S1(t) - k3 \cdot S2(t)}
#' \deqn{Q(t) = k1 \cdot S1(t) + k3 \cdot S2(t)}

# ── Internal ODE function (not exported) ────────────────────────────────────

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
#' This is the main simulation function — give it parameters and rainfall,
#' and it returns a complete daily time series.
#'
#' @param k1 Numeric. Surface runoff coefficient \[1/day\].
#'   Higher k1 = faster surface response (e.g., urban: 0.3, forest: 0.08).
#' @param k2 Numeric. Percolation coefficient \[1/day\].
#'   Higher k2 = more water moves to the lower tank.
#' @param k3 Numeric. Baseflow coefficient \[1/day\].
#'   Higher k3 = faster baseflow recession.
#' @param times Numeric vector. Time steps (typically 0, 1, 2, ... n_days-1).
#' @param precip_vec Numeric vector. Daily precipitation \[mm/day\],
#'   same length as \code{times}.
#' @param S1_0 Numeric. Initial storage in upper tank \[mm\]. Default 0.
#' @param S2_0 Numeric. Initial storage in lower tank \[mm\]. Default 0.
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{time}{Time step index (0, 1, 2, ...)}
#'     \item{S1}{Upper tank storage \[mm\]}
#'     \item{S2}{Lower tank storage \[mm\]}
#'     \item{Q1}{Surface runoff \[mm/day\]}
#'     \item{Q2}{Baseflow \[mm/day\]}
#'     \item{Q_total}{Total discharge Q1 + Q2 \[mm/day\]}
#'     \item{P}{Precipitation at each time step \[mm/day\]}
#'   }
#'
#' @examples
#' # Simple 1-year simulation
#' times  <- 0:364
#' precip <- c(rep(0, 50), rep(10, 20), rep(0, 295))
#' sim    <- run_two_tank(k1 = 0.3, k2 = 0.1, k3 = 0.02, times, precip)
#' plot(sim$time, sim$Q_total, type = "l", xlab = "Day", ylab = "Q (mm/day)")
#'
#' @importFrom deSolve ode
#' @export
run_two_tank <- function(k1, k2, k3, times, precip_vec,
                         S1_0 = 0, S2_0 = 0) {

  # ── Friendly input validation ──
  if (!is.numeric(k1) || !is.numeric(k2) || !is.numeric(k3)) {
    stop("Parameters k1, k2, k3 must be numeric values.")
  }
  if (k1 <= 0 || k2 <= 0 || k3 <= 0) {
    stop("Parameters k1, k2, k3 must be positive (> 0).\n",
         "  You provided: k1=", k1, "  k2=", k2, "  k3=", k3)
  }
  if (length(times) != length(precip_vec)) {
    stop("'times' and 'precip_vec' must have the same length.\n",
         "  times has ", length(times), " elements, ",
         "precip_vec has ", length(precip_vec), " elements.")
  }
  if (any(precip_vec < 0)) {
    n_neg <- sum(precip_vec < 0)
    stop("Precipitation cannot be negative.\n",
         "  Found ", n_neg, " negative value(s) in precip_vec.")
  }

  precip_fun <- approxfun(times, precip_vec, method = "constant", rule = 2)
  pars  <- c(k1 = k1, k2 = k2, k3 = k3)
  state <- c(S1 = S1_0, S2 = S2_0)

  out <- deSolve::ode(
    y      = state,
    times  = times,
    func   = two_tank_ode,
    parms  = pars,
    precip_fun = precip_fun,
    method = "lsoda"
  )

  df <- as.data.frame(out)
  names(df) <- c("time", "S1", "S2", "Q1", "Q2", "Q_total", "P")
  return(df)
}


#' Print Model Summary
#'
#' Prints a human-readable summary of a two-tank simulation result.
#'
#' @param sim Data.frame. Output from \code{\link{run_two_tank}}.
#' @param label Character. Name to display. Default "Simulation".
#'
#' @export
print_model_summary <- function(sim, label = "Simulation") {

  n <- nrow(sim)
  cat("\n")
  cat("  ╔══════════════════════════════════════════════╗\n")
  cat(sprintf("  ║  %s  — Model Summary\n", label))
  cat("  ╠══════════════════════════════════════════════╣\n")
  cat(sprintf("  ║  Duration           : %d days\n", n))
  cat(sprintf("  ║  Total rainfall     : %.1f mm\n", sum(sim$P)))
  cat(sprintf("  ║  Total discharge    : %.1f mm\n", sum(sim$Q_total)))
  cat(sprintf("  ║  Peak discharge     : %.2f mm/day (day %d)\n",
              max(sim$Q_total), which.max(sim$Q_total) - 1))
  cat(sprintf("  ║  Peak surface (Q1)  : %.2f mm/day\n", max(sim$Q1)))
  cat(sprintf("  ║  Peak baseflow (Q2) : %.2f mm/day\n", max(sim$Q2)))
  cat(sprintf("  ║  Baseflow index     : %.1f %%\n",
              sum(sim$Q2) / sum(sim$Q_total) * 100))
  cat(sprintf("  ║  Runoff coefficient : %.3f\n",
              sum(sim$Q_total) / sum(sim$P)))
  cat("  ╚══════════════════════════════════════════════╝\n\n")
}
