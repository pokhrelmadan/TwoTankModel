#' @title Performance Metrics for Hydrological Model Evaluation
#' @description Functions to compute goodness-of-fit statistics between
#'   observed and simulated discharge, with human-readable interpretation.


#' Nash-Sutcliffe Efficiency (NSE)
#'
#' Measures overall model fit. Emphasises high flows (peaks).
#'
#' \strong{Interpretation:}
#' \itemize{
#'   \item 1.0       = Perfect match
#'   \item 0.7 – 1.0 = Good
#'   \item 0.5 – 0.7 = Acceptable
#'   \item < 0.5     = Poor
#'   \item < 0       = Worse than using the observed mean
#' }
#'
#' @param obs Numeric vector. Observed discharge.
#' @param sim Numeric vector. Simulated discharge (same length as obs).
#' @return Numeric scalar. NSE value.
#' @export
calc_nse <- function(obs, sim) {
  if (length(obs) != length(sim))
    stop("obs and sim must have the same length (",
         length(obs), " vs ", length(sim), ").")
  1 - sum((obs - sim)^2) / sum((obs - mean(obs))^2)
}


#' Kling-Gupta Efficiency (KGE)
#'
#' Balanced metric: correlation + variability bias + mean bias.
#' Often preferred over NSE for calibration (Gupta et al., 2009).
#'
#' \strong{Interpretation:}
#' \itemize{
#'   \item 1.0       = Perfect
#'   \item > 0.5     = Acceptable
#'   \item < -0.41   = Worse than using the mean
#' }
#'
#' @param obs Numeric vector. Observed discharge.
#' @param sim Numeric vector. Simulated discharge.
#' @return Numeric scalar. KGE value.
#' @export
calc_kge <- function(obs, sim) {
  if (length(obs) != length(sim))
    stop("obs and sim must have the same length.")
  r     <- cor(sim, obs)
  alpha <- sd(sim) / sd(obs)
  beta  <- mean(sim) / mean(obs)
  1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2)
}


#' Log Nash-Sutcliffe Efficiency (LogNSE)
#'
#' NSE on log-transformed flows. Emphasises low flows and baseflow.
#' Use alongside NSE for a complete picture.
#'
#' @param obs Numeric vector. Observed discharge.
#' @param sim Numeric vector. Simulated discharge.
#' @return Numeric scalar. LogNSE value.
#' @export
calc_lognse <- function(obs, sim) {
  if (length(obs) != length(sim))
    stop("obs and sim must have the same length.")
  eps <- 0.01 * mean(obs[obs > 0])
  lo  <- log(obs + eps)
  ls  <- log(sim + eps)
  1 - sum((lo - ls)^2) / sum((lo - mean(lo))^2)
}


#' Root Mean Square Error (RMSE)
#'
#' Average error magnitude in the same units as discharge \[mm/day\].
#' Lower = better. Compare against the standard deviation of observations.
#'
#' @param obs Numeric vector. Observed discharge.
#' @param sim Numeric vector. Simulated discharge.
#' @return Numeric scalar. RMSE value.
#' @export
calc_rmse <- function(obs, sim) {
  if (length(obs) != length(sim))
    stop("obs and sim must have the same length.")
  sqrt(mean((obs - sim)^2))
}


#' Percent Bias (PBIAS)
#'
#' Systematic tendency to over- or under-predict.
#'
#' \strong{Interpretation:}
#' \itemize{
#'   \item Positive = model over-predicts (too much water)
#'   \item Negative = model under-predicts (too little water)
#'   \item Ideal = 0\%
#'   \item |PBIAS| < 10\% = Good
#'   \item |PBIAS| < 25\% = Acceptable
#' }
#'
#' @param obs Numeric vector. Observed discharge.
#' @param sim Numeric vector. Simulated discharge.
#' @return Numeric scalar. PBIAS in percent.
#' @export
calc_pbias <- function(obs, sim) {
  if (length(obs) != length(sim))
    stop("obs and sim must have the same length.")
  sum(sim - obs) / sum(obs) * 100
}


#' Compute All Performance Metrics At Once
#'
#' Convenience function that returns all five metrics in a tidy list.
#'
#' @param obs Numeric vector. Observed discharge.
#' @param sim Numeric vector. Simulated discharge.
#' @param print Logical. Print a formatted summary. Default TRUE.
#' @return A named list with NSE, KGE, LogNSE, RMSE, PBIAS.
#'
#' @examples
#' \dontrun{
#' m <- calc_all_metrics(Q_obs, sim$Q_total)
#' }
#'
#' @export
calc_all_metrics <- function(obs, sim, print = TRUE) {
  m <- list(
    NSE    = calc_nse(obs, sim),
    KGE    = calc_kge(obs, sim),
    LogNSE = calc_lognse(obs, sim),
    RMSE   = calc_rmse(obs, sim),
    PBIAS  = calc_pbias(obs, sim)
  )

  if (print) {
    # Rating helper
    rate_nse <- function(x) {
      if (x >= 0.7) return("Good")
      if (x >= 0.5) return("Acceptable")
      if (x >= 0)   return("Poor")
      return("Unacceptable")
    }
    rate_pbias <- function(x) {
      if (abs(x) < 10)  return("Good")
      if (abs(x) < 25)  return("Acceptable")
      return("Poor")
    }

    cat("\n  ┌─────────────────────────────────────────────┐\n")
    cat("  │         MODEL PERFORMANCE REPORT            │\n")
    cat("  ├────────────┬──────────┬─────────────────────┤\n")
    cat("  │ Metric     │ Value    │ Rating              │\n")
    cat("  ├────────────┼──────────┼─────────────────────┤\n")
    cat(sprintf("  │ NSE        │ %7.4f  │ %-20s│\n", m$NSE, rate_nse(m$NSE)))
    cat(sprintf("  │ KGE        │ %7.4f  │ %-20s│\n", m$KGE, rate_nse(m$KGE)))
    cat(sprintf("  │ LogNSE     │ %7.4f  │ %-20s│\n", m$LogNSE, rate_nse(m$LogNSE)))
    cat(sprintf("  │ RMSE       │ %7.4f  │ mm/day              │\n", m$RMSE))
    cat(sprintf("  │ PBIAS      │ %6.2f%%  │ %-20s│\n", m$PBIAS, rate_pbias(m$PBIAS)))
    cat("  └────────────┴──────────┴─────────────────────┘\n\n")
  }

  invisible(m)
}
