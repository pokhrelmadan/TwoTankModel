#' @title Utility Functions
#' @description Helpers for data loading, rainfall generation, and
#'   unit conversions.


#' Load Your Own Data
#'
#' Convenience function to load daily rainfall and optional observed
#' discharge from CSV files. Handles common date formats and validates
#' the data before returning.
#'
#' @param rainfall_file Character. Path to rainfall CSV. Must contain
#'   columns 'date' and 'P' (or 'P_mm' or 'precipitation').
#' @param discharge_file Character. Optional. Path to observed discharge CSV.
#'   Must contain 'date' and 'Q' (or 'Q_mm_day' or 'discharge').
#' @param area_km2 Numeric. Optional. If provided AND discharge is in m³/s,
#'   converts to mm/day automatically. Default NULL (no conversion).
#' @param date_format Character. Date format string. Default "\%Y-\%m-\%d".
#'
#' @return A list with: dates, times, precip, Q_obs (if provided), n_days.
#'
#' @examples
#' \dontrun{
#' # Rainfall only
#' data <- load_data("daily_rain.csv")
#'
#' # Rainfall + observed discharge
#' data <- load_data("rain.csv", "discharge.csv")
#'
#' # With unit conversion from m³/s
#' data <- load_data("rain.csv", "discharge_m3s.csv", area_km2 = 150)
#' }
#'
#' @export
load_data <- function(rainfall_file, discharge_file = NULL,
                      area_km2 = NULL, date_format = "%Y-%m-%d") {

  cat("\n  Loading data...\n")

  # ── Load rainfall ──
  if (!file.exists(rainfall_file))
    stop("Rainfall file not found: ", rainfall_file)

  rain_df <- read.csv(rainfall_file, stringsAsFactors = FALSE)

  # Find date column
  date_col <- intersect(tolower(names(rain_df)), c("date", "time", "day"))[1]
  if (is.na(date_col)) stop("Cannot find a 'date' column in ", rainfall_file)
  date_col <- names(rain_df)[tolower(names(rain_df)) == date_col]

  # Find precipitation column
  p_col <- intersect(tolower(names(rain_df)),
                     c("p", "p_mm", "precipitation", "rainfall", "precip"))[1]
  if (is.na(p_col)) stop("Cannot find a precipitation column (P, P_mm, precipitation, rainfall) in ",
                          rainfall_file)
  p_col <- names(rain_df)[tolower(names(rain_df)) == p_col]

  dates  <- as.Date(rain_df[[date_col]], format = date_format)
  precip <- as.numeric(rain_df[[p_col]])
  n_days <- length(dates)
  times  <- 0:(n_days - 1)

  if (any(is.na(dates))) stop("Could not parse dates. Check date_format parameter.")
  if (any(is.na(precip))) warning("Found NA values in precipitation. Replacing with 0.")
  precip[is.na(precip)] <- 0

  cat(sprintf("    Rainfall: %s to %s (%d days)\n",
              min(dates), max(dates), n_days))
  cat(sprintf("    Total = %.0f mm | Max = %.1f mm/day\n", sum(precip), max(precip)))

  result <- list(dates = dates, times = times, precip = precip,
                 n_days = n_days, Q_obs = NULL)

  # ── Load discharge (optional) ──
  if (!is.null(discharge_file)) {
    if (!file.exists(discharge_file))
      stop("Discharge file not found: ", discharge_file)

    q_df <- read.csv(discharge_file, stringsAsFactors = FALSE)

    q_col <- intersect(tolower(names(q_df)),
                       c("q", "q_mm_day", "q_mm", "discharge",
                         "flow", "q_m3s", "streamflow"))[1]
    if (is.na(q_col)) stop("Cannot find a discharge column in ", discharge_file)
    q_col <- names(q_df)[tolower(names(q_df)) == q_col]

    Q_obs <- as.numeric(q_df[[q_col]])

    # Convert m³/s to mm/day if area provided
    if (!is.null(area_km2)) {
      cat(sprintf("    Converting Q from m³/s to mm/day (area = %.1f km²)\n", area_km2))
      Q_obs <- m3s_to_mmday(Q_obs, area_km2)
    }

    if (length(Q_obs) != n_days)
      stop("Discharge data has ", length(Q_obs), " values but rainfall has ",
           n_days, ". They must match.")

    Q_obs[is.na(Q_obs)] <- 0
    Q_obs[Q_obs < 0]    <- 0

    cat(sprintf("    Discharge: mean = %.2f mm/day | max = %.2f mm/day\n",
                mean(Q_obs), max(Q_obs)))

    result$Q_obs <- Q_obs
  }

  cat("    ✓ Data loaded successfully.\n\n")
  return(result)
}


#' Generate Synthetic Daily Rainfall
#'
#' Creates a realistic daily rainfall time series with seasonal variation
#' and large storm events. Useful for testing and demonstrations.
#'
#' @param n_days Integer. Number of days.
#' @param dates Date vector. Must have length \code{n_days}.
#' @param seed Integer. Random seed. Default 42.
#' @param n_storms Integer. Number of large storms to inject. Default 5.
#' @param storm_range Numeric(2). Min/max storm depth \[mm\]. Default c(30, 80).
#'
#' @return Numeric vector of daily precipitation \[mm/day\].
#'
#' @examples
#' dates  <- seq(as.Date("2024-01-01"), as.Date("2024-12-31"), by = "day")
#' precip <- generate_daily_rainfall(length(dates), dates)
#' plot(dates, precip, type = "h", ylab = "P (mm/day)")
#'
#' @export
generate_daily_rainfall <- function(n_days, dates, seed = 42,
                                    n_storms = 5,
                                    storm_range = c(30, 80)) {
  if (length(dates) != n_days)
    stop("dates must have exactly ", n_days, " elements.")

  set.seed(seed)
  P   <- numeric(n_days)
  doy <- as.numeric(format(dates, "%j"))

  rain_prob <- 0.15 + 0.35 * sin(pi * (doy - 60) / 180)
  rain_prob <- pmin(pmax(rain_prob, 0.05), 0.60)

  for (i in 1:n_days) {
    if (runif(1) < rain_prob[i]) {
      mean_depth <- max(8 + 20 * sin(pi * (doy[i] - 60) / 180), 3)
      P[i] <- rexp(1, rate = 1 / mean_depth)
    }
  }

  wet_season <- which(doy > 120 & doy < 280)
  if (length(wet_season) >= n_storms && n_storms > 0) {
    storm_days <- sample(wet_season, n_storms)
    P[storm_days] <- P[storm_days] + runif(n_storms, storm_range[1], storm_range[2])
  }

  return(round(P, 1))
}


#' Convert m³/s to mm/day
#'
#' @param Q_m3s Numeric vector. Discharge in m³/s.
#' @param area_km2 Numeric. Catchment area in km².
#' @return Numeric vector. Discharge in mm/day.
#' @export
m3s_to_mmday <- function(Q_m3s, area_km2) {
  Q_m3s * 86400 / (area_km2 * 1e6) * 1000
}


#' Convert mm/day to m³/s
#'
#' @param Q_mm Numeric vector. Discharge in mm/day.
#' @param area_km2 Numeric. Catchment area in km².
#' @return Numeric vector. Discharge in m³/s.
#' @export
mmday_to_m3s <- function(Q_mm, area_km2) {
  Q_mm * area_km2 * 1e6 / (86400 * 1000)
}


#' Export All Results to CSV
#'
#' Saves calibration results, behavioural sets, and best simulation
#' to CSV files in the specified directory.
#'
#' @param cal_result List. Output from \code{\link{calibrate_montecarlo}}.
#' @param uncertainty List. Output from \code{\link{extract_uncertainty}}.
#'   Optional. Default NULL.
#' @param dates Date vector. Optional. Default NULL.
#' @param Q_obs Numeric vector. Optional observed Q. Default NULL.
#' @param output_dir Character. Output directory. Default ".".
#' @param prefix Character. Filename prefix. Default "twotank".
#'
#' @export
export_results <- function(cal_result, uncertainty = NULL,
                           dates = NULL, Q_obs = NULL,
                           output_dir = ".", prefix = "twotank") {

  cat("\n  Exporting results to CSV...\n")

  # All MC samples
  f1 <- file.path(output_dir, paste0(prefix, "_mc_samples.csv"))
  write.csv(cal_result$samples, f1, row.names = FALSE)
  cat(sprintf("    ✓ %s (%d samples)\n", basename(f1), nrow(cal_result$samples)))

  # Behavioural sets
  if (!is.null(uncertainty)) {
    f2 <- file.path(output_dir, paste0(prefix, "_behavioural.csv"))
    write.csv(uncertainty$behavioural, f2, row.names = FALSE)
    cat(sprintf("    ✓ %s (%d sets)\n", basename(f2), nrow(uncertainty$behavioural)))
  }

  # Best simulation time series
  sim_df <- cal_result$best_sim
  if (!is.null(dates)) sim_df$date <- dates
  if (!is.null(Q_obs)) sim_df$Q_obs <- Q_obs
  if (!is.null(uncertainty)) {
    sim_df$Q_lower  <- uncertainty$Q_lower
    sim_df$Q_upper  <- uncertainty$Q_upper
    sim_df$Q_median <- uncertainty$Q_median
  }
  f3 <- file.path(output_dir, paste0(prefix, "_best_simulation.csv"))
  write.csv(sim_df, f3, row.names = FALSE)
  cat(sprintf("    ✓ %s (%d days)\n", basename(f3), nrow(sim_df)))

  cat("  Done.\n\n")
}
