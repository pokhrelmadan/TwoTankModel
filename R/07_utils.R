#' @title Utility Functions

#' Load Data from CSV Files (Config-Aware)
#'
#' Smart loader that reads rainfall and optional discharge from CSV files
#' and auto-detects column names.
#'
#' @param rainfall_file Path to rainfall CSV (columns: date, P).
#' @param discharge_file Path to discharge CSV. Optional.
#' @param discharge_units Character. "m3s" or "mm_day". Default "m3s".
#' @param area_km2 Numeric. Catchment area (required if discharge_units="m3s"
#'   and you want to compare in mm/day instead — but usually leave Q as m³/s).
#' @param date_format Character. Default "%Y-%m-%d".
#'
#' @return List: dates, times, precip, Q_obs, n_days, units, area_km2.
#' @export
#' @param area_km2 Numeric. Catchment area.
#' @param date_format Character. Date format. Default "auto" — tries common
#'   formats automatically: "%Y-%m-%d", "%m/%d/%Y", "%d/%m/%Y", "%Y/%m/%d".
#'
#' @return List: dates, times, precip, Q_obs, n_days, units, area_km2.
#' @export
load_data <- function(rainfall_file, discharge_file = NULL,
                      discharge_units = "m3s",
                      area_km2 = NULL,
                      date_format = "auto") {

  # Helper: parse dates trying common formats if "auto"
  parse_dates <- function(x, fmt) {
    if (fmt == "auto") {
      formats_to_try <- c("%Y-%m-%d", "%m/%d/%Y", "%d/%m/%Y", "%Y/%m/%d",
                          "%d-%m-%Y", "%m-%d-%Y")
      for (f in formats_to_try) {
        d <- as.Date(x, format = f)
        if (sum(!is.na(d)) > 0.9 * length(x)) {
          cat(sprintf("    Date format detected: %s\n", f))
          return(d)
        }
      }
      stop("Could not auto-detect date format. Set DATE_FORMAT in config.R, e.g. '%m/%d/%Y'")
    } else {
      as.Date(x, format = fmt)
    }
  }
  cat("\n  Loading data...\n")

  if (!file.exists(rainfall_file))
    stop("Rainfall file not found: ", rainfall_file)

  rain_df <- read.csv(rainfall_file, stringsAsFactors = FALSE)

  date_col <- intersect(tolower(names(rain_df)), c("date", "time", "day"))[1]
  if (is.na(date_col)) stop("No date column in ", rainfall_file)
  date_col <- names(rain_df)[tolower(names(rain_df)) == date_col]

  p_col <- intersect(tolower(names(rain_df)),
                     c("p", "p_mm", "precipitation", "rainfall", "precip"))[1]
  if (is.na(p_col)) stop("No precipitation column in ", rainfall_file)
  p_col <- names(rain_df)[tolower(names(rain_df)) == p_col]

  dates  <- parse_dates(rain_df[[date_col]], date_format)
  precip <- as.numeric(rain_df[[p_col]])
  n_days <- length(dates)
  times  <- 0:(n_days - 1)

  if (any(is.na(dates))) stop("Could not parse all dates. ",
			      sum(is.na(dates)), " NA dates found.")
  precip[is.na(precip)] <- 0
  precip[precip < 0]    <- 0

  cat(sprintf("    Rainfall: %s to %s (%d days)\n",
              min(dates), max(dates), n_days))
  cat(sprintf("    Total = %.0f mm | Max = %.1f mm/day\n",
              sum(precip), max(precip)))

  result <- list(dates = dates, times = times, precip = precip,
                 n_days = n_days, Q_obs = NULL,
                 units = discharge_units, area_km2 = area_km2)

  if (!is.null(discharge_file)) {
    if (!file.exists(discharge_file))
      stop("Discharge file not found: ", discharge_file)

    q_df <- read.csv(discharge_file, stringsAsFactors = FALSE)
    q_col <- intersect(tolower(names(q_df)),
                       c("q", "q_m3s", "q_mm_day", "q_mm",
                         "discharge", "flow", "streamflow"))[1]
    if (is.na(q_col)) stop("No discharge column in ", discharge_file)
    q_col <- names(q_df)[tolower(names(q_df)) == q_col]

    Q_obs <- as.numeric(q_df[[q_col]])
    Q_obs[is.na(Q_obs)] <- 0
    Q_obs[Q_obs < 0]    <- 0

    if (length(Q_obs) != n_days)
      stop("Discharge has ", length(Q_obs), " values but rainfall has ",
           n_days, ".")

    unit_lbl <- if (discharge_units == "m3s") "m³/s" else "mm/day"
    cat(sprintf("    Discharge (%s): mean=%.3f, max=%.3f\n",
                unit_lbl, mean(Q_obs), max(Q_obs)))

    if (discharge_units == "m3s" && is.null(area_km2))
      warning("Discharge is in m³/s but area_km2 not provided. ",
              "Simulated Q cannot be converted for comparison!")

    result$Q_obs <- Q_obs
  }

  cat("    ✓ Data loaded.\n\n")
  return(result)
}


#' Generate Synthetic Daily Rainfall (for testing)
#' @export
generate_daily_rainfall <- function(n_days, dates, seed = 42,
                                    n_storms = 5, storm_range = c(30, 80)) {
  if (length(dates) != n_days) stop("dates must have ", n_days, " elements.")
  set.seed(seed)
  P   <- numeric(n_days)
  doy <- as.numeric(format(dates, "%j"))
  rain_prob <- pmin(pmax(0.15 + 0.35 * sin(pi * (doy - 60) / 180), 0.05), 0.60)
  for (i in 1:n_days) {
    if (runif(1) < rain_prob[i])
      P[i] <- rexp(1, rate = 1 / max(8 + 20 * sin(pi * (doy[i] - 60) / 180), 3))
  }
  wet <- which(doy > 120 & doy < 280)
  if (length(wet) >= n_storms && n_storms > 0) {
    P[sample(wet, n_storms)] <- P[sample(wet, n_storms)] +
                                runif(n_storms, storm_range[1], storm_range[2])
  }
  return(round(P, 1))
}


#' Convert m³/s to mm/day
#' @export
m3s_to_mmday <- function(Q_m3s, area_km2) {
  Q_m3s * 86400 / (area_km2 * 1e6) * 1000
}

#' Convert mm/day to m³/s
#' @export
mmday_to_m3s <- function(Q_mm, area_km2) {
  Q_mm * area_km2 * 1e6 / (86400 * 1000)
}


#' Create Run Folder
#'
#' Creates a timestamped results folder for a calibration run.
#'
#' @param results_dir Character. Parent directory. Default "results".
#' @param run_name Character. Subfolder name. Default "run_01".
#' @param overwrite Logical. If FALSE (default), appends timestamp.
#'
#' @return Character. Full path to created folder.
#' @export
create_run_folder <- function(results_dir = "results",
                              run_name = "run_01",
                              overwrite = FALSE) {
  run_name <- gsub(" ", "_", run_name)
  path <- file.path(results_dir, run_name)

  if (dir.exists(path) && !overwrite) {
    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    path <- file.path(results_dir, paste0(run_name, "_", ts))
  }

  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(path, "plots"), showWarnings = FALSE)
  dir.create(file.path(path, "csv"),   showWarnings = FALSE)

  cat(sprintf("  Results folder: %s\n\n", path))
  return(path)
}


#' Save Parameters to Text File
#' @export
save_parameters <- function(cal_result, path) {
  p <- cal_result$best_params
  m <- cal_result$best_metrics

  lines <- c(
    "===================================================",
    "TWO-TANK MODEL — CALIBRATED PARAMETERS",
    "===================================================",
    "",
    sprintf("Timestamp    : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("Samples      : %d", cal_result$n_samples),
    sprintf("Runtime      : %.2f seconds", cal_result$elapsed),
    sprintf("Objective    : %s", cal_result$objective),
    "",
    "Best Parameters [1/day]:",
    sprintf("  k1 (surface runoff) = %.6f", p[1]),
    sprintf("  k2 (percolation)    = %.6f", p[2]),
    sprintf("  k3 (baseflow)       = %.6f", p[3]),
    "",
    "Tank Residence Times:",
    sprintf("  Upper tank ≈ %.2f days", 1 / (p[1] + p[2])),
    sprintf("  Lower tank ≈ %.2f days", 1 / p[3]),
    "",
    "Performance Metrics:",
    sprintf("  NSE    = %.4f", m$NSE),
    sprintf("  KGE    = %.4f", m$KGE),
    sprintf("  LogNSE = %.4f", m$LogNSE),
    sprintf("  RMSE   = %.4f", m$RMSE),
    sprintf("  PBIAS  = %.2f %%", m$PBIAS),
    ""
  )

  writeLines(lines, path)
  cat(sprintf("    ✓ %s\n", basename(path)))
}


#' Export All Results to CSV
#' @export
export_results <- function(cal_result, uncertainty = NULL, sensitivity = NULL,
                           dates = NULL, Q_obs = NULL, precip = NULL,
                           output_dir = ".", prefix = "twotank") {

  cat("  Exporting CSVs...\n")

  # 1. All MC samples
  f1 <- file.path(output_dir, paste0(prefix, "_mc_samples.csv"))
  write.csv(cal_result$samples, f1, row.names = FALSE)
  cat(sprintf("    ✓ %s\n", basename(f1)))

  # 2. Behavioural sets
  if (!is.null(uncertainty)) {
    f2 <- file.path(output_dir, paste0(prefix, "_behavioural.csv"))
    write.csv(uncertainty$behavioural, f2, row.names = FALSE)
    cat(sprintf("    ✓ %s (%d sets)\n", basename(f2),
                nrow(uncertainty$behavioural)))
  }

  # 3. Sensitivity
  if (!is.null(sensitivity)) {
    f_s <- file.path(output_dir, paste0(prefix, "_sensitivity.csv"))
    write.csv(sensitivity, f_s, row.names = FALSE)
    cat(sprintf("    ✓ %s\n", basename(f_s)))
  }

  # 4. Best simulation time series
  sim <- cal_result$best_sim
  out_df <- data.frame(
    time = sim$time,
    P_mm = sim$P,
    S1_mm = sim$S1, S2_mm = sim$S2,
    Q1_mm = sim$Q1, Q2_mm = sim$Q2, Q_total_mm = sim$Q_total
  )
  if (!is.null(dates)) out_df$date <- dates
  if (!is.null(cal_result$area_km2)) {
    out_df$Q1_m3s      <- sim$Q1_m3s
    out_df$Q2_m3s      <- sim$Q2_m3s
    out_df$Q_total_m3s <- sim$Q_total_m3s
  }
  if (!is.null(Q_obs)) out_df$Q_obs <- Q_obs
  if (!is.null(uncertainty)) {
    out_df$Q_lower  <- uncertainty$Q_lower
    out_df$Q_upper  <- uncertainty$Q_upper
    out_df$Q_median <- uncertainty$Q_median
  }

  # Reorder columns nicely
  if ("date" %in% names(out_df))
    out_df <- out_df[, c("date", setdiff(names(out_df), "date"))]

  f3 <- file.path(output_dir, paste0(prefix, "_best_simulation.csv"))
  write.csv(out_df, f3, row.names = FALSE)
  cat(sprintf("    ✓ %s\n", basename(f3)))

  cat("\n")
}
