#' @title Visualisation Functions
#' @description Publication-ready plots for the two-tank model including
#'   hydrographs, uncertainty envelopes, diagnostics, and analysis.
 
 
# ── Shared theme ─────────────────────────────────────────────────────────────
 
theme_hydro <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 9, colour = "grey40"),
      legend.position  = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1)
    )
}
 
 
# ── Smart date-axis scaling based on data length ───────────────────────────
#   Picks reasonable tick intervals and label formats to avoid cramming.
smart_date_scale <- function(dates) {
  n_days  <- length(dates)
  n_years <- as.numeric(diff(range(dates))) / 365.25
 
  if (n_years <= 1) {
    breaks <- "1 month";   labels <- "%b"
  } else if (n_years <= 3) {
    breaks <- "3 months";  labels <- "%b %Y"
  } else if (n_years <= 10) {
    breaks <- "1 year";    labels <- "%Y"
  } else if (n_years <= 25) {
    breaks <- "2 years";   labels <- "%Y"
  } else {
    breaks <- "5 years";   labels <- "%Y"
  }
 
  ggplot2::scale_x_date(date_breaks = breaks, date_labels = labels)
}
 
 
#' Plot Daily Rainfall Hyetograph
#'
#' @param dates Date vector.
#' @param precip Numeric vector. Daily precipitation [mm/day].
#' @return A ggplot object.
#' @export
plot_rainfall <- function(dates, precip) {
  df <- data.frame(date = dates, P = precip)
  ggplot2::ggplot(df, ggplot2::aes(date, P)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue", width = 1) +
    ggplot2::scale_y_reverse(expand = ggplot2::expansion(mult = c(0.05, 0))) +
    smart_date_scale(dates) +
    ggplot2::labs(title = "Daily Precipitation", x = NULL, y = "P (mm/day)") +
    theme_hydro()
}
 
 
#' Plot Observed vs Simulated Hydrograph
#'
#' @param dates Date vector.
#' @param Q_obs Numeric vector. Observed discharge.
#' @param sim Data.frame. Output from \code{\link{run_two_tank}}.
#' @param nse Numeric. NSE value to display. Default NULL (calculated).
#' @return A ggplot object.
#' @export
plot_hydrograph <- function(dates, Q_obs, sim, nse = NULL) {
  # Use m³/s columns if available (matches Q_obs units)
  use_m3s <- "Q_total_m3s" %in% names(sim)
  Q_sim   <- if (use_m3s) sim$Q_total_m3s else sim$Q_total
  y_label <- if (use_m3s) "Q (m³/s)" else "Q (mm/day)"
 
  if (is.null(nse)) nse <- calc_nse(Q_obs, Q_sim)
  df <- data.frame(
    date = rep(dates, 2),
    Q    = c(Q_obs, Q_sim),
    Source = rep(c("Observed", "Simulated"), each = length(dates))
  )
  ggplot2::ggplot(df, ggplot2::aes(date, Q, colour = Source)) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::scale_colour_manual(values = c("Observed" = "black",
                                             "Simulated" = "red")) +
    smart_date_scale(dates) +
    ggplot2::labs(title = sprintf("Observed vs Simulated (NSE = %.4f)", nse),
                  x = NULL, y = y_label, colour = NULL) +
    theme_hydro()
}
 
 
#' Plot Uncertainty Envelope
#'
#' @param dates Date vector.
#' @param Q_obs Numeric vector. Observed discharge.
#' @param sim Data.frame. Best simulation from \code{\link{run_two_tank}}.
#' @param uncertainty List. Output from \code{\link{extract_uncertainty}}.
#' @return A ggplot object.
#' @export
plot_uncertainty_envelope <- function(dates, Q_obs, sim, uncertainty) {
  use_m3s <- "Q_total_m3s" %in% names(sim)
  Q_sim   <- if (use_m3s) sim$Q_total_m3s else sim$Q_total
  y_label <- if (use_m3s) "Q (m³/s)" else "Q (mm/day)"
 
  df <- data.frame(
    date   = dates,
    Obs    = Q_obs,
    Best   = Q_sim,
    Median = uncertainty$Q_median,
    Lower  = uncertainty$Q_lower,
    Upper  = uncertainty$Q_upper
  )
  ggplot2::ggplot(df, ggplot2::aes(date)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Lower, ymax = Upper),
                         fill = "lightblue", alpha = 0.5) +
    ggplot2::geom_line(ggplot2::aes(y = Obs, colour = "Observed"), linewidth = 0.6) +
    ggplot2::geom_line(ggplot2::aes(y = Best, colour = "Best (MC)"), linewidth = 0.6) +
    ggplot2::geom_line(ggplot2::aes(y = Median, colour = "Median (MC)"),
                       linewidth = 0.4, linetype = "dashed") +
    ggplot2::scale_colour_manual(values = c("Observed" = "black",
                                             "Best (MC)" = "red",
                                             "Median (MC)" = "blue")) +
    smart_date_scale(dates) +
    ggplot2::labs(
      title = "Hydrograph with 90% Prediction Uncertainty Envelope",
      subtitle = sprintf("%.1f%% of observations within bounds | %d behavioural sets",
                         uncertainty$containment_pct, uncertainty$n_behavioural),
      x = NULL, y = y_label, colour = NULL
    ) +
    theme_hydro()
}
 
 
#' Plot Observed vs Simulated Scatter
#'
#' @param Q_obs Numeric vector. Observed discharge.
#' @param sim Data.frame. Output from \code{\link{run_two_tank}}.
#' @return A ggplot object.
#' @export
plot_scatter <- function(Q_obs, sim) {
  use_m3s <- "Q_total_m3s" %in% names(sim)
  Q_sim   <- if (use_m3s) sim$Q_total_m3s else sim$Q_total
  unit    <- if (use_m3s) "m³/s" else "mm/day"
 
  nse <- calc_nse(Q_obs, Q_sim)
  df  <- data.frame(Obs = Q_obs, Sim = Q_sim)
  ggplot2::ggplot(df, ggplot2::aes(Obs, Sim)) +
    ggplot2::geom_point(alpha = 0.4, size = 1.5, colour = "steelblue") +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", colour = "red") +
    ggplot2::labs(title = sprintf("Observed vs Simulated (NSE = %.4f)", nse),
                  x = sprintf("Observed Q (%s)", unit),
                  y = sprintf("Simulated Q (%s)", unit)) +
    ggplot2::coord_equal() +
    theme_hydro()
}
 
 
#' Plot Flow Components (Stacked Area)
#'
#' @param dates Date vector.
#' @param sim Data.frame. Output from \code{\link{run_two_tank}}.
#' @return A ggplot object.
#' @export
plot_components <- function(dates, sim) {
  use_m3s <- "Q_total_m3s" %in% names(sim)
  Q1 <- if (use_m3s) sim$Q1_m3s else sim$Q1
  Q2 <- if (use_m3s) sim$Q2_m3s else sim$Q2
  unit <- if (use_m3s) "m³/s" else "mm/day"
 
  bfi <- sum(Q2) / (sum(Q1) + sum(Q2)) * 100
  df <- data.frame(
    date = rep(dates, 2),
    Q    = c(Q1, Q2),
    Flow = rep(c("Q1 Surface Runoff", "Q2 Baseflow"), each = length(dates))
  )
  ggplot2::ggplot(df, ggplot2::aes(date, Q, fill = Flow)) +
    ggplot2::geom_area(alpha = 0.7, position = "stack") +
    ggplot2::scale_fill_manual(values = c("Q1 Surface Runoff" = "tomato",
                                           "Q2 Baseflow" = "dodgerblue3")) +
    smart_date_scale(dates) +
    ggplot2::labs(title = sprintf("Flow Components (BFI = %.1f%%)", bfi),
                  x = NULL, y = sprintf("Q (%s)", unit)) +
    theme_hydro()
}
 
 
#' Plot Tank Storage Dynamics
#'
#' @param dates Date vector.
#' @param sim Data.frame. Output from \code{\link{run_two_tank}}.
#' @return A ggplot object.
#' @export
plot_storage <- function(dates, sim) {
  df <- data.frame(
    date = rep(dates, 2),
    S    = c(sim$S1, sim$S2),
    Tank = rep(c("S1 Upper Tank", "S2 Lower Tank"), each = length(dates))
  )
  ggplot2::ggplot(df, ggplot2::aes(date, S, colour = Tank)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::scale_colour_manual(values = c("S1 Upper Tank" = "orange",
                                             "S2 Lower Tank" = "purple")) +
    smart_date_scale(dates) +
    ggplot2::labs(title = "Tank Storage Dynamics",
                  x = NULL, y = "Storage (mm)") +
    theme_hydro()
}
 
 
#' Plot Dotty Plots for Parameter Identifiability
#'
#' @param cal_result List. Output from \code{\link{calibrate_montecarlo}}.
#' @param nse_threshold Numeric. Behavioural threshold line. Default 0.5.
#' @return A ggplot object.
#' @export
plot_dotty <- function(cal_result, nse_threshold = 0.5) {
  s <- cal_result$samples
  has_k4 <- "k4" %in% names(s) && any(s$k4 > 0)
 
  mc_long <- rbind(
    data.frame(Parameter = "k1 (surface runoff)", Value = s$k1, NSE = s$NSE),
    data.frame(Parameter = "k2 (percolation)",    Value = s$k2, NSE = s$NSE),
    data.frame(Parameter = "k3 (baseflow)",       Value = s$k3, NSE = s$NSE)
  )
  if (has_k4) {
    mc_long <- rbind(mc_long,
      data.frame(Parameter = "k4 (ET)", Value = s$k4, NSE = s$NSE))
  }
 
  ggplot2::ggplot(mc_long, ggplot2::aes(Value, NSE)) +
    ggplot2::geom_point(alpha = 0.15, size = 0.6, colour = "grey40") +
    ggplot2::geom_hline(yintercept = nse_threshold,
                        linetype = "dashed", colour = "red") +
    ggplot2::facet_wrap(~Parameter, scales = "free_x", ncol = 2) +
    ggplot2::labs(
      title = "Parameter Identifiability (Dotty Plots)",
      subtitle = "Sharp peak = well identified | Flat = equifinal (uncertain)",
      x = "Parameter Value", y = "NSE"
    ) +
    theme_hydro()
}
 
 
#' Plot NSE Distribution Histogram
#'
#' @param cal_result List. Output from \code{\link{calibrate_montecarlo}}.
#' @param nse_threshold Numeric. Behavioural threshold. Default 0.5.
#' @return A ggplot object.
#' @export
plot_nse_histogram <- function(cal_result, nse_threshold = 0.5) {
  s <- cal_result$samples
  best_nse <- cal_result$best_metrics$NSE
  ggplot2::ggplot(s, ggplot2::aes(NSE)) +
    ggplot2::geom_histogram(bins = 80, fill = "steelblue",
                            colour = "white", alpha = 0.8) +
    ggplot2::geom_vline(xintercept = nse_threshold,
                        linetype = "dashed", colour = "red") +
    ggplot2::geom_vline(xintercept = best_nse, colour = "gold", linewidth = 1.2) +
    ggplot2::labs(
      title = "NSE Distribution Across Monte Carlo Samples",
      subtitle = sprintf("Gold = best (%.4f) | Red = threshold (%.2f)",
                         best_nse, nse_threshold),
      x = "NSE", y = "Count"
    ) +
    theme_hydro()
}
 
 
#' Plot Parameter Correlation (Behavioural Sets)
#'
#' @param uncertainty List. Output from \code{\link{extract_uncertainty}}.
#' @param x_param Character. Parameter for x-axis: "k1", "k2", or "k3".
#' @param y_param Character. Parameter for y-axis: "k1", "k2", or "k3".
#' @return A ggplot object.
#' @export
plot_param_correlation <- function(uncertainty, x_param = "k1", y_param = "k2") {
  behav <- uncertainty$behavioural
  labels <- c(k1 = "k1 (surface runoff)",
              k2 = "k2 (percolation)",
              k3 = "k3 (baseflow)")
  ggplot2::ggplot(behav, ggplot2::aes(.data[[x_param]], .data[[y_param]],
                                       colour = NSE)) +
    ggplot2::geom_point(alpha = 0.5, size = 1.2) +
    ggplot2::scale_colour_viridis_c(option = "plasma") +
    ggplot2::labs(
      title = sprintf("%s vs %s (behavioural sets)", labels[x_param], labels[y_param]),
      x = labels[x_param], y = labels[y_param]
    ) +
    theme_hydro()
}
 
 
#' Plot Flow Duration Curve
#'
#' @param Q_obs Numeric vector. Observed discharge.
#' @param sim Data.frame. Output from \code{\link{run_two_tank}}.
#' @return A ggplot object.
#' @export
plot_flow_duration <- function(Q_obs, sim) {
  use_m3s <- "Q_total_m3s" %in% names(sim)
  Q_sim   <- if (use_m3s) sim$Q_total_m3s else sim$Q_total
  unit    <- if (use_m3s) "m³/s" else "mm/day"
 
  n <- length(Q_obs)
  fdc_obs <- sort(Q_obs, decreasing = TRUE)
  fdc_sim <- sort(Q_sim, decreasing = TRUE)
  exceed  <- (1:n) / n * 100
 
  df <- data.frame(
    Exceedance = rep(exceed, 2),
    Q      = c(fdc_obs, fdc_sim),
    Source = rep(c("Observed", "Simulated"), each = n)
  )
  ggplot2::ggplot(df, ggplot2::aes(Exceedance, Q, colour = Source)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_colour_manual(values = c("Observed" = "black",
                                             "Simulated" = "red")) +
    ggplot2::labs(title = "Flow Duration Curve",
                  x = "Exceedance Probability (%)",
                  y = sprintf("Q (%s, log scale)", unit)) +
    theme_hydro()
}
 
 
#' Generate All Plots and Save to PDF
#'
#' Creates a multi-page PDF with all diagnostic plots and also
#' saves individual PNGs.
#'
#' @param dates Date vector.
#' @param precip Numeric vector. Daily precipitation.
#' @param Q_obs Numeric vector. Observed discharge.
#' @param cal_result List. Output from \code{\link{calibrate_montecarlo}}.
#' @param uncertainty List. Output from \code{\link{extract_uncertainty}}.
#' @param output_dir Character. Directory for output files. Default ".".
#' @param prefix Character. Filename prefix. Default "twotank".
#'
#' @return Invisible list of all ggplot objects.
#' @export
plot_all <- function(dates, precip, Q_obs, cal_result, uncertainty,
                     output_dir = ".", prefix = "twotank") {
 
  sim <- cal_result$best_sim
  has_k4 <- "k4" %in% names(uncertainty$behavioural) &&
            any(uncertainty$behavioural$k4 > 0)
 
  p1  <- plot_rainfall(dates, precip)
  p2  <- plot_uncertainty_envelope(dates, Q_obs, sim, uncertainty)
  p3  <- plot_hydrograph(dates, Q_obs, sim)
  p4  <- plot_scatter(Q_obs, sim)
  p5  <- plot_components(dates, sim)
  p6  <- plot_storage(dates, sim)
  p7  <- plot_dotty(cal_result)
  p8  <- plot_nse_histogram(cal_result)
  p9  <- plot_param_correlation(uncertainty, "k1", "k2")
  p10 <- plot_param_correlation(uncertainty, "k1", "k3")
  p11 <- plot_param_correlation(uncertainty, "k2", "k3")
  p12 <- plot_flow_duration(Q_obs, sim)
 
  # Extra k4 correlation plots if ET is enabled
  if (has_k4) {
    p13 <- plot_param_correlation(uncertainty, "k1", "k4")
    p14 <- plot_param_correlation(uncertainty, "k2", "k4")
    p15 <- plot_param_correlation(uncertainty, "k3", "k4")
  }
 
  # Save combined PDF
  pdf_path <- file.path(output_dir, paste0(prefix, "_results.pdf"))
  pdf_height <- if (has_k4) 32 else 28
  grDevices::pdf(pdf_path, width = 12, height = pdf_height)
  if (has_k4) {
    gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p12,
                             ncol = 1,
                             heights = c(0.7, 1.2, 1, 1, 1, 1, 1.3, 0.8, 1))
  } else {
    gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p12,
                             ncol = 1,
                             heights = c(0.7, 1.2, 1, 1, 1, 1, 1, 0.8, 1))
  }
  grDevices::dev.off()
 
  # Save individual PNGs
  ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_rainfall.png")),
                  p1, width = 12, height = 3, dpi = 150)
  ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_envelope.png")),
                  p2, width = 12, height = 5, dpi = 150)
  ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_hydrograph.png")),
                  p3, width = 12, height = 5, dpi = 150)
  ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_scatter.png")),
                  p4, width = 7, height = 7, dpi = 150)
  ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_components.png")),
                  p5, width = 12, height = 4, dpi = 150)
  ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_storage.png")),
                  p6, width = 12, height = 4, dpi = 150)
  dotty_h <- if (has_k4) 7 else 5
  ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_dotty.png")),
                  p7, width = 12, height = dotty_h, dpi = 150)
  ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_nse_hist.png")),
                  p8, width = 10, height = 4, dpi = 150)
  ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_corr_k1k2.png")),
                  p9, width = 7, height = 6, dpi = 150)
  ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_corr_k1k3.png")),
                  p10, width = 7, height = 6, dpi = 150)
  ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_corr_k2k3.png")),
                  p11, width = 7, height = 6, dpi = 150)
  ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_fdc.png")),
                  p12, width = 10, height = 5, dpi = 150)
  if (has_k4) {
    ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_corr_k1k4.png")),
                    p13, width = 7, height = 6, dpi = 150)
    ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_corr_k2k4.png")),
                    p14, width = 7, height = 6, dpi = 150)
    ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_corr_k3k4.png")),
                    p15, width = 7, height = 6, dpi = 150)
  }
 
  n_plots <- if (has_k4) 15 else 12
  cat(sprintf("  Saved: %s + %d PNGs in %s\n", basename(pdf_path), n_plots, output_dir))
 
  out_list <- list(rainfall = p1, envelope = p2, hydrograph = p3,
                   scatter = p4, components = p5, storage = p6,
                   dotty = p7, nse_hist = p8,
                   corr_k1k2 = p9, corr_k1k3 = p10, corr_k2k3 = p11,
                   fdc = p12)
  if (has_k4) {
    out_list$corr_k1k4 <- p13
    out_list$corr_k2k4 <- p14
    out_list$corr_k3k4 <- p15
  }
  invisible(out_list)
}
