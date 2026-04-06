# TwoTankModel 
 
A user-friendly two-tank conceptual rainfall–runoff model for daily streamflow simulation with Monte Carlo calibration in R.
 
## Installation
 
```r
# Install from GitHub (one line!)
install.packages("devtools")
devtools::install_github("pokhrelmadan/TwoTankModel")
```
Lightest way to install from GitHub
```r
Rscript -e 'install.packages("remotes", repos="https://cloud.r-project.org"); remotes::install_github("pokhrelmadan/TwoTankModel")'
```
 
Then load it:
```r
library(TwoTankModel)
```
 
## Model Structure
 
```
     Precipitation P(t)  [mm/day]
           │
           ▼
  ┌────────────────────┐
  │    Upper Tank       │──► Q1(t) = k1·S1  (Surface Runoff & Interflow)
  │    Storage: S1      │
  └────────┬───────────┘         Total: Q(t) = Q1(t) + Q2(t)
           │ k2·S1
           │ (Percolation)
           ▼
  ┌────────────────────┐
  │    Lower Tank       │──► Q2(t) = k3·S2  (Baseflow)
  │    Storage: S2      │
  └────────────────────┘
```
 
**Governing equations:**
 
| Equation | Description |
|----------|-------------|
| dS1/dt = P(t) - k1·S1 - k2·S1 | Upper tank water balance |
| dS2/dt = k2·S1 - k3·S2 | Lower tank water balance |
| Q(t) = k1·S1 + k3·S2 | Total discharge |
 
**Parameters:**
 
| Parameter | Description | Unit | Typical Range |
|-----------|-------------|------|---------------|
| k1 | Surface runoff coefficient | 1/day | 0.01 – 0.80 |
| k2 | Percolation coefficient | 1/day | 0.01 – 0.50 |
| k3 | Baseflow coefficient | 1/day | 0.001 – 0.15 |
 
## Quick Start
 
### Option 1: Demo (no data needed)
 
```r
library(TwoTankModel)
 
# Generate synthetic data
dates  <- seq(as.Date("2024-01-01"), as.Date("2024-12-31"), by = "day")
times  <- 0:(length(dates) - 1)
precip <- generate_daily_rainfall(length(dates), dates)
 
# Run a simulation
sim <- run_two_tank(k1 = 0.30, k2 = 0.10, k3 = 0.025, times, precip)
print_model_summary(sim)
plot(dates, sim$Q_total, type = "l", xlab = "Date", ylab = "Q (mm/day)")
```
 
### Option 2: Full workflow with calibration
 
```r
library(TwoTankModel)
 
# 1. Load your data
data <- load_data("daily_rainfall.csv", "daily_discharge.csv")
 
# 2. Calibrate (Monte Carlo with Latin Hypercube Sampling)
cal <- calibrate_montecarlo(data$times, data$precip, data$Q_obs, n_samples = 5000)
 
# 3. Uncertainty analysis
unc <- extract_uncertainty(cal, data$times, data$precip, data$Q_obs)
 
# 4. Performance report
calc_all_metrics(data$Q_obs, cal$best_sim$Q_total)
 
# 5. Generate all plots
plots <- plot_all(data$dates, data$precip, data$Q_obs, cal, unc)
```
 
### Option 3: Run the guided workflow script
 
```r
# From the package directory:
source("inst/examples/run_model.R")
```
 
Edit the **USER SETTINGS** at the top of `run_model.R` to point to your data.
 
## Input Data Format
 
**Rainfall CSV** (required):
```csv
date,P
2024-01-01,0.0
2024-01-02,5.3
2024-01-03,12.1
```
 
**Discharge CSV** (for calibration):
```csv
date,Q
2024-01-01,0.5
2024-01-02,0.8
2024-01-03,2.1
```
 
> **Tip:** If your discharge is in m³/s, use `m3s_to_mmday(Q, area_km2)` to convert, or set `area_km2` in `load_data()` for automatic conversion.
 
## Package Structure
 
```
TwoTankModel/
├── R/
│   ├── 01_model.R          # Core: run_two_tank(), print_model_summary()
│   ├── 02_metrics.R        # NSE, KGE, LogNSE, RMSE, PBIAS, calc_all_metrics()
│   ├── 03_calibration.R    # calibrate_montecarlo()
│   ├── 04_uncertainty.R    # extract_uncertainty()
│   ├── 05_analysis.R       # extract_metrics(), monthly_summary(),
│   │                       # check_mass_balance(), compare_watersheds()
│   ├── 06_plots.R          # 12 plot functions + plot_all()
│   └── 07_utils.R          # load_data(), generate_daily_rainfall(),
│                            # m3s_to_mmday(), export_results()
├── inst/examples/
│   ├── run_model.R          # Full guided workflow
│   └── quick_start.R        # 30-second demo
├── tests/
│   └── test_twotank.R       # Unit tests
├── DESCRIPTION
├── NAMESPACE
├── LICENSE
└── README.md
```
 
## All Functions
 
### Model
 
| Function | Description |
|----------|-------------|
| `run_two_tank(k1, k2, k3, times, precip)` | Run the two-tank simulation |
| `print_model_summary(sim)` | Pretty-print simulation results |
 
### Performance Metrics
 
| Function | Description |
|----------|-------------|
| `calc_nse(obs, sim)` | Nash-Sutcliffe Efficiency |
| `calc_kge(obs, sim)` | Kling-Gupta Efficiency |
| `calc_lognse(obs, sim)` | Log-NSE (baseflow focus) |
| `calc_rmse(obs, sim)` | Root Mean Square Error |
| `calc_pbias(obs, sim)` | Percent Bias |
| `calc_all_metrics(obs, sim)` | All 5 metrics with ratings table |
 
### Calibration & Uncertainty
 
| Function | Description |
|----------|-------------|
| `calibrate_montecarlo(times, precip, Q_obs)` | Monte Carlo calibration with LHS |
| `extract_uncertainty(cal, times, precip, Q_obs)` | Behavioural sets + 90% prediction bounds |
 
### Analysis
 
| Function | Description |
|----------|-------------|
| `extract_metrics(sim, dates)` | Peak Q, volume, BFI, runoff coefficient |
| `monthly_summary(sim, dates)` | Monthly water balance table |
| `check_mass_balance(sim)` | Verify P = Q + ΔS |
| `compare_watersheds(times, precip, params_A, params_B)` | Side-by-side watershed comparison |
 
### Plots (12 types)
 
| Function | Description |
|----------|-------------|
| `plot_all(...)` | Generate all plots + save to PDF/PNG |
| `plot_rainfall(dates, precip)` | Daily hyetograph |
| `plot_hydrograph(dates, Q_obs, sim)` | Observed vs simulated |
| `plot_uncertainty_envelope(dates, Q_obs, sim, unc)` | With 90% prediction bounds |
| `plot_scatter(Q_obs, sim)` | 1:1 scatter plot |
| `plot_components(dates, sim)` | Stacked Q1 + Q2 |
| `plot_storage(dates, sim)` | Tank storage dynamics |
| `plot_dotty(cal)` | Parameter identifiability |
| `plot_nse_histogram(cal)` | NSE distribution |
| `plot_param_correlation(unc, "k1", "k2")` | Parameter interactions |
| `plot_flow_duration(Q_obs, sim)` | Flow duration curve |
 
### Utilities
 
| Function | Description |
|----------|-------------|
| `load_data(rain_file, q_file)` | Smart CSV loader |
| `generate_daily_rainfall(n, dates)` | Synthetic rainfall |
| `m3s_to_mmday(Q, area_km2)` | Convert m³/s to mm/day |
| `mmday_to_m3s(Q, area_km2)` | Convert mm/day to m³/s |
| `export_results(cal, unc)` | Save results to CSV |
 
## Parameter Guide
 
| Watershed Type | k1 | k2 | k3 |
|---------------|-----|-----|-----|
| Urban (impervious) | 0.30 – 0.50 | 0.02 – 0.08 | 0.005 – 0.02 |
| Agricultural | 0.15 – 0.30 | 0.10 – 0.20 | 0.015 – 0.04 |
| Forest (natural) | 0.05 – 0.15 | 0.15 – 0.30 | 0.02 – 0.05 |
 
**Physical meaning:** Tank residence time ≈ 1/k days. For example, k1 = 0.3 means the upper tank drains 30% of its storage per day (~3 day residence time).
 
## Performance Ratings
 
| Rating | NSE | KGE | |PBIAS| |
|--------|-----|-----|----|
| Very Good | > 0.75 | > 0.75 | < 10% |
| Good | 0.65 – 0.75 | 0.50 – 0.75 | 10 – 15% |
| Acceptable | 0.50 – 0.65 | 0.0 – 0.50 | 15 – 25% |
| Poor | < 0.50 | < 0.0 | > 25% |
 
Based on Moriasi et al. (2007).
 
## Two-Watershed Comparison
 
Compare hydrologic response under identical rainfall:
 
```r
result <- compare_watersheds(
  times, precip,
  params_A = c(k1 = 0.35, k2 = 0.05, k3 = 0.01),   # Urban
  params_B = c(k1 = 0.08, k2 = 0.25, k3 = 0.04),   # Forest
  label_A  = "Urban Watershed",
  label_B  = "Forest Watershed",
  dates    = dates
)
```
 
## Building from Source
 
```r
# Generate roxygen documentation
devtools::document("path/to/TwoTankModel")
 
# Run tests
devtools::test("path/to/TwoTankModel")
 
# Check for issues
devtools::check("path/to/TwoTankModel")
 
# Build tar.gz
devtools::build("path/to/TwoTankModel")
```
 
## Dependencies
 
| Package | Purpose |
|---------|---------|
| `deSolve` | ODE solver for tank equations |
| `lhs` | Latin Hypercube Sampling |
| `ggplot2` | Publication-ready plots |
| `gridExtra` | Multi-panel plot layouts |
 
All dependencies install automatically.
 
## License
 
MIT
 
## Citation
 
If you use this package in your research, please cite:
 
```
Pokhrel, M. (2025). TwoTankModel: Two-Tank Conceptual Rainfall-Runoff Model
with Monte Carlo Calibration. R package version 1.0.0.
https://github.com/pokhrelmadan/TwoTankModel
```
