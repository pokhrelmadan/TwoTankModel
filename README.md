# TwoTankModel

A user-friendly two-tank conceptual rainfall–runoff model for daily streamflow simulation with Monte Carlo calibration.

## Model Structure

```
     Precipitation P(t)  [mm/day]
           │
           ▼
  ┌────────────────────┐
  │    Upper Tank       │──► Q1(t) = k1·S1  (Surface Runoff)
  │  (Surface & Inter-  │
  │       flow)         │          Total: Q(t) = Q1(t) + Q2(t)
  └────────┬───────────┘
           │ k2·S1  (Percolation)
           ▼
  ┌────────────────────┐
  │    Lower Tank       │──► Q2(t) = k3·S2  (Baseflow)
  │    (Baseflow)       │
  └────────────────────┘
```

## Quick Start (< 1 minute)

```r
# 1. Set working directory to the package folder
setwd("path/to/TwoTankModel")

# 2. Run the quick start demo
source("inst/examples/quick_start.R")
```

That's it! Check your folder for plots and results.

## Full Workflow (with your own data)

```r
setwd("path/to/TwoTankModel")
source("inst/examples/run_model.R")
```

Edit the **USER SETTINGS** at the top of `run_model.R`:
1. Set `use_own_data <- TRUE`
2. Point to your CSV files
3. Run the script

### Your CSV format

**Rainfall file** (required):
```csv
date,P
2024-01-01,0.0
2024-01-02,5.3
2024-01-03,12.1
```

**Discharge file** (for calibration):
```csv
date,Q
2024-01-01,0.5
2024-01-02,0.8
2024-01-03,2.1
```

If your discharge is in m³/s instead of mm/day, set `catchment_area <- 150` (your area in km²) and the script converts automatically.

## Package Structure

```
TwoTankModel/
├── R/
│   ├── 01_model.R          Core model engine
│   │     └── run_two_tank()          — Run a simulation
│   │     └── print_model_summary()   — Pretty-print results
│   │
│   ├── 02_metrics.R        Performance statistics
│   │     └── calc_nse()              — Nash-Sutcliffe Efficiency
│   │     └── calc_kge()              — Kling-Gupta Efficiency
│   │     └── calc_lognse()           — Log-NSE (baseflow focus)
│   │     └── calc_rmse()             — Root Mean Square Error
│   │     └── calc_pbias()            — Percent Bias
│   │     └── calc_all_metrics()      — All 5 at once with ratings
│   │
│   ├── 03_calibration.R    Monte Carlo calibration
│   │     └── calibrate_montecarlo()  — LHS sampling + evaluation
│   │
│   ├── 04_uncertainty.R    Uncertainty analysis
│   │     └── extract_uncertainty()   — Behavioural sets + envelope
│   │
│   ├── 05_analysis.R       Hydrograph analysis
│   │     └── extract_metrics()       — Peak Q, volume, BFI
│   │     └── monthly_summary()       — Monthly water balance
│   │     └── check_mass_balance()    — P = Q + ΔS verification
│   │     └── compare_watersheds()    — Side-by-side comparison
│   │
│   ├── 06_plots.R          12 diagnostic plots
│   │     └── plot_all()              — Generate everything
│   │     └── plot_rainfall()         — Hyetograph
│   │     └── plot_hydrograph()       — Obs vs Sim
│   │     └── plot_uncertainty_envelope()
│   │     └── plot_scatter()          — 1:1 plot
│   │     └── plot_components()       — Q1 + Q2 stacked
│   │     └── plot_storage()          — S1, S2 dynamics
│   │     └── plot_dotty()            — Parameter identifiability
│   │     └── plot_nse_histogram()    — NSE distribution
│   │     └── plot_param_correlation()— Parameter interactions
│   │     └── plot_flow_duration()    — FDC
│   │
│   └── 07_utils.R          Helpers
│         └── load_data()             — Smart CSV loader
│         └── generate_daily_rainfall()
│         └── m3s_to_mmday()          — Unit conversion
│         └── mmday_to_m3s()          — Unit conversion
│         └── export_results()        — Save to CSV
│
├── inst/examples/
│   ├── run_model.R          Full guided workflow (start here!)
│   └── quick_start.R        30-second demo
│
├── tests/
│   └── test_twotank.R       Unit tests
│
├── DESCRIPTION              Package metadata
├── NAMESPACE                Exports
├── LICENSE                  MIT
└── README.md                This file
```

## Parameters Guide

| Parameter | Description | Typical Range | Urban | Forest |
|-----------|-------------|---------------|-------|--------|
| k1 | Surface runoff [1/day] | 0.01 – 0.80 | 0.30 – 0.50 | 0.05 – 0.15 |
| k2 | Percolation [1/day] | 0.01 – 0.50 | 0.02 – 0.08 | 0.15 – 0.30 |
| k3 | Baseflow [1/day] | 0.001 – 0.15 | 0.005 – 0.02 | 0.02 – 0.05 |

**Physical meaning:**
- Tank residence time ≈ 1/k days
- k1 = 0.3 means upper tank drains 30% of storage per day (~3 day residence)
- k3 = 0.02 means lower tank drains 2% per day (~50 day residence)

## Performance Rating Guide

| Metric | Very Good | Good | Acceptable | Poor |
|--------|-----------|------|------------|------|
| NSE | > 0.75 | 0.65–0.75 | 0.50–0.65 | < 0.50 |
| KGE | > 0.75 | 0.50–0.75 | 0.0–0.50 | < 0.0 |
| PBIAS | < ±10% | ±10–15% | ±15–25% | > ±25% |

(Based on Moriasi et al., 2007)

## Installing as an R Package

```r
# Generate roxygen documentation
devtools::document("path/to/TwoTankModel")

# Check for issues
devtools::check("path/to/TwoTankModel")

# Install
devtools::install("path/to/TwoTankModel")

# Then use anywhere:
library(TwoTankModel)
sim <- run_two_tank(0.3, 0.1, 0.02, 0:364, precip)
```

## Dependencies

- `deSolve` — ODE solver
- `lhs` — Latin Hypercube Sampling
- `ggplot2` — Plotting
- `gridExtra` — Multi-panel layouts

All install automatically when needed.

## License

MIT
