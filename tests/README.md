# TwoTankModel

A user-friendly two-tank conceptual rainfall–runoff model for daily streamflow simulation with Monte Carlo calibration in R.

**Key features:**
- Config-file-driven (edit `config.R`, run `main.R`)
- Native m³/s comparison (no unit conversion of observations needed)
- Automatic organized result folders per run
- Monte Carlo calibration with Latin Hypercube Sampling
- Uncertainty analysis with 90% prediction bounds
- Sensitivity analysis
- 12 publication-ready plots

## Quick Start

```bash
# 1. Clone the repository
git clone https://github.com/pokhrelmadan/TwoTankModel.git
cd TwoTankModel

# 2. Generate sample data (optional — skip if you have real data)
Rscript generate_sample_data.R

# 3. Edit config.R with your settings (see below)

# 4. Run the model
Rscript main.R
```

Results appear in `results/run_01/`.

## Installation as an R Package

```r
install.packages("remotes", repos = "https://cloud.r-project.org")
remotes::install_github("pokhrelmadan/TwoTankModel")
library(TwoTankModel)
```

## Model Structure

```
     Precipitation P(t)  [mm/day]
           │
           ▼
  ┌────────────────────┐
  │    Upper Tank       │──► Q1(t) = k1·S1  (Surface Runoff)
  │    Storage: S1      │
  └────────┬───────────┘         Total: Q(t) = Q1(t) + Q2(t)
           │ k2·S1
           ▼
  ┌────────────────────┐
  │    Lower Tank       │──► Q2(t) = k3·S2  (Baseflow)
  │    Storage: S2      │
  └────────────────────┘
```

**Equations:**
- dS1/dt = P(t) - k1·S1 - k2·S1
- dS2/dt = k2·S1 - k3·S2
- Q(t) = k1·S1 + k3·S2

## Configuration (`config.R`)

All settings live in one file. The main options:

```r
# ── Input data ──
RAINFALL_FILE       <- "data/rainfall.csv"       # columns: date, P
DISCHARGE_FILE      <- "data/discharge.csv"      # columns: date, Q
DISCHARGE_UNITS     <- "m3s"                     # "m3s" or "mm_day"
CATCHMENT_AREA_KM2  <- 150                       # required if m3s

# ── Calibration ──
N_MC_SAMPLES        <- 5000                      # 1000 quick, 5000 standard
K1_BOUNDS           <- c(0.01, 0.80)             # surface runoff
K2_BOUNDS           <- c(0.01, 0.50)             # percolation
K3_BOUNDS           <- c(0.001, 0.15)            # baseflow
NSE_THRESHOLD       <- 0.5                       # for behavioural sets
OBJECTIVE           <- "NSE"                     # "NSE", "KGE", or "LogNSE"

# ── Sensitivity ──
RUN_SENSITIVITY     <- TRUE
SENSITIVITY_PERTURB <- 0.20                      # ±20%

# ── Output ──
RUN_NAME            <- "run_01"                  # subfolder name
RESULTS_DIR         <- "results"                 # parent folder
```

## How Units Work

The model works in mm/day internally, but for comparison with observed discharge:

- **If your gauge data is in m³/s** (typical): set `DISCHARGE_UNITS <- "m3s"` and provide `CATCHMENT_AREA_KM2`. The model converts its simulated Q to m³/s and compares there. **This is the recommended approach** — you keep the original gauge units for reporting.

- **If your gauge data is in mm/day**: set `DISCHARGE_UNITS <- "mm_day"` and leave `CATCHMENT_AREA_KM2 <- NULL`. Comparison happens in mm/day.

## Input Data Format

**`data/rainfall.csv`:**
```csv
date,P
2024-01-01,0.0
2024-01-02,5.3
2024-01-03,12.1
```

**`data/discharge.csv`** (in m³/s):
```csv
date,Q
2024-01-01,0.5
2024-01-02,0.8
2024-01-03,2.1
```

The loader auto-detects column names like `P`, `P_mm`, `precipitation`, `Q`, `Q_m3s`, `discharge`, etc.

## Output Structure

After running `main.R`, results appear in:

```
results/
└── run_01/
    ├── config_used.R             ← Copy of your config (reproducibility)
    ├── parameters.txt            ← Best k1, k2, k3 + metrics
    ├── plots/
    │   ├── run_01_results.pdf    ← All plots in one PDF
    │   ├── run_01_envelope.png   ← Hydrograph + uncertainty bounds
    │   ├── run_01_scatter.png    ← Obs vs Sim scatter
    │   ├── run_01_dotty.png      ← Parameter identifiability
    │   ├── run_01_nse_hist.png   ← NSE distribution
    │   └── ... (8 more PNGs)
    └── csv/
        ├── run_01_mc_samples.csv         ← All N_MC_SAMPLES parameter sets
        ├── run_01_behavioural.csv        ← Acceptable sets (NSE > threshold)
        ├── run_01_sensitivity.csv        ← Sensitivity results
        ├── run_01_best_simulation.csv    ← Daily P, Q_obs, Q_sim, S1, S2
        └── monthly_summary.csv           ← Monthly water balance
```

Run it again with a different `RUN_NAME` to create a new folder and compare scenarios.

## Package Structure

```
TwoTankModel/
├── config.R                   # EDIT THIS
├── main.R                     # RUN THIS
├── generate_sample_data.R     # Creates test data
│
├── R/                         # All package functions
│   ├── 01_model.R             # run_two_tank()
│   ├── 02_metrics.R           # NSE, KGE, LogNSE, RMSE, PBIAS
│   ├── 03_calibration.R       # calibrate_montecarlo()
│   ├── 04_uncertainty.R       # extract_uncertainty()
│   ├── 05_analysis.R          # extract_metrics, monthly_summary
│   ├── 06_plots.R             # 12 plots + plot_all()
│   ├── 07_utils.R             # load_data, save_parameters, etc.
│   └── 08_sensitivity.R       # sensitivity_analysis()
│
├── data/                      # Your input CSVs go here
├── results/                   # Auto-generated run folders
├── tests/                     # Unit tests
│
├── DESCRIPTION
├── NAMESPACE
├── LICENSE
└── README.md
```

## Parameter Guide

| Watershed Type | k1 | k2 | k3 |
|---------------|-----|-----|-----|
| Urban         | 0.30 – 0.50 | 0.02 – 0.08 | 0.005 – 0.02 |
| Agricultural  | 0.15 – 0.30 | 0.10 – 0.20 | 0.015 – 0.04 |
| Forest        | 0.05 – 0.15 | 0.15 – 0.30 | 0.02 – 0.05 |

Tank residence time ≈ 1/k days.

## Performance Ratings

| Rating | NSE | KGE | \|PBIAS\| |
|--------|-----|-----|----------|
| Very Good | > 0.75 | > 0.75 | < 10% |
| Good | 0.65 – 0.75 | 0.50 – 0.75 | 10 – 15% |
| Acceptable | 0.50 – 0.65 | 0.0 – 0.50 | 15 – 25% |
| Poor | < 0.50 | < 0.0 | > 25% |

(Moriasi et al., 2007)

## License

MIT

## Citation

```
Pokhrel, M. (2025). TwoTankModel: Two-Tank Conceptual Rainfall-Runoff Model
with Monte Carlo Calibration. R package version 1.0.0.
https://github.com/pokhrelmadan/TwoTankModel
```
