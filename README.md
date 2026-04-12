# TwoTankModel
 
A user-friendly two-tank conceptual rainfall–runoff model for daily streamflow simulation with Monte Carlo calibration in R.
 
## Installation
 
```r
install.packages("remotes", repos = "https://cloud.r-project.org")
remotes::install_github("pokhrelmadan/TwoTankModel")
```
 
## Quick Start (3 commands)
 
```bash
# 1. Create a config template in your current folder
Rscript -e 'TwoTankModel::create_config_template()'
 
# 2. Edit config.R — set your data file paths and run name
 
# 3. Run the complete pipeline
Rscript -e 'TwoTankModel::run_model("config.R")'
```
 
That's it. Results land in `results/<RUN_NAME>/` with all plots, CSVs, and parameters.
 
## Running from the Command Line
 
After installing the package, you can run the complete pipeline entirely from the terminal without opening R:
 
```bash
# Step 1: Navigate to your project folder
cd /path/to/my_project
 
# Step 2: Create a config template (first time only)
Rscript -e 'TwoTankModel::create_config_template()'
 
# Step 3: Edit config.R in any text editor
nano config.R          # or vim, VS Code, etc.
 
# Step 4: Put your data in the data/ folder
#   data/rainfall.csv   (columns: date, P)
#   data/discharge.csv  (columns: date, Q)
 
# Step 5: Run the model
Rscript -e 'TwoTankModel::run_model("config.R")'
```
 
**Running with a custom config name:**
 
```bash
Rscript -e 'TwoTankModel::run_model("my_watershed_config.R")'
```
 
**Running multiple scenarios back-to-back:**
 
```bash
Rscript -e 'TwoTankModel::run_model("urban_config.R")'
Rscript -e 'TwoTankModel::run_model("forest_config.R")'
```
 
Each run uses a different `RUN_NAME` in its config and creates a separate subfolder under `results/`, so you can compare scenarios without overwriting previous results.
 
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
 
**Governing equations:**
 
| Equation | Description |
|----------|-------------|
| dS1/dt = P(t) - k1·S1 - k2·S1 | Upper tank water balance |
| dS2/dt = k2·S1 - k3·S2 | Lower tank water balance |
| Q(t) = k1·S1 + k3·S2 | Total discharge |
 
## Configuration (`config.R`)
 
Edit these settings in your config file:
 
```r
# Input data
RAINFALL_FILE      <- "data/rainfall.csv"      # columns: date, P
DISCHARGE_FILE     <- "data/discharge.csv"     # columns: date, Q
DISCHARGE_UNITS    <- "m3s"                    # "m3s" or "mm_day"
CATCHMENT_AREA_KM2 <- 150                      # required if m3s
DATE_FORMAT        <- "auto"                   # or "%m/%d/%Y", "%Y-%m-%d"
 
# Calibration
N_MC_SAMPLES       <- 5000                     # 1000 quick / 5000 standard
K1_BOUNDS          <- c(0.01, 0.80)
K2_BOUNDS          <- c(0.01, 0.50)
K3_BOUNDS          <- c(0.001, 0.15)
OBJECTIVE          <- "NSE"                    # "NSE", "KGE", "LogNSE"
 
# Output
RUN_NAME           <- "run_01"
RESULTS_DIR        <- "results"
```
 
## How Units Work
 
The model works in mm/day internally, but comparison with observations happens in your chosen units:
 
- **Gauge data in m³/s** (typical): set `DISCHARGE_UNITS <- "m3s"` and provide `CATCHMENT_AREA_KM2`. The simulated Q gets converted to m³/s for comparison. **Recommended** — keeps original gauge units.
 
- **Data in mm/day**: set `DISCHARGE_UNITS <- "mm_day"`, leave `CATCHMENT_AREA_KM2 <- NULL`.
 
## Date Formats
 
The loader auto-detects common date formats by default (`DATE_FORMAT <- "auto"`). Supported formats:
 
| Format | Example |
|--------|---------|
| `%Y-%m-%d` | 2024-01-15 |
| `%m/%d/%Y` | 1/15/2024 (US) |
| `%d/%m/%Y` | 15/1/2024 (European) |
| `%Y/%m/%d` | 2024/01/15 |
| `%d-%m-%Y` | 15-01-2024 |
| `%m-%d-%Y` | 01-15-2024 |
 
If auto-detection fails or your dates are ambiguous, specify the format explicitly:
 
```r
DATE_FORMAT <- "%m/%d/%Y"
```
 
## Input Data Format
 
**`rainfall.csv`:**
```csv
date,P
2024-01-01,0.0
2024-01-02,5.3
2024-01-03,12.1
```
 
**`discharge.csv`** (m³/s):
```csv
date,Q
2024-01-01,0.5
2024-01-02,0.8
2024-01-03,2.1
```
 
Column names are auto-detected (`P`, `P_mm`, `precipitation` for rainfall; `Q`, `Q_m3s`, `discharge` for flow).
 
## Output Structure
 
After `run_model()`:
 
```
results/
└── run_01/
    ├── parameters.txt             ← best k1, k2, k3 + metrics
    ├── config_used.R              ← your config (reproducibility)
    ├── plots/
    │   ├── run_01_results.pdf     ← all plots
    │   ├── run_01_envelope.png    ← hydrograph + uncertainty
    │   ├── run_01_scatter.png     ← obs vs sim
    │   ├── run_01_dotty.png       ← identifiability
    │   └── … (9 more PNGs)
    └── csv/
        ├── run_01_mc_samples.csv          ← all parameter sets
        ├── run_01_behavioural.csv         ← acceptable sets
        ├── run_01_sensitivity.csv         ← sensitivity results
        ├── run_01_best_simulation.csv     ← daily time series
        └── monthly_summary.csv            ← monthly balance
```
 
## Using from R Interactively
 
```r
library(TwoTankModel)
 
# One-command pipeline
result <- run_model("config.R")
 
# Or call functions individually
data <- load_data("data/rainfall.csv", "data/discharge.csv",
                  discharge_units = "m3s", area_km2 = 150)
 
cal  <- calibrate_montecarlo(data$times, data$precip, data$Q_obs,
                             n_samples = 5000, area_km2 = 150)
 
unc  <- extract_uncertainty(cal, data$times, data$precip, data$Q_obs)
sens <- sensitivity_analysis(cal$best_params, data$times, data$precip,
                             data$Q_obs, area_km2 = 150)
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
