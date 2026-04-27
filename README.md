# TwoTankModel
 
A two-tank conceptual rainfall–runoff model for daily streamflow simulation with nonlinear storage–discharge relationships, Monte Carlo calibration (LHS), uncertainty analysis, and sensitivity analysis in R. 
## Installation
 
```r
install.packages("remotes", repos = "https://cloud.r-project.org")
remotes::install_github("pokhrelmadan/TwoTankModel")
```
 
## Quick Start
 
```bash
# 1. Create a config template
Rscript -e 'TwoTankModel::create_config_template()'
 
# 2. Edit config.R with your data paths and settings
 
# 3. Run the complete pipeline
Rscript -e 'TwoTankModel::run_model("config.R")'
```
 
Results are saved in `results/<RUN_NAME>/` with plots, CSVs, and calibrated parameters.
 
## Running from the Command Line
 
```bash
cd /path/to/my_project
 
# First time — create config and put data in data/ folder
Rscript -e 'TwoTankModel::create_config_template()'
nano config.R
 
# Run the model
Rscript -e 'TwoTankModel::run_model("config.R")'
 
# Run multiple scenarios
Rscript -e 'TwoTankModel::run_model("config_urban.R")'
Rscript -e 'TwoTankModel::run_model("config_forest.R")'
```
 
Each run creates a separate results folder based on `RUN_NAME` in the config file.
 
## Model Structure
 
```
     Precipitation P(t)  [mm/day]
           │
           ▼
  ┌────────────────────┐
  │    Upper Tank       │──► Q1 = k1·S1^b1   (Surface Runoff, nonlinear)
  │    Storage: S1      │──► ET = k4·S1       (Evapotranspiration)
  └────────┬───────────┘
           │ k2·S1        Total: Q(t) = Q1(t) + Q2(t)
           │ (Percolation)
           ▼
  ┌────────────────────┐
  │    Lower Tank       │──► Q2 = k3·S2       (Baseflow, linear)
  │    Storage: S2      │
  └────────────────────┘
```
 
### Governing Equations
 
```
dS1/dt = P(t) - k1·S1^b1 - k2·S1 - k4·S1
dS2/dt = k2·S1 - k3·S2
Q(t)   = k1·S1^b1 + k3·S2
ET(t)  = k4·S1
```
 
When b1 = 1 (default), the model reduces to a standard linear two-tank model. When b1 > 1, the surface runoff responds nonlinearly to storage, producing sharper peaks during wet periods and lower flows during dry periods (Botter et al., 2009).
 
### Parameters
 
| Parameter | Description | Unit | Typical Range |
|-----------|-------------|------|---------------|
| k1 | Surface runoff coefficient | varies with b1 | 0.001 – 1.5 |
| k2 | Percolation coefficient | 1/day | 0.001 – 0.30 |
| k3 | Baseflow coefficient | 1/day | 0.001 – 0.05 |
| k4 | Evapotranspiration coefficient | 1/day | 0 – 0.20 |
| b1 | Nonlinear exponent | dimensionless | 1.0 – 5.0 |
 
## Configuration (`config.R`)
 
```r
# ── Input data ──
RAINFALL_FILE      <- "data/rainfall.csv"
DISCHARGE_FILE     <- "data/discharge.csv"
DISCHARGE_UNITS    <- "m3s"                    # "m3s" or "mm_day"
CATCHMENT_AREA_KM2 <- 150                      # required if m3s
DATE_FORMAT        <- "auto"                   # or "%m/%d/%Y", "%Y-%m-%d"
 
# ── Calibration ──
N_MC_SAMPLES       <- 10000
K1_BOUNDS          <- c(0.001, 1.50)
K2_BOUNDS          <- c(0.001, 0.30)
K3_BOUNDS          <- c(0.001, 0.05)
K4_BOUNDS          <- c(0.005, 0.20)           # set c(0, 0) to disable ET
B1_BOUNDS          <- c(1.0, 3.0)              # set c(1, 1) for linear mode
OBJECTIVE          <- "NSE"                    # "NSE", "KGE", or "LogNSE"
CAL_MONTHS         <- NULL                     # NULL=all, c(4:10)=Apr-Oct
 
# ── Output ──
RUN_NAME           <- "run_01"
RESULTS_DIR        <- "results"
```
 
### Calibration Period
 
Use `CAL_MONTHS` to calibrate on specific months only. This is useful for catchments where the model cannot represent certain processes (e.g., snowmelt in winter):
 
```r
CAL_MONTHS <- NULL         # Full year (default)
CAL_MONTHS <- c(4:10)      # April through October only
CAL_MONTHS <- c(6, 7, 8, 9)  # Summer only
```
 
The model still runs on the full record, but NSE is computed only on the selected months.
 
## How Units Work
 
- **Gauge data in m³/s** (typical): set `DISCHARGE_UNITS <- "m3s"` and provide `CATCHMENT_AREA_KM2`. The simulated Q is converted to m³/s for comparison.
- **Data in mm/day**: set `DISCHARGE_UNITS <- "mm_day"`, leave `CATCHMENT_AREA_KM2 <- NULL`.
## Date Formats
 
The loader auto-detects common formats with `DATE_FORMAT <- "auto"`:
 
| Format | Example |
|--------|---------|
| `%Y-%m-%d` | 2024-01-15 |
| `%m/%d/%Y` | 1/15/2024 |
| `%d/%m/%Y` | 15/1/2024 |
| `%Y/%m/%d` | 2024/01/15 |
 
## Input Data Format
 
**`rainfall.csv`:**
```csv
date,P
2024-01-01,0.0
2024-01-02,5.3
```
 
**`discharge.csv`** (m³/s):
```csv
date,Q
2024-01-01,0.5
2024-01-02,0.8
```
 
Column names are auto-detected (`P`, `P_mm`, `precipitation`; `Q`, `Q_m3s`, `discharge`, etc.).
 
## Output Structure
 
```
results/
└── run_01/
    ├── parameters.txt
    ├── config_used.R
    ├── plots/
    │   ├── run_01_results.pdf
    │   ├── run_01_hydrograph.png
    │   ├── run_01_envelope.png
    │   ├── run_01_scatter.png
    │   ├── run_01_dotty.png
    │   ├── run_01_components.png
    │   └── ... (up to 18 PNGs)
    └── csv/
        ├── run_01_mc_samples.csv
        ├── run_01_behavioural.csv
        ├── run_01_sensitivity.csv
        ├── run_01_best_simulation.csv
        └── monthly_summary.csv
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
                             n_samples = 10000, area_km2 = 150,
                             b1_range = c(1.0, 3.0),
                             k4_range = c(0.005, 0.20))
 
unc  <- extract_uncertainty(cal, data$times, data$precip, data$Q_obs)
sens <- sensitivity_analysis(cal$best_params, data$times, data$precip,
                             data$Q_obs, area_km2 = 150)
```
 
## Performance Ratings
 
| Rating | NSE | KGE | |PBIAS| |
|--------|-----|-----|---------|
| Very Good | > 0.75 | > 0.75 | < 10% |
| Good | 0.65 – 0.75 | 0.50 – 0.75 | 10 – 15% |
| Acceptable | 0.50 – 0.65 | 0.0 – 0.50 | 15 – 25% |
| Poor | < 0.50 | < 0.0 | > 25% |
 
Based on Moriasi et al. (2007).
 
## Package Structure
 
```
TwoTankModel/
├── R/
│   ├── 01_model.R             # run_two_tank() with Q1=k1·S1^b1
│   ├── 02_metrics.R           # NSE, KGE, LogNSE, RMSE, PBIAS
│   ├── 03_calibration.R       # Monte Carlo + LHS, parallel, CAL_MONTHS
│   ├── 04_uncertainty.R       # Behavioural sets, 90% prediction bounds
│   ├── 05_analysis.R          # Hydrograph metrics, monthly balance
│   ├── 06_plots.R             # 12–18 diagnostic plots
│   ├── 07_utils.R             # Data loading, unit conversion
│   ├── 08_sensitivity.R       # OAT perturbation analysis
│   └── 09_run.R               # run_model(), create_config_template()
├── DESCRIPTION
├── NAMESPACE
├── LICENSE                     # MIT
└── README.md
```
 
## References
 
### Tank Model
 
- Sugawara, M. and Funiyuki, F. (1956). A method of revision of the river discharge by means of a rainfall model. *Collection of Research Papers about Forecasting Hydrologic Variables*. Symp. Darcy, Dijon, IAHS Publication No. 51, pp. 71–76.
- Sugawara, M. (1961). Automatic calibration of the tank model. *Hydrological Sciences Bulletin*, 24(3), 375–388.
- Sugawara, M., Watanabe, I., Ozaki, E., and Katsuyame, Y. (1984). *Reference Manual for the Tank Model*. National Research Center for Disaster Prevention, Tokyo, Japan.
- Sugawara, M., Watanabe, E., Ozaki, E., and Katsuyama, Y. (1984). *Tank Model with Snow Component*. National Research Center for Disaster Prevention, Science and Technology Agency, Japan.
- Lee, Y.H., Singh, V.P., et al. (2020). A review of Tank Model and its applicability to various Korean catchment conditions. *Water*, 12(12), 3588. https://doi.org/10.3390/w12123588
### Nonlinear Storage–Discharge Relationships
 
- Botter, G., Porporato, A., Rodriguez-Iturbe, I., and Rinaldo, A. (2009). Nonlinear storage-discharge relations and catchment streamflow regimes. *Water Resources Research*, 45, W10427. https://doi.org/10.1029/2008WR007658
- Wittenberg, H. (1999). Baseflow recession and recharge as nonlinear storage processes. *Hydrological Processes*, 13, 715–726.
- Gan, R. and Luo, Y. (2013). Using the nonlinear aquifer storage–discharge relationship to simulate the base flow of glacier- and snowmelt-dominated basins in northwest China. *Hydrology and Earth System Sciences*, 17, 3577–3586. https://doi.org/10.5194/hess-17-3577-2013
### Model Evaluation Metrics
 
- Nash, J.E. and Sutcliffe, J.V. (1970). River flow forecasting through conceptual models: Part I — A discussion of principles. *Journal of Hydrology*, 10(3), 282–290. https://doi.org/10.1016/0022-1694(70)90255-6
- Moriasi, D.N., Arnold, J.G., Van Liew, M.W., Bingner, R.L., Harmel, R.D., and Veith, T.L. (2007). Model evaluation guidelines for systematic quantification of accuracy in watershed simulations. *Transactions of the ASABE*, 50(3), 885–900.
- Gupta, H.V., Kling, H., Yilmaz, K.K., and Martinez, G.F. (2009). Decomposition of the mean squared error and NSE performance criteria: Implications for improving hydrological modelling. *Journal of Hydrology*, 377(1–2), 80–91. https://doi.org/10.1016/j.jhydrol.2009.08.003
### Monte Carlo Calibration and Latin Hypercube Sampling
 
- McKay, M.D., Beckman, R.J., and Conover, W.J. (1979). A comparison of three methods for selecting values of input variables in the analysis of output from a computer code. *Technometrics*, 21(2), 239–245.
- Beven, K. and Binley, A. (1992). The future of distributed models: Model calibration and uncertainty prediction. *Hydrological Processes*, 6(3), 279–298. https://doi.org/10.1002/hyp.3360060305
## License
 
MIT
 
## Citation
 
```
Pokhrel, M. (2025). TwoTankModel: Two-Tank Conceptual Rainfall-Runoff Model
with Nonlinear Storage-Discharge Relationships and Monte Carlo Calibration.
R package version 1.0.0.
https://github.com/pokhrelmadan/TwoTankModel
```
