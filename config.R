################################################################################
#  CONFIG.R — USER CONFIGURATION FILE
#  ────────────────────────────────────
#  Edit this file to match your data and setup. Do not edit any other files.
#  Then run:  Rscript main.R
################################################################################


# ═══════════════════════════════════════════════════════════════════════════════
# 1. INPUT DATA
# ═══════════════════════════════════════════════════════════════════════════════

# Path to your daily rainfall CSV
#   Required columns: date, P (or P_mm, precipitation)
#   Example row:      2024-01-01, 5.3
RAINFALL_FILE <- "data/rainfall.csv"

# Path to your daily observed discharge CSV
#   Required columns: date, Q (or Q_m3s, discharge, flow)
#   Discharge can be in m³/s (recommended) or mm/day
#   Example row:      2024-01-01, 2.5
DISCHARGE_FILE <- "data/discharge.csv"

# Units of observed discharge in the file above
#   Use "m3s"    if discharge is in m³/s  (recommended for most gauge data)
#   Use "mm_day" if discharge is already in mm/day
DISCHARGE_UNITS <- "m3s"

# Catchment area in km² (REQUIRED if DISCHARGE_UNITS = "m3s")
#   This is used to convert simulated Q (mm/day) to m³/s for comparison.
#   Set to NULL if DISCHARGE_UNITS = "mm_day".
CATCHMENT_AREA_KM2 <- 150

# Date format in your CSV files
#   Use "auto" to auto-detect (tries common formats)
#   Or specify explicitly:
#     "%Y-%m-%d"  for 2024-01-15
#     "%m/%d/%Y"  for 1/15/2024   (US format)
#     "%d/%m/%Y"  for 15/1/2024   (European format)
DATE_FORMAT <- "auto"

# ═══════════════════════════════════════════════════════════════════════════════
# 2. CALIBRATION SETTINGS
# ═══════════════════════════════════════════════════════════════════════════════

# Number of Monte Carlo samples
#   Use 1000  for quick testing
#   Use 5000  for standard analysis (recommended)
#   Use 10000+ for publication
N_MC_SAMPLES <- 5000

# Parameter search bounds [1/day]
K1_BOUNDS <- c(0.01, 0.80)     # Surface runoff coefficient
K2_BOUNDS <- c(0.01, 0.50)     # Percolation coefficient
K3_BOUNDS <- c(0.001, 0.15)    # Baseflow coefficient

K4_BOUNDS <- c(0.001, 0.15)    # ET coefficient
# Nonlinear surface runoff exponent
#   Q1 = k1 * S1^b1
#   B1_BOUNDS <- c(1, 1)       # linear mode (default)
#   B1_BOUNDS <- c(1.0, 3.0)   # nonlinear mode (sharper peaks)
B1_BOUNDS <- c(1.0, 3.0)
# Threshold for behavioural (acceptable) parameter sets
#   Standard: NSE >= 0.5
NSE_THRESHOLD <- 0.5

# Objective function for calibration
#   Options: "NSE", "KGE", "LogNSE"
OBJECTIVE <- "NSE"


# ═══════════════════════════════════════════════════════════════════════════════
# 3. SENSITIVITY ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════════

# Enable sensitivity analysis?
RUN_SENSITIVITY <- TRUE

# Perturbation size (± fraction of parameter value)
#   0.20 means ±20% around each parameter
SENSITIVITY_PERTURB <- 0.20


# ═══════════════════════════════════════════════════════════════════════════════
# 4. OUTPUT SETTINGS
# ═══════════════════════════════════════════════════════════════════════════════

# Name of this run (creates a subfolder under results/)
#   Use a descriptive name; spaces become underscores
RUN_NAME <- "run_01"

# Parent results folder (each run creates a subfolder here)
RESULTS_DIR <- "results"

# What to save (set any to FALSE to skip)
SAVE_PLOTS        <- TRUE    # All diagnostic plots (PDF + PNGs)
SAVE_CSV_RESULTS  <- TRUE    # Calibration samples, best simulation
SAVE_PARAMETERS   <- TRUE    # Best parameters as a .txt file
SAVE_CONFIG_COPY  <- TRUE    # Copy of this config.R to the results folder


# ═══════════════════════════════════════════════════════════════════════════════
# 5. ADVANCED (usually leave as default)
# ═══════════════════════════════════════════════════════════════════════════════

# Random seed for reproducibility
RANDOM_SEED <- 42

# Initial tank storage (usually 0, but can warm up with a value)
INITIAL_S1 <- 0
INITIAL_S2 <- 0

# Print verbose messages during run?
VERBOSE <- TRUE


################################################################################
#  END OF CONFIG — DO NOT EDIT BELOW THIS LINE
################################################################################
