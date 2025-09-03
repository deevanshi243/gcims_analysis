# GC-IMS Data Analysis (R + Python via reticulate)

This repository contains my Master's thesis analysis for **gas chromatography–ion mobility spectrometry (GC‑IMS)** data, including preprocessing (smoothing, baseline estimation, alignment), peak detection, and machine‑learning/statistical modeling.

## Quick start

```bash
# clone (replace with your URL)
git clone https://github.com/deevanshi243/gcims-analysis.git
cd gcims-analysis

# install packages
install.packages(c(
  "reticulate",
  "tools",
  "dbscan",
  "ggplot2",
  "reshape2",
  "dplyr",
  "readr",
  "caret",
  "tidyr",
  "MASS",
  "pROC",
  "randomForest",
  "e1071",
  "xgboost",
  "GGally",
  "tibble"
))


# Python deps
python -m pip install -r py/requirements.txt

# Configure reticulate to use this Python in R
# # If using system Python:
use_python("/usr/bin/python3", required = TRUE)
#(Optional) create a Python virtual environment and install Python dependencies
python -m venv .venv
source .venv/bin/activate  # (Windows: .venv\Scripts\activate)
python -m pip install -r py/requirements.txt

# Run the setup script to configure reticulate
Rscript scripts/00_setup.R

# Run the 1st main analysis pipeline
Rscript scripts/Smoothing_DWT_Adaptive_w_Padding.R
```

## Project structure

```
R/                 # Reusable R functions (wrappers, utilities)
py/                # Minimal Python helpers (used via reticulate)
scripts/           # Reproducible R scripts (entry points)
notebooks/         # R Markdown / Quarto exploratory notebooks
data/              # (Not versioned) raw/interim/processed data placeholders
figures/           # Exported plots
reports/           # Manuscript/slides outputs
tests/             # Unit tests (testthat)
.github/workflows/ # CI (GitHub Actions)
```
    

## License

Released under the MIT License (see `LICENSE`).
