# GC-IMS Data Analysis (R + Python via reticulate)

This repository contains my Master's thesis analysis for **gas chromatography–ion mobility spectrometry (GC‑IMS)** data, including preprocessing (smoothing, baseline estimation, alignment), peak detection, and machine‑learning/statistical modeling.

## Quick start

```bash
# clone (replace with your URL)
git clone https://github.com/<your-username>/gcims-analysis.git
cd gcims-analysis

# (Recommended) restore R dependencies if you've used renv
# R -e 'install.packages("renv"); renv::restore()'

# Python deps
python -m pip install -r py/requirements.txt

# Configure reticulate to use this Python in R
# See scripts/00_setup.R for an example
```

## Project structure

```
R/                 # Reusable R functions (wrappers, utilities)
py/                # Minimal Python helpers (used via reticulate)
scripts/           # Reproducible R scripts (entry points)
notebooks/         # R Markdown / Quarto exploratory notebooks
data/              # (Not versioned) raw/interim/processed data placeholders
models/            # Trained model artifacts (usually ignored or kept small)
figures/           # Exported plots
reports/           # Manuscript/slides outputs
tests/             # Unit tests (testthat)
.github/workflows/ # CI (GitHub Actions)
```

## Reproduce the pipeline (example)

In R:

```r
source("scripts/00_setup.R")     # set Python path for reticulate, load packages
source("scripts/01_preprocess.R")# smoothing, baseline, alignment
source("scripts/02_model.R")     # training + evaluation
```

## Data

- Do **not** commit large or sensitive data. Place files locally under `data/`.
- Provide a small, anonymized sample dataset and a data dictionary if possible; otherwise document where to obtain the data.

## License

Released under the MIT License (see `LICENSE`).
