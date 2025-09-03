source("scripts/00_setup.R")
source("R/peak_detection.R")

# ---- placeholders ----
# 1) Load raw GC-IMS data (user-provided) from data/raw/
# 2) Apply smoothing / filtering
# 3) Baseline estimation
# 4) Alignment across samples
# 5) Peak detection using local_maxima_peaks_R for 1D slices or adapt to 2D

# Example (toy 1D signal):
set.seed(42)
x <- stats::filter(rnorm(200), rep(1/5, 5))
x <- as.numeric(x)
res <- local_maxima_peaks_R(x, prominence = 0.5)
print(res$peaks)
