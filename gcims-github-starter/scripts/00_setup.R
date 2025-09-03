# Configure reticulate to use your Python (already has packages from py/requirements.txt)
if (!requireNamespace("reticulate", quietly = TRUE)) install.packages("reticulate")
py <- Sys.which("python")
reticulate::use_python(py, required = TRUE)

# Load common R packages
pkgs <- c("tidyverse", "data.table", "matrixStats", "ggplot2", "testthat")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install)

message("Using Python: ", reticulate::py_config()$python)
