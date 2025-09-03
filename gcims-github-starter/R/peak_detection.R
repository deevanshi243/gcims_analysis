#' Local maxima peak detection (wrapper)
#'
#' @param x numeric vector representing a 1D signal
#' @param ... passed to scipy.signal::find_peaks via Python helper
#' @return list with `peaks` and `properties`
#' @examples
#' set.seed(1); x <- c(rnorm(100), 5, rnorm(100))
#' local_maxima_peaks_R(x, prominence = 1)
local_maxima_peaks_R <- function(x, ...) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Please install the 'reticulate' package")
  }
  reticulate::source_python("py/peak_detection.py")
  local_maxima_peaks(x, ...)
}
