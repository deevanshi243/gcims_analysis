library(testthat)

test_that("peak helper returns peaks and properties", {
  source("scripts/00_setup.R")
  source("R/peak_detection.R")
  x <- c(0,1,0,1,0)
  out <- local_maxima_peaks_R(x, prominence = 0)
  expect_true(is.list(out))
  expect_true("peaks" %in% names(out))
})
