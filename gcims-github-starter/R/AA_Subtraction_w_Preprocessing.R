# =============================================================================
#
#   --- FINAL ROBUST WORKFLOW: DENOISE AA, AVERAGE, & SUBTRACT (v3) ---
#
# Description:
# This definitive script is designed to handle complex data formats, including
# multiple initial metadata columns. It provides a complete, robust pipeline.
#
#   1. Intelligently parses files to separate metadata from intensity data.
#   2. Validates each file to ensure it's a 2D matrix before processing.
#   3. Denoises each valid AA file using the Discrete Wavelet Transform (DWT).
#   4. Averages the denoised AA matrices to create a master background.
#   5. Loops through sample files, subtracts the background, and saves the
#      final, cleanly formatted matrix with a single DriftTime column.
#
# =============================================================================

# --- 1. Load Required Libraries ---
# install.packages(c("waveslim", "dplyr", "readr", "tools"))
library(waveslim)
library(dplyr)
library(readr)
library(tools)

# =============================================================================
# PART A: DWT HELPER FUNCTIONS (UNCHANGED)
# =============================================================================

soft_threshold <- function(data, threshold) {
  sign(data) * pmax(0, abs(data) - threshold)
}

pad_matrix_for_dwt <- function(mat, J) {
  target_multiple <- 2^J
  current_rows <- nrow(mat)
  target_rows <- ceiling(current_rows / target_multiple) * target_multiple
  rows_to_add <- target_rows - current_rows
  if (rows_to_add > 0) {
    padding_rows <- mat[1:rows_to_add, , drop = FALSE]
    padding_rows <- padding_rows[rev(seq_len(nrow(padding_rows))), , drop = FALSE]
    mat <- rbind(mat, padding_rows)
  }
  current_cols <- ncol(mat)
  target_cols <- ceiling(current_cols / target_multiple) * target_multiple
  cols_to_add <- target_cols - current_cols
  if (cols_to_add > 0) {
    padding_cols <- mat[, 1:cols_to_add, drop = FALSE]
    padding_cols <- padding_cols[, rev(seq_len(ncol(padding_cols))), drop = FALSE]
    mat <- cbind(mat, padding_cols)
  }
  return(mat)
}

# =============================================================================
# PART B: ROBUST DWT WORKER FUNCTION (HEAVILY MODIFIED)
# =============================================================================

denoise_single_matrix <- function(file_path, wavelet, levels, threshold_fraction, separator) {
  
  cat("  - Attempting to process:", basename(file_path), "\n")
  
  raw_data <- tryCatch({
    read.csv(file_path, check.names = FALSE, sep = separator) 
  }, error = function(e) {
    warning("    -> CRITICAL: Could not read file:", basename(file_path), "\n       Reason: ", e$message)
    return(NULL)
  })
  
  if (is.null(raw_data)) return(NULL)
  
  # --- INTELLIGENT METADATA DETECTION ---
  # Find the first column whose name is a pure number. That's our data start.
  # suppressWarnings is used because non-numeric names will become NA, which is what we want.
  numeric_headers <- suppressWarnings(as.numeric(colnames(raw_data)))
  first_data_col <- min(which(!is.na(numeric_headers)))
  
  if (is.infinite(first_data_col)) {
    warning(paste0("    -> SKIPPING: '", basename(file_path), "' has no numeric column headers."))
    return(NULL)
  }
  
  # --- ROBUST SLICING ---
  # Keep the first column as our primary DriftTime vector
  drift_time_vector <- raw_data[, 1]
  # The intensity matrix starts at the first numeric column header
  intensity_matrix <- as.matrix(raw_data[, first_data_col:ncol(raw_data)])
  
  intensity_matrix[intensity_matrix < 0] <- 0
  original_rows <- nrow(intensity_matrix)
  original_cols <- ncol(intensity_matrix)
  
  # --- DWT Processing ---
  padded_matrix <- pad_matrix_for_dwt(intensity_matrix, J = levels)
  dwt_result <- dwt.2d(padded_matrix, wf = wavelet, J = levels)
  
  for (name in names(dwt_result)) {
    if (grepl("LH|HL|HH", name)) {
      coeff_matrix <- dwt_result[[name]]
      adaptive_thresh <- threshold_fraction * max(abs(coeff_matrix), na.rm = TRUE)
      dwt_result[[name]] <- soft_threshold(coeff_matrix, adaptive_thresh)
    }
  }
  
  denoised_padded <- idwt.2d(dwt_result)
  denoised_matrix <- denoised_padded[1:original_rows, 1:original_cols, drop = FALSE]
  
  cat("    -> SUCCESS: Denoising complete for", basename(file_path), "\n")
  # Return a list containing both the denoised matrix AND the drift time vector
  return(list(drift_times = drift_time_vector, matrix = denoised_matrix))
}


# =============================================================================
# PART C: MAIN SCRIPT EXECUTION
# =============================================================================

# -----------------------------------------------------------------------------
# ❗ --- 1. CONFIGURE YOUR SCRIPT HERE --- ❗
# -----------------------------------------------------------------------------

air_dir <- "C:/Users/deevanshi.walia/Desktop/Try 2.0/AA"
sample_dir <- "C:/Users/deevanshi.walia/Desktop/Gruppe 1/Gruppe 1 neu/With DT"
output_dir <- "C:/Users/deevanshi.walia/Desktop/Try 2.0/Subtracted_AA_Grp_1"

dwt_wavelet <- "haar"
dwt_levels <- 5
dwt_threshold <- 0.05
file_separator <- ","
file_pattern <- "^_[0-9].*\\.csv$" # Finds files starting with a number

# -----------------------------------------------------------------------------
# ❗ --- SCRIPT EXECUTION STARTS HERE --- ❗
# -----------------------------------------------------------------------------

# --- STEP A: DENOISE AND AVERAGE AMBIENT AIR SPECTRA ---
cat("--- STEP A: Processing Ambient Air (AA) files ---\n")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

air_files <- list.files(air_dir, full.names = TRUE, pattern = file_pattern)
if (length(air_files) == 0) stop(paste("No data files matching pattern '", file_pattern, "' found in 'air_dir'."))

# The result is now a list of lists, e.g., list(list(drift_times, matrix), list(drift_times, matrix), ...)
denoised_air_results_list <- lapply(air_files, function(file) {
  denoise_single_matrix(file, dwt_wavelet, dwt_levels, dwt_threshold, file_separator)
})

denoised_air_results_list <- denoised_air_results_list[!sapply(denoised_air_results_list, is.null)]
if (length(denoised_air_results_list) == 0) stop("All AA files were skipped or failed processing. Check warnings.")

# Extract just the matrices for averaging
denoised_matrices_only <- lapply(denoised_air_results_list, `[[`, "matrix")

air_array <- simplify2array(denoised_matrices_only)
average_denoised_air <- apply(air_array, c(1, 2), mean)
cat("\n✅  Successfully created average background from", length(denoised_air_results_list), "valid AA files.\n")

# To save, grab the coordinates and headers from the first validly processed file
first_valid_result <- denoised_air_results_list[[1]]
drift_times_for_saving <- first_valid_result$drift_times
retention_time_headers <- colnames(first_valid_result$matrix)

average_air_df <- as.data.frame(average_denoised_air)
colnames(average_air_df) <- retention_time_headers
average_air_df <- cbind(DriftTime = drift_times_for_saving, average_air_df)
write.csv(average_air_df, file.path(output_dir, "MASTER_AVERAGE_AIR_BACKGROUND.csv"), row.names = FALSE)

file_pattern_samples<-"^[0-9].*\\.csv$"
# --- STEP B: SUBTRACT AVERAGE BACKGROUND FROM SAMPLES ---
cat("\n--- STEP B: Subtracting background from sample files ---\n")
sample_files <- list.files(sample_dir, full.names = TRUE, pattern = file_pattern_samples)
if (length(sample_files) == 0) stop(paste("No data files matching pattern '", file_pattern, "' found in 'sample_dir'."))

for (file_path in sample_files) {
  sample_name <- tools::file_path_sans_ext(basename(file_path))
  cat("  - Correcting:", sample_name, "\n")
  
  sample_data_full <- read.csv(file_path, check.names = FALSE, sep = file_separator)
  
  # Use the same logic to robustly slice the sample data
  numeric_headers_sample <- suppressWarnings(as.numeric(colnames(sample_data_full)))
  first_data_col_sample <- min(which(!is.na(numeric_headers_sample)))
  
  drift_times_vector <- sample_data_full[, 1]
  sample_intensity_matrix <- as.matrix(sample_data_full[, first_data_col_sample:ncol(sample_data_full)])
  
  if (!all(dim(sample_intensity_matrix) == dim(average_denoised_air))) {
    warning(paste("    -> SKIPPING:", sample_name, "due to dimension mismatch."))
    next
  }
  
  corrected_intensity_matrix <- sample_intensity_matrix - average_denoised_air
  corrected_intensity_matrix[corrected_intensity_matrix < 0] <- 0
  
  # Re-combine into the final, clean format
  corrected_df <- as.data.frame(corrected_intensity_matrix)
  colnames(corrected_df) <- colnames(average_denoised_air) # Use headers from the processed data
  corrected_df_final <- cbind(DriftTime = drift_times_vector, corrected_df)
  
  output_file_path <- file.path(output_dir, paste0(sample_name, "_background_corrected.csv"))
  write.csv(corrected_df_final, file = output_file_path, row.names = FALSE)
}

cat("\n🎉 --- Workflow complete! Corrected files are in:", output_dir, "---\n")

# =============================================================================
#
#   --- SCRIPT: GENERATE BEFORE & AFTER HEATMAPS FOR QC ---
#
# Description:
# This script creates visual comparisons for the background subtraction step.
# It finds all original sample files and their corresponding corrected
# files, and for each pair, it generates two heatmaps in a dedicated
# subfolder for easy comparison.
#
# =============================================================================

# --- 1. Load Required Libraries ---
# install.packages(c("ggplot2", "reshape2", "tools", "dplyr", "readr"))
library(ggplot2)
library(reshape2)
library(tools)
library(dplyr)
library(readr)

# =============================================================================
# PART A: HELPER FUNCTION TO CREATE A STANDARDIZED HEATMAP
# =============================================================================

#' Creates a high-quality heatmap from a GC-IMS data frame.
#'
#' @param data_df A data frame with the first column as DriftTime and the rest as intensities.
#' @param plot_title The title for the plot.
#' @param file_separator The separator used when reading the CSV (e.g., ",").
#' @return A ggplot object.

create_gcims_heatmap <- function(data_df, plot_title) {
  first_col_name <- colnames(data_df)[1]
  # Melt the data into a long format for ggplot
  # This converts the wide matrix into three columns: DriftTime, RetentionTime, Intensity
  data_long <- data_df %>%
    reshape2::melt(
      id.vars = first_col_name,
      variable.name = "RetentionTime",
      value.name = "Intensity"
    ) %>%
    # Ensure RetentionTime is treated as a numeric value for the axis
    mutate(RetentionTime = as.numeric(as.character(RetentionTime)))
  colnames(data_long)[1] <- "DriftTime"
  # Create the plot
  p <- ggplot(data_long, aes(x = DriftTime, y = RetentionTime, fill = Intensity)) +
    geom_tile() +
    scale_fill_viridis_c(option = "plasma") + # Use a visually appealing color scale
    scale_x_continuous(expand = c(0, 0)) +   # Remove padding on axes
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      title = plot_title,
      x = "Drift Time (ms)",
      y = "Retention Time (s)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      panel.grid = element_blank() # Remove grid lines for a cleaner look
    )
  
  return(p)
}


# =============================================================================
# PART B: MAIN SCRIPT EXECUTION
# =============================================================================

# -----------------------------------------------------------------------------
# ❗ --- 1. CONFIGURE YOUR SCRIPT HERE --- ❗
# -----------------------------------------------------------------------------

# --- Set Directories ---
# Folder containing the ORIGINAL sample CSV files (the "before" data)
original_samples_dir <- "C:/Users/deevanshi.walia/Desktop/Gruppe 1/Gruppe 1 neu/With DT"

# Folder containing the BACKGROUND-CORRECTED sample files (the "after" data)
corrected_samples_dir <- "C:/Users/deevanshi.walia/Desktop/Try 2.0/Subtracted_AA_Grp_1"

# A NEW parent folder where all the output plot subfolders will be created
output_plots_dir <- "C:/Users/deevanshi.walia/Desktop/Try 2.0/QC_Heatmap_Comparisons_Grp_1"

# --- Set File Reading Parameters ---
file_separator <- ","

# -----------------------------------------------------------------------------
# ❗ --- SCRIPT EXECUTION STARTS HERE --- ❗
# -----------------------------------------------------------------------------

cat("--- Starting Heatmap Generation Workflow ---\n")
dir.create(output_plots_dir, showWarnings = FALSE, recursive = TRUE)

# Find all the original sample files
original_files <- list.files(original_samples_dir, full.names = TRUE, pattern = "\\.csv$")
if (length(original_files) == 0) stop("No original sample files found in the specified directory.")

cat("Found", length(original_files), "original samples to process.\n\n")

# Loop through each original sample file
for (before_path in original_files) {
  
  # --- a. Construct file names and paths ---
  sample_name_base <- tools::file_path_sans_ext(basename(before_path))
  cat("-> Processing sample:", sample_name_base, "\n")
  
  # Construct the expected name of the corresponding "after" file
  after_filename <- paste0(sample_name_base, "_background_corrected.csv")
  after_path <- file.path(corrected_samples_dir, after_filename)
  
  # --- b. Check if both "before" and "after" files exist ---
  if (!file.exists(after_path)) {
    warning(paste("  - SKIPPING:", sample_name_base, "- Could not find corresponding corrected file at:", after_path))
    next # Skip to the next file in the loop
  }
  
  # --- c. Read both data files ---
  # Use a tryCatch block for safe reading
  data_before <- tryCatch(read.csv(before_path, sep = file_separator, check.names = FALSE), error = function(e) NULL)
  data_after <- tryCatch(read.csv(after_path, sep = file_separator, check.names = FALSE), error = function(e) NULL)
  
  if (is.null(data_before) || is.null(data_after)) {
    warning(paste("  - SKIPPING:", sample_name_base, "- Error reading one of the files."))
    next
  }
  
  # --- d. Create plots ---
  cat("  - Generating heatmaps...\n")
  plot_before <- create_gcims_heatmap(data_before, plot_title = paste(sample_name_base, "\n(Original Data)"))
  plot_after <- create_gcims_heatmap(data_after, plot_title = paste(sample_name_base, "\n(After Background Subtraction)"))
  
  # --- e. Create a unique subfolder and save the plots ---
  sample_output_subdir <- file.path(output_plots_dir, sample_name_base)
  dir.create(sample_output_subdir, showWarnings = FALSE, recursive = TRUE)
  
  path_before_plot <- file.path(sample_output_subdir, "01_heatmap_before_correction.png")
  path_after_plot <- file.path(sample_output_subdir, "02_heatmap_after_correction.png")
  
  ggsave(path_before_plot, plot_before, width = 10, height = 7, dpi = 300)
  ggsave(path_after_plot, plot_after, width = 10, height = 7, dpi = 300)
  
  cat("  - ✅ Saved plots to subfolder:", sample_name_base, "\n\n")
}

cat("🎉 --- Workflow complete! All heatmaps have been generated. ---\n")
