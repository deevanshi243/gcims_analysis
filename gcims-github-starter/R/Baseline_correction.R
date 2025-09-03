
#
# =============================================================================

# --- 1. Load Required Libraries ---
# install.packages(c("dplyr", "readr", "ggplot2", "reshape2", "tools"))
library(dplyr)
library(readr)
library(ggplot2)
library(reshape2)
library(tools)

# =============================================================================
# PART A: HELPER FUNCTION FOR BASELINE CORRECTION (MODIFIED)
# =============================================================================
# This function is now simplified as it works directly on the matrix.

#' Performs baseline correction on a numeric matrix.
#'
#' @param intensity_matrix A numeric matrix (rows=drift, cols=retention).
#' @param noise_rt_start_idx The starting COLUMN index for the noise region.
#' @param noise_rt_end_idx The ending COLUMN index for the noise region.
#' @param noise_dt_start_idx The starting ROW index for the noise region.
#' @param noise_dt_end_idx The ending ROW index for the noise region.
#' @return A new data frame with the baseline corrected.

baseline_correct_matrix <- function(intensity_matrix,
                                    noise_rt_start_idx, noise_rt_end_idx,
                                    noise_dt_start_idx, noise_dt_end_idx) {
  
  # --- Step 1: Isolate the Noise Region using matrix indexing ---
  # Ensure the indices are within the matrix bounds
  max_row <- nrow(intensity_matrix)
  max_col <- ncol(intensity_matrix)
  
  effective_dt_start <- min(noise_dt_start_idx, max_row)
  effective_dt_end   <- min(noise_dt_end_idx, max_row)
  effective_rt_start <- min(noise_rt_start_idx, max_col)
  effective_rt_end   <- min(noise_rt_end_idx, max_col)
  
  noise_block <- intensity_matrix[effective_dt_start:effective_dt_end,
                                  effective_rt_start:effective_rt_end]
  
  # --- Step 2: Calculate the Mean Noise Level ---
  mean_noise_level <- mean(noise_block, na.rm = TRUE)
  cat(paste0("    - Calculated Mean Noise Baseline: ", round(mean_noise_level, 4), "\n"))
  
  # --- Step 3: Subtract the Noise Level from the whole matrix ---
  corrected_matrix <- intensity_matrix - mean_noise_level
  
  # Set any resulting negative values to 0 (optional but good practice)
  corrected_matrix[corrected_matrix < 0] <- 0
  
  return(corrected_matrix)
}


# =============================================================================
# PART B: MAIN PROCESSING FUNCTION
# =============================================================================

process_denoised_files <- function(input_dir, output_dir) {
  
  # --- Setup ---
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  # Find the denoised CSV files
  csv_files <- list.files(path = input_dir, pattern = "_processed\\.csv$", full.names = TRUE, recursive = TRUE)
  if (length(csv_files) == 0) stop("No denoised CSV files found. Check input directory and file naming.")
  cat("Found", length(csv_files), "denoised files to process.\n\n")
  
  # --- Loop Through Each File ---
  for (file_path in csv_files) {
    
    # Get a clean name for the output file
    sample_name <- tools::file_path_sans_ext(basename(file_path)) # The folder name is the sample name
    cat("--- Processing Sample:", sample_name, "---\n")
    
    # --- a. Load Denoised Data ---
    # use check.names=FALSE to handle the numeric column headers correctly
    denoised_data <- read.csv(file_path, check.names = FALSE)
    
    # Separate metadata (DriftTime) from the intensity matrix
    drift_times <- denoised_data[, 1]
    intensity_matrix <- as.matrix(denoised_data[, -1])
    
    # --- b. Apply Baseline Correction ---
    cat("    - Applying baseline correction...\n")
    
    # DEFINE NOISE REGION using ROW and COLUMN INDICES
    # Based on our analysis of your heatmaps:
    
    # Retention Time (COLUMNS): Use the last part of the run.
    # The retention times are columns 210 in your data. Let's use 180-210.
    noise_col_start <- 180 
    noise_col_end   <- ncol(intensity_matrix) # Use all columns until the end
    
    # Drift Time (ROWS): Use the high drift time region.
    # Your drift time goes up to 20ms. Let's find the rows for 12ms to 20ms.
    # Assuming drift times are evenly spaced, we can approximate the row indices.
    # A safer way is to find the row indices corresponding to the values.
    noise_row_start <- which.min(abs(drift_times - 12)) # Find row index closest to 12ms
    noise_row_end   <- nrow(intensity_matrix)           # Go to the last row
    
    corrected_matrix <- baseline_correct_matrix(
      intensity_matrix = intensity_matrix,
      noise_rt_start_idx = noise_col_start,
      noise_rt_end_idx = noise_col_end,
      noise_dt_start_idx = noise_row_start,
      noise_dt_end_idx = noise_row_end
    )
    
    # --- c. Prepare Data for Plotting and Saving ---
    # Re-combine with DriftTime column for saving
    final_corrected_df <- as.data.frame(corrected_matrix)
    final_corrected_df <- cbind(DriftTime = drift_times, final_corrected_df)
    
    # Melt for plotting
    data_long <- reshape2::melt(final_corrected_df, id.vars = "DriftTime", variable.name = "RetentionTime", value.name = "Intensity")
    # Convert retention time (factor) to numeric for plotting
    data_long$RetentionTime <- as.numeric(as.character(data_long$RetentionTime))
    
    # --- d. Create and Save Final Heatmap ---
    p_final <- ggplot(data_long, aes(x = DriftTime, y = RetentionTime, fill = Intensity)) +
      geom_tile() +
      scale_fill_viridis_c(option = "plasma") +
      labs(
        title = paste(sample_name, "- Final Denoised & Baseline Corrected"),
        x = "Drift Time (ms)",
        y = "Retention Time (s)"
      ) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    
    # Create a specific output folder for this sample
    specific_output_dir <- file.path(output_dir, sample_name)
    if (!dir.exists(specific_output_dir)) dir.create(specific_output_dir, recursive = TRUE)
    
    output_plot_path <- file.path(specific_output_dir, "03_final_corrected_heatmap.png")
    ggsave(output_plot_path, p_final, width = 10, height = 7, dpi = 300)
    
    # --- e. Save the Final Corrected Data Matrix ---
    output_csv_path <- file.path(specific_output_dir, "04_final_corrected_matrix.csv")
    write.csv(final_corrected_df, output_csv_path, row.names = FALSE)
    
    cat("    - ✅ Saved final plot and CSV to:", specific_output_dir, "\n\n")
  }
  cat("🎉 --- Batch baseline correction complete! ---\n")
}
 


# =============================================================================
#                                 RUN THE SCRIPT
# =============================================================================

# ❗ --- YOU ONLY NEED TO CHANGE THESE TWO LINES --- ❗

# 1. Set the path to the folder containing your raw CSV files.
#    (Use forward slashes `/` or double backslashes `\\`)
input_directory <- "C:/Users/deevanshi.walia/Desktop/Try 2.0/Group 2 Sub AA/Aligned_and_RIP_Removed"

# 2. Set the path to the folder where you want to save the output images.
output_directory <- "C:/Users/deevanshi.walia/Desktop/Try 2.0/Group 2 Sub AA/Baseline_Estimated_Grp_2"


process_denoised_files(input_dir = input_directory, output_dir = output_directory)

