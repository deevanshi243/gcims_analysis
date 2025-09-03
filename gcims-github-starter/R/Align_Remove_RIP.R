# -----------------------------------------------------------------------------
# Script for RIP-Based Alignment, Excision, and Visual Validation
# -----------------------------------------------------------------------------

# --- 1. Load Required Libraries ---
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(scales) # For pretty_breaks on plot axes

# =============================================================================
# 1. HELPER FUNCTIONS
# =============================================================================

#' Applies a circular shift to a vector.
apply_shift_r <- function(x, shift) {
  n <- length(x)
  if (n == 0 || shift == 0) return(x)
  shift <- shift %% n
  if (shift > 0) {
    return(c(x[(n - shift + 1):n], x[1:(n - shift)]))
  } else if (shift < 0) {
    shift <- abs(shift)
    return(c(x[(shift + 1):n], x[1:shift]))
  }
  return(x)
}

# =============================================================================
# 2. MAIN WORKFLOW FUNCTION
# =============================================================================

#' Aligns multiple GC-IMS samples using the RIP, then removes the RIP.
align_and_excise_rip <- function(input_dir,
                                 align_rip_window = c(5.0, 6.0),
                                 excise_rip_window = c(0, 5.7)) {
  
  # --- STEP 1: LOAD ALL DATA ---
  cat("STEP 1: Loading all denoised data files...\n")
  denoised_files <- list.files(input_dir, pattern = "_background_corrected\\.csv$", full.names = TRUE)
  if (length(denoised_files) == 0) stop("No denoised CSV files found.")
  
  all_data_list <- lapply(denoised_files, function(f) read_csv(f, show_col_types = FALSE))
  
  if (length(unique(sapply(all_data_list, function(df) paste(dim(df), collapse="x")))) > 1) {
    stop("Input files do not have consistent dimensions. Please check your denoising output.")
  }
  
  sample_names <- gsub("_background_corrected.csv", "", basename(denoised_files))
  drift_times <- all_data_list[[1]][[1]]
  retention_times <- as.numeric(colnames(all_data_list[[1]])[-1])
  
  intensity_array <- array(
    dim = c(length(all_data_list), nrow(all_data_list[[1]]), ncol(all_data_list[[1]]) - 1),
    dimnames = list(Sample = sample_names, DriftTime = drift_times, RetentionTime = retention_times)
  )
  for(i in 1:length(all_data_list)) { intensity_array[i,,] <- as.matrix(all_data_list[[i]][,-1]) }
  
  total_ion_matrix <- apply(intensity_array, c(1, 2), sum)
  
  # --- STEP 2: ALIGN SAMPLES USING THE RIP ---
  cat("\nSTEP 2: Aligning all samples based on the RIP...\n")
  align_cols <- which(drift_times >= align_rip_window[1] & drift_times <= align_rip_window[2])
  if (length(align_cols) < 5) stop("Alignment RIP window is too narrow.")
  
  master_reference_rip <- colMeans(total_ion_matrix[, align_cols], na.rm = TRUE)
  
  shifts <- sapply(1:nrow(total_ion_matrix), function(i) {
    ccf_result <- ccf(master_reference_rip, total_ion_matrix[i, align_cols], lag.max = length(align_cols) - 1, plot = FALSE)
    return(ccf_result$lag[which.max(ccf_result$acf)])
  })
  
  names(shifts) <- sample_names
  cat("    - Calculated alignment shifts (in data points):\n")
  print(shifts)
  
  aligned_array <- intensity_array
  for(i in 1:nrow(aligned_array)) {
    for(j in 1:dim(aligned_array)[3]) {
      aligned_array[i, , j] <- apply_shift_r(aligned_array[i, , j], shifts[i])
    }
  }
  
  # --- STEP 3: EXCISE THE RIP FROM ALIGNED DATA ---
  cat("\nSTEP 3: Excising RIP from all aligned samples...\n")
  rows_to_cut <- which(drift_times >= excise_rip_window[1] & drift_times <= excise_rip_window[2])
  if (length(rows_to_cut) == 0) {
    cat("    - No rows found in the excision window. Skipping removal.\n")
    final_array <- aligned_array
    final_drift_times <- drift_times
  } else {
    final_array <- aligned_array[, -rows_to_cut, ]
    final_drift_times_original <- drift_times[-rows_to_cut]
    cat("STEP 4: Re-calculating final continuous Drift Time axis with interpolation...\n")
    
    # Get original drift time step
    dt_step <- mean(diff(drift_times), na.rm = TRUE)
    
    # Interpolate linearly between the new min and max
    final_drift_times <- seq(
      from = min(final_drift_times_original),
      to   = max(final_drift_times_original),
      length.out = length(final_drift_times_original)
    )
  }
  
  cat("✅  Processing complete.\n")
  return(list(
    final_aligned_array = final_array, 
    final_drift_times = final_drift_times,
    original_retention_times = retention_times
  ))
}

# =============================================================================
# 3. BATCH WORKFLOW AND SAVING / PLOTTING
# =============================================================================

# --- Define Your Directories ---
input_dir <- "C:/Users/deevanshi.walia/Desktop/Try 2.0/Covid_New Grouping Analysis/Subtracted_AA_Grp_2"
output_dir <- "C:/Users/deevanshi.walia/Desktop/Try 2.0/Covid_New Grouping Analysis/Aligned_and_RIP_Removed_Grp_2"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- Run the full workflow ---
processed_data <- align_and_excise_rip(
  input_dir = input_dir,
  align_rip_window = c(5.0, 6.0),
  excise_rip_window = c(0, 5.7)
)

# --- Save final files AND generate plots for each sample ---
cat("\nSaving final processed files and generating validation heatmaps...\n")
final_array <- processed_data$final_aligned_array
final_dt <- processed_data$final_drift_times
retention_times <- processed_data$original_retention_times

for (i in 1:dim(final_array)[1]) {
  sample_name <- dimnames(final_array)$Sample[i]
  cat("  - Processing sample:", sample_name, "\n")
  
  # --- Reconstruct the 2D data frame for this sample ---
  sample_matrix <- final_array[i, , ]
  sample_df <- as.data.frame(sample_matrix)
  colnames(sample_df) <- retention_times
  final_df_to_save <- cbind(DriftTime = final_dt, sample_df)
  
  # --- 1. Save the processed CSV ---
  output_csv_path <- file.path(output_dir, paste0(sample_name, "_processed.csv"))
  write_csv(final_df_to_save, output_csv_path)
  
  # --- 2. Generate and Save the Validation Heatmap ---
  # Melt data into long format for ggplot
  plot_df <- final_df_to_save %>%
    pivot_longer(cols = -DriftTime, names_to = "RetentionTime", values_to = "Intensity") %>%
    mutate(RetentionTime = as.numeric(RetentionTime))
  
  p <- ggplot(plot_df, aes(x = DriftTime, y = RetentionTime, fill = Intensity)) +
    geom_tile() +
    scale_fill_viridis_c(option = "plasma", na.value = "grey80") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    labs(
      title = paste("Aligned & RIP-Removed -", sample_name),
      x = "Drift Time (ms)",
      y = "Retention Time (s)"
    ) +
    theme_minimal(base_family = "sans") +
    theme(
      panel.grid.major = element_blank(),  # remove major grids
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),  # white background
      plot.background = element_rect(fill = "white")
    )
  
  # Save the plot
  output_png_path <- file.path(output_dir, paste0(sample_name, "_heatmap.png"))
  ggsave(
    filename = output_png_path,
    plot = p,
    width = 10,
    height = 8,
    dpi = 150 # Use a lower DPI for faster saving during batch processing
  )
}

cat("\n🎉 --- All samples have been aligned, RIP-removed, saved, and plotted. ---\n")

