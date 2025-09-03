# -----------------------------------------------------------------------------
# Required Libraries
# -----------------------------------------------------------------------------
# install.packages(c("waveslim", "ggplot2", "reshape2", "tools"))
library(waveslim)
library(ggplot2)
library(reshape2)
library(tools) # Needed for file_path_sans_ext

#raw_data<-read.csv("C:/Users/deevanshi.walia/Desktop/Try 2.0/Group 2 Sub AA/With DT/12186 09-11-04 (pos)_.csv",check.names=FALSE, header=TRUE, sep=",")


soft_threshold <- function(data, threshold) {
  sign(data) * pmax(0, abs(data) - threshold)
}

process_and_validate_dwt <- function(raw_data, 
                                     wavelet = "haar", 
                                     levels = 5,
                                     threshold_fraction = 0.01,
                                     output_dir = "Denoised Heatmaps") {
  
  # --- Step 1: Data Preparation ---
  cat("⚙️  Preparing data...\n")
  raw_data[raw_data < 0] <- 0
  
  # NOTE: The line below removes the second column. This is specific to your
  # data format. If other files have a different structure, you may need to adjust this.
  #raw_data <- raw_data[, -2] 
  
  drift_times <- raw_data[[1]]
  numeric_matrix <- as.matrix(raw_data[, -1])
  retention_times_str <- colnames(numeric_matrix)
  retention_times_numeric <- as.numeric(gsub("^X", "", retention_times_str))
  
  # Create the output directory if it doesn't already exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # --- Step 2: Loop Through Decomposition Levels ---
  for (lvl in levels) {
    cat("▶️  Processing at decomposition level:", lvl, "\n")
    
    # ... (DWT, Thresholding, and Reconstruction are unchanged) ...
    dwt_result <- dwt.2d(numeric_matrix, wf = wavelet, J = lvl)
    for (name in names(dwt_result)) {
      if (grepl("LH|HL|HH", name)) {
        coeff_matrix <- dwt_result[[name]]
        adaptive_thresh <- 0.01 * max(abs(coeff_matrix), na.rm = TRUE)
        dwt_result[[name]] <- soft_threshold(coeff_matrix, adaptive_thresh)
      }
    }
    denoised_smoothed <- idwt.2d(dwt_result)
    
    # --- Data Cropping and Reshaping (Corrected for Robustness) ---
    min_rows <- min(nrow(numeric_matrix), nrow(denoised_smoothed))
    min_cols <- min(ncol(numeric_matrix), ncol(denoised_smoothed))
    
    numeric_matrix_cropped <- numeric_matrix[1:min_rows, 1:min_cols]
    denoised_smoothed_cropped <- denoised_smoothed[1:min_rows, 1:min_cols]
    drift_times_cropped <- drift_times[1:min_rows]
    #retention_times_cropped <- retention_times_numeric[1:min_cols]
    #retention_times_str_cropped <- retention_times_str[1:min_cols]
    retention_times_numeric_cropped <- retention_times_numeric[1:min_cols]
    
    # 1. Save the Denoised Matrix
    denoised_df <- as.data.frame(denoised_smoothed_cropped)
    # --- FIX ---: Use the cropped column names
    colnames(denoised_df) <- retention_times_numeric_cropped 
    denoised_df_to_save <- cbind(DriftTime = drift_times_cropped, denoised_df)
    denoised_csv_file <- file.path(output_dir, "denoised_matrix.csv")
    write.csv(denoised_df_to_save, file = denoised_csv_file, row.names = FALSE)
    cat("  - Saved denoised matrix to:", denoised_csv_file, "\n")
    
  
    # --- Plotting Section ---
    # Denoised Plot
    data_long_denoised <- reshape2::melt(denoised_smoothed_cropped)
    colnames(data_long_denoised) <- c("DriftTime", "RetentionTime", "Intensity")
    data_long_denoised$DriftTime <- drift_times_cropped
    #data_long_denoised$RetentionTime <- retention_times_cropped
    
    p_denoised <- ggplot(data_long_denoised, aes(x = DriftTime, y = RetentionTime, fill = Intensity)) +
      geom_tile() +
      scale_fill_viridis_c(option = "plasma") +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
      labs(title = paste("Denoised (Level", lvl, ", Thresh", threshold_fraction, ")"), x = "Drift Time (ms)", y = "Retention Time (s)")+
      theme_minimal(base_family = "sans") +
      theme(
        panel.grid.major = element_blank(),  # remove major grids
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),  # white background
        plot.background = element_rect(fill = "white"),   # white plot background
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        legend.text = element_text(color = "black")
  
      )
      
    
    out_file <- file.path(output_dir, paste0("GCIMS_denoised_level_", lvl, ".png"))
    ggsave(out_file, p_denoised, width = 10, height = 8, dpi = 600)
    cat("✅  Saved denoised map to:", out_file, "\n")
    
    # 2. Save the Difference Matrix (for Validation)
    difference_matrix <- numeric_matrix_cropped - denoised_smoothed_cropped
    diff_df <- as.data.frame(difference_matrix)
    # --- FIX ---: Use the cropped column names
    colnames(diff_df) <- retention_times_numeric_cropped
    diff_df_to_save <- cbind(DriftTime = drift_times_cropped, diff_df)
    diff_csv_file <- file.path(output_dir, "difference_matrix.csv")
    write.csv(diff_df_to_save, file = diff_csv_file, row.names = FALSE)
    cat("  - Saved difference matrix to:", diff_csv_file, "\n")
    
    
    # Difference Map Plot
    data_long_diff <- reshape2::melt(difference_matrix)
    colnames(data_long_diff) <- c("DriftTime", "RetentionTime_Index", "Intensity_Change")
    data_long_diff$DriftTime <- drift_times_cropped
    retention_times_numeric_long <- as.numeric(gsub("^X", "", data_long_diff$RetentionTime_Index))
    data_long_diff$RetentionTime_Index <- retention_times_numeric_long
    
    p_diff <- ggplot(data_long_diff, aes(x = DriftTime, y = RetentionTime_Index, fill = Intensity_Change)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "black", high = "red", midpoint = 0) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
      labs(title = paste("Removed Signal (Difference Map) - Level", lvl, "Thresh", threshold_fraction), x = "Drift Time(ms)", y = "Retention Time (s)")+
      theme_minimal(base_family = "sans") +
      theme(
        panel.grid.major = element_blank(),  # remove major grids
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),  # white background
        plot.background = element_rect(fill = "white"),   # white plot background
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        legend.text = element_text(color = "black")
      )
      
    
    diff_file <- file.path(output_dir, paste0("GCIMS_DIFFERENCE_level_", lvl, ".png"))
    ggsave(diff_file, p_diff, width = 10, height = 8, dpi = 600)
    cat("📊  Saved difference map to:", diff_file, "\n")
    
    # --- Statistics ---
    cor_val <- cor(as.vector(numeric_matrix_cropped), as.vector(denoised_smoothed_cropped))
    cat("🔍  Correlation to original at level", lvl, ":", round(cor_val, 5), "\n")
  }
}


# =============================================================================
# PART 2: MANAGER SCRIPT (This is the new batch processing loop)
# =============================================================================

# --- 1. Define Your Main Directories ---
# The folder where all your raw .csv files are located
input_data_dir <- "C:/Users/deevanshi.walia/Desktop/Try 2.0/Group 2 Sub AA/With DT" 
# A parent folder where all results will be saved
main_output_dir <- "C:/Users/deevanshi.walia/Desktop/Try 2.0/Group 2 Sub AA/Denoised Heatmaps"

# --- 2. Find All CSV Files to Process ---
# The `pattern` argument ensures we only get files ending in .csv
# `full.names = TRUE` gives us the full path to each file, which is essential
csv_files <- list.files(path = input_data_dir, pattern = "\\.csv$", full.names = TRUE)

cat("Found", length(csv_files), "CSV files to process.\n\n")

# --- 4. Loop Through Each File and Process It ---
for (i in 1:length(csv_files)) {
  
  current_file_path <- csv_files[i]
  
  # --- Create a unique output directory for this file's results ---
  file_name_base <- file_path_sans_ext(basename(current_file_path))
  specific_output_dir <- file.path(main_output_dir, file_name_base)
  
  # --- Print Progress Message ---
 
  cat("  Processing file", i, "of", length(csv_files), ":", basename(current_file_path), "\n")
 
  
  # --- Read the current CSV file ---
  current_raw_data <- read.csv(current_file_path, header = TRUE, sep = ",", check.names=FALSE)
  
  # --- Call the Worker Function ---
  # Pass the loaded data and the unique output directory to the function
  process_and_validate_dwt(
    raw_data = current_raw_data,
    levels = c(5), # Set your desired level(s)
    threshold_fraction = 0.01, # Set your optimal threshold
    output_dir = specific_output_dir
  )
  
  cat("\n") # Add a blank line for readability
}

cat("🎉 --- Batch processing complete. All files processed. ---\n")

