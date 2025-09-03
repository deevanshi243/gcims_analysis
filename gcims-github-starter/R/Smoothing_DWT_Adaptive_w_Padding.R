# =============================================================================
# PART 1: HELPER FUNCTION FOR PADDING
# =============================================================================
# This function pads a matrix by reflecting its edges, ensuring its dimensions
# are divisible by 2^J, which is required by dwt.2d.
pad_matrix_for_dwt <- function(mat, J) {
  target_multiple <- 2^J
  
  # --- Pad Rows ---
  current_rows <- nrow(mat)
  # Calculate the next multiple of target_multiple
  target_rows <- ceiling(current_rows / target_multiple) * target_multiple
  rows_to_add <- target_rows - current_rows
  
  if (rows_to_add > 0) {
    # Create the padding by reflecting the last few rows
    padding_rows <- mat[(current_rows - rows_to_add + 1):current_rows, , drop = FALSE]
    # Reverse the order of the reflected rows
    padding_rows <- padding_rows[rev(seq_len(nrow(padding_rows))), , drop = FALSE]
    mat <- rbind(mat, padding_rows)
  }
  
  # --- Pad Columns ---
  current_cols <- ncol(mat)
  target_cols <- ceiling(current_cols / target_multiple) * target_multiple
  cols_to_add <- target_cols - current_cols
  
  if (cols_to_add > 0) {
    # Create padding by reflecting the last few columns
    padding_cols <- mat[, (current_cols - cols_to_add + 1):current_cols, drop = FALSE]
    # Reverse the order of the reflected columns
    padding_cols <- padding_cols[, rev(seq_len(ncol(padding_cols))), drop = FALSE]
    mat <- cbind(mat, padding_cols)
  }
  
  return(mat)
}

# =============================================================================
# PART 2: THE MAIN WORKER FUNCTION (with padding)
# =============================================================================
process_and_validate_dwt <- function(raw_data, 
                                     wavelet = "haar", 
                                     levels = c(5),
                                     threshold_fraction = 0.01,
                                     output_dir = "./DWT_Results") {
  
  # --- Step 1: Data Preparation ---
  cat("⚙️  Preparing data...\n")
  raw_data[raw_data < 0] <- 0
  # raw_data <- raw_data[, -2] # You have this commented, which is fine.
  
  drift_times <- raw_data[[1]]
  numeric_matrix <- as.matrix(raw_data[, -1])
  retention_times_str <- colnames(numeric_matrix)
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # --- Step 2: Loop Through Decomposition Levels ---
  for (lvl in levels) {
    cat("▶️  Processing at decomposition level:", lvl, "\n")
    
    # --- NEW: Pad the matrix before DWT to prevent data loss ---
    cat("    - Padding matrix to be compatible with Level", lvl, "...\n")
    numeric_matrix_padded <- pad_matrix_for_dwt(numeric_matrix, J = lvl)
    
    # --- DWT on the PADDED matrix ---
    dwt_result <- dwt.2d(numeric_matrix_padded, wf = wavelet, J = lvl)
    
    # ... (Thresholding loop is unchanged) ...
    for (name in names(dwt_result)) {
      if (grepl("LH|HL|HH", name)) {
        coeff_matrix <- dwt_result[[name]]
        adaptive_thresh <- threshold_fraction * max(abs(coeff_matrix), na.rm = TRUE)
        dwt_result[[name]] <- soft_threshold(coeff_matrix, adaptive_thresh)
      }
    }
  }
}
    