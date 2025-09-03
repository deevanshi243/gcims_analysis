library(data.table)
library(waveslim)
library(viridis)
library(ggplot2)
library(reshape2)

# Soft threshold function
threshold_array <- function(arr, threshold) {
  if (is.null(arr) || !is.numeric(arr)) return(arr)
  thresh_val <- threshold * max(abs(arr), na.rm = TRUE)
  return(ifelse(abs(arr) > thresh_val, arr, 0))
}

# Main processing function
process_gc_ims_and_save <- function(raw_data, wavelet = "haar", levels = c(1, 3, 4, 5), threshold = 0.01, output_dir = "C:/Users/deevanshi.walia/Desktop/Try 2.0/Try for DWT") {
  
  # Clean data: remove empty/unnamed columns, convert to numeric
  #raw_data <- raw_data[, !is.na(colnames(raw_data)) & colnames(raw_data) != ""]
  ##raw_data <- as.data.frame(raw_data)
  #numeric_data <- as.data.frame(lapply(raw_data[, -1], function(x) as.numeric(as.character(x))))
  #numeric_data[is.na(raw_data)] <- 0
  #numeric_matrix <- as.matrix(numeric_data)
  numeric_matrix <- as.matrix(raw_data[, -1])
  numeric_matrix[numeric_matrix < 0] <- 0
  
  # Create output dir if it doesn't exist
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  for (lvl in levels) {
    cat("▶️ Processing at level:", lvl, "\n")
    
    # 2D wavelet decomposition
    dwt_result <- dwt.2d(numeric_matrix, wf = wavelet, J = lvl)
    
    # Apply threshold to detail coefficients at each level
    for (j in 1:lvl) {
      dwt_result$d[[j]]$LH <- threshold_array(dwt_result$d[[j]]$LH, threshold)
      dwt_result$d[[j]]$HL <- threshold_array(dwt_result$d[[j]]$HL, threshold)
      dwt_result$d[[j]]$HH <- threshold_array(dwt_result$d[[j]]$HH, threshold)
    }
    
    # Reconstruct denoised matrix
    denoised_smoothed <- idwt.2d(dwt_result)
    drift_times <- raw_data[[1]]
    
    # Prepare for heatmap: melt long-format dataframe
    data_long <- reshape2::melt(denoised_smoothed)
    data_long$Var1 <- drift_times
    colnames(data_long) <- c("DriftTime", "RetentionTime", "Intensity")
   
    
    # Create heatmap
    p <- ggplot(data_long, aes(x = DriftTime, y = RetentionTime, fill = Intensity)) +
      geom_tile() +
      scale_fill_viridis_c(option = "plasma", na.value = "grey80") +
      scale_x_continuous(breaks = seq(1,20.48,by=0.2),minor_breaks = seq(1, 20.48, by = 0.1))+
      scale_y_continuous(breaks=seq(1,210,by=5), minor_breaks = seq(1, 213, by = 0.5))+
      theme_minimal(base_family = "sans") +
      theme( 
        plot.title = element_text(color = "white", size = 12, hjust = 0.5),
        axis.title.x = element_text(color = "white"),
        axis.title.y = element_text(color = "white"),
        axis.text.x = element_text(color = "white", angle = 90, hjust = 1, size = 6),
        axis.text.y = element_text(color = "white", size = 6),
        plot.background = element_rect(fill = "black", color = "white"),
        panel.background = element_rect(fill = "black", color = "white"),
        legend.title = element_text(color = "white"),
        legend.text  = element_text(color = "white"),
        legend.background = element_rect(fill = "black", color = NA),
        legend.key = element_rect(fill = "white", color = NA),
        panel.grid = element_blank()
      )+
      labs(title = paste("Denoised Smoothed GC-IMS (Level",lvl, ")"), x = "DriftTime", y = "RetentionTime")
    
    # Save heatmap image
    out_file <- file.path(output_dir, paste0("GCIMS_denoised_level_", lvl, ".png"))
    ggsave(out_file, p, width = 12, height = 10, dpi = 600)
    
    cat("✅ Saved heatmap at level", lvl, "to", out_file, "\n")
  }
}

# Example usage:
raw_data <- "C:/Users/deevanshi.walia/Desktop/Try 2.0/Try for DWT/1714 09-26-42 (pos)_.csv"
raw_data <- read.csv("1714 09-26-42 (pos)_.csv", header = TRUE,sep=",")
raw_data<-raw_data[,-1]
process_gc_ims_and_save(raw_data)

######check that the variance retained is a value close to 1 

original_var <- var(as.vector(numeric_matrix))
denoised_var <- var(as.vector(denoised_smoothed))

variance_retained <- denoised_var / original_var
print(variance_retained)

####SNR
 signal <- numeric_matrix[numeric_matrix > 0.01]
 noise <- numeric_matrix[numeric_matrix <= 0.01]
 
   snr_before <- mean(signal) / sd(noise)
 
   denoised_signal <- denoised_smoothed[denoised_smoothed > 0.01]
 denoised_noise <- denoised_smoothed[denoised_smoothed <= 0.01]
 
   snr_after <- mean(denoised_signal) / sd(denoised_noise)
 
  snr_before; snr_after
##[1] 0.7401071
 ## [1] 0.7419239

  
  
  # Find common dimensions
  common_nrow <- min(nrow(numeric_matrix), nrow(denoised_smoothed))
  common_ncol <- min(ncol(numeric_matrix), ncol(denoised_smoothed))
  
  # Subset both matrices to common dimensions
  numeric_matrix_cropped <- numeric_matrix[1:common_nrow, 1:common_ncol]
  denoised_smoothed_cropped <- denoised_smoothed[1:common_nrow, 1:common_ncol]
  
  # Now compute correlation
  cor(as.vector(numeric_matrix_cropped), as.vector(denoised_smoothed_cropped))
  
  
   