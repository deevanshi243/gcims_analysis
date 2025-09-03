## =============================================================================
#
#   --- MASTER DATA PROCESSING PIPELINE: FROM HEATMAPS TO FEATURE MATRICES ---
#
# Workflow:
#   1. Sets up Python integration using 'reticulate'.
#   2. Calls an external Python script to perform 2D peak picking on all samples.
#   3. Aggregates all detected peaks from all samples.
#   4. Performs DBSCAN clustering to define features from peaks.
#   5. Automatically defines Regions of Interest (ROIs) from cluster boundaries.
#   6. Generates and saves four distinct feature matrices:
#      - Binary (Presence/Absence)
#      - Maximum Intensity
#      - Median Intensity
#      - ROI Integration (Peak Volume)
#   7. Saves ROI validation heatmaps for quality control.
#

# =============================================================================

# --- 1. LOAD REQUIRED LIBRARIES ---

library(reticulate)
library(tools)
library(dbscan)
library(ggplot2)
library(reshape2)
library(dplyr); library(readr); library(caret); library(tidyr); library(MASS)
library(pROC); library(ggplot2); library(randomForest); library(e1071)
library(xgboost); library(GGally); library(tibble)

# =============================================================================
# PART A: CONFIGURATION
# ============================================================================

# --- A1: Python and Script Configuration ---
python_executable_path <- "C:/Users/deevanshi.walia/Desktop/new/Scripts/python.exe"

# Define the full path to your Python script.
python_script_path <- "C:/Users/deevanshi.walia/Desktop/Try 2.0/Onset of Symptoms/peak_picker.py" # Assuming it's here now

# --- A quick check to make sure both files actually exist ---
if (!file.exists(python_executable_path)) stop(paste("PYTHON NOT FOUND at:", python_executable_path))
if (!file.exists(python_script_path)) stop(paste("PYTHON SCRIPT NOT FOUND at:", python_script_path))


# --- A2: Directory Configuration ---
# Parent directory containing your group folders
# These folders should contain your FINAL preprocessed heatmap CSVs in groups for the intial code (everything before ROI matrix generation).
#change the base dir to the folder with all the samples in one folder for the ROI matrices
base_data_dir <-"C:/Users/deevanshi.walia/Desktop/Try 2.0/Covid_New Grouping Analysis/Preprocessed_AA_All_Samples"
group_folders <- c("Days_1_3", "Days_4_7")
#C:/Users/deevanshi.walia/Desktop/Try 2.0/Covid_New Grouping Analysis/Preprocessed_AA_All_Samples
# Parent directory for all outputs
output_base_dir <- "C:/Users/deevanshi.walia/Desktop/Try 2.0/Onset of Symptoms/Final_Analysis_Output_Peaks_Clustering_9_0.05"

# Specific output subdirectories (will be created automatically)
peaks_output_dir <- file.path(output_base_dir, "Peak_Lists")
clustering_output_dir <- file.path(output_base_dir, "Clustering_Results")
feature_matrix_dir <- file.path(output_base_dir, "Feature_Matrices")
roi_heatmap_dir <- file.path(output_base_dir, "ROI_Heatmaps")

# Create all output directories
sapply(c(peaks_output_dir, clustering_output_dir, feature_matrix_dir, roi_heatmap_dir),
       dir.create, showWarnings = FALSE, recursive = TRUE)

# --- A3: Parameter Configuration ---
# Peak Picking Parameters (to be passed to the Python script)
PEAK_MIN_DISTANCE <- 9
PEAK_THRESHOLD_ABS <- 0.05

# DBSCAN Clustering Parameters
DBSCAN_EPS <- 0.08  # Epsilon: search radius (on SCALED data)
DBSCAN_MIN_PTS <- 9    # Minimum points to form a cluster

# ROI Integration Parameters
ROI_RT_PADDING <- 2.0   # Padding in seconds
ROI_DT_PADDING <- 0.1   # Padding in milliseconds



# PART B: PEAK PICKING (CALLING PYTHON FROM R)
# =============================================================================
# =============================================================================
# PART B: PEAK PICKING (CALLING PYTHON FROM R - ROBUST METHOD)
# =============================================================================
cat("--- PART B: RUNNING PYTHON PEAK PICKING SCRIPT ---\n")

# NOTE: We no longer use `source_python` for this script.
for (group in group_folders) {
  cat(paste0("  - Processing group: ", group, "\n"))
  
  input_group_dir <- file.path(base_data_dir, group)
  output_group_dir <- file.path(peaks_output_dir, group)
  
  # Now, we build the command using the FULL path to the python.exe
  command <- paste(
    shQuote(python_executable_path), # <-- Use the full path to the correct python
    shQuote(python_script_path),
    shQuote(input_group_dir),
    shQuote(output_group_dir),
    PEAK_MIN_DISTANCE,
    PEAK_THRESHOLD_ABS
  )
  cat("    - Executing Python command:\n     ", command, "\n")
  
  # The system() command will now use your fully-featured virtual environment.
  python_output <- system(command, intern = TRUE)
  
  # Print the captured output
  cat("    --- Python Script Output ---\n")
  cat(paste0("      ", python_output), sep = "\n")
  cat("    --- End of Python Output ---\n\n")
}
cat("✅  Peak picking complete for all groups.\n")



# =============================================================================
# PART C: AGGREGATE PEAKS AND PERFORM DBSCAN CLUSTERING
# =============================================================================
cat("\n--- PART C: AGGREGATING PEAKS AND RUNNING DBSCAN ---\n")

all_peak_files <- list.files(peaks_output_dir, recursive = TRUE, pattern = "^peaks_.*\\.csv$", full.names = TRUE)
if (length(all_peak_files) == 0) stop("FATAL ERROR: No peak list files were generated or found.")

all_peaks <- all_peak_files %>%
  lapply(function(f) {
    group_name <- basename(dirname(f))
    sample_name <- tools::file_path_sans_ext(basename(f)) %>% gsub("peaks_", "", .)
    read_csv(f, show_col_types = TRUE) %>% mutate(SampleName = sample_name, GroupName = group_name)
  }) %>% bind_rows()

peaks_for_clustering <- all_peaks %>% dplyr::select(RetentionTime, DriftTime)%>%scale()
dbscan_result <- dbscan(peaks_for_clustering, eps = DBSCAN_EPS, minPts = DBSCAN_MIN_PTS)
all_peaks$Cluster <- dbscan_result$cluster
write_csv(all_peaks, file.path(clustering_output_dir, "master_peak_list_with_clusters.csv"))
cat("✅  DBSCAN complete. Found", max(all_peaks$Cluster), "clusters.\n")

# --- Generate DBSCAN Visualization Plot ---
dbscan_plot_data <- as.data.frame(peaks_for_clustering)
dbscan_plot_data$Cluster <- as.factor(all_peaks$Cluster)
dbscan_scaled_plot <- ggplot(dbscan_plot_data, aes(x = DriftTime, y = RetentionTime, color = Cluster)) +
  geom_point(alpha = 0.6, size = 1.5) +
  labs(title = "DBSCAN Clustering Results", x="Drift Time", y="Retention Time") +
  theme_minimal() + guides(color = guide_legend(override.aes = list(size = 3)))
ggsave(file.path(clustering_output_dir, "DBSCAN_cluster_mapD.png"), dbscan_scaled_plot, width=10, height=10, dpi=1200, bg="white")
cat("✅  DBSCAN cluster map saved.\n")

# =============================================================================
# PART D: GENERATE ALL FOUR FEATURE MATRICES
# =============================================================================
cat("\n--- PART D: GENERATING FEATURE MATRICES ---\n")

master_metadata <- all_peaks %>% distinct(SampleName, GroupName)
peaks_no_noise <- all_peaks %>% filter(Cluster != 0)

# --- 1, 2, 3: Binary, Max, Median Matrices ---
cat("  - Generating peak-list based matrices (Binary, Max, Median)...\n")
feature_matrix_binary <- peaks_no_noise %>% distinct(SampleName, Cluster) %>% mutate(Present=1) %>%
  pivot_wider(names_from=Cluster, values_from=Present, names_prefix="Cluster_") %>%
  right_join(master_metadata, by="SampleName") %>% mutate(across(starts_with("Cluster_"), ~replace_na(., 0)))
write_csv(feature_matrix_binary, file.path(feature_matrix_dir, "feature_matrix_BINARY.csv"))

feature_matrix_max <- peaks_no_noise %>% group_by(SampleName, Cluster) %>% summarise(Intensity=max(Intensity),.groups='drop') %>%
  pivot_wider(names_from=Cluster, values_from=Intensity, names_prefix="Cluster_") %>%
  right_join(master_metadata, by="SampleName") %>% mutate(across(starts_with("Cluster_"), ~replace_na(., 0)))
write_csv(feature_matrix_max, file.path(feature_matrix_dir, "feature_matrix_MAX_INTENSITY.csv"))

feature_matrix_median <- peaks_no_noise %>% group_by(SampleName, Cluster) %>% summarise(Intensity=median(Intensity),.groups='drop') %>%
  pivot_wider(names_from=Cluster, values_from=Intensity, names_prefix="Cluster_") %>%
  right_join(master_metadata, by="SampleName") %>% mutate(across(starts_with("Cluster_"), ~replace_na(., 0)))
write_csv(feature_matrix_median, file.path(feature_matrix_dir, "feature_matrix_MEDIAN_INTENSITY.csv"))

# --- 4. ROI Integration Matrix ---
cat("  - Generating ROI Integration Matrix...\n")
roi_definitions <- peaks_no_noise %>% group_by(Cluster) %>%
  summarise(RT_min=min(RetentionTime)-ROI_RT_PADDING, RT_max=max(RetentionTime)+ROI_RT_PADDING,
            DT_min=min(DriftTime)-ROI_DT_PADDING, DT_max=max(DriftTime)+ROI_DT_PADDING, .groups='drop')
feature_rois <- setNames(split(roi_definitions, seq(nrow(roi_definitions))), paste0("Cluster_", roi_definitions$Cluster))
roi_plot_df <- bind_rows(feature_rois, .id = "FeatureName")

all_sample_features_roi <- list()
all_heatmap_files <- list.files(path=base_data_dir, recursive=TRUE, pattern="\\.csv$", full.names=TRUE)

for (file_path in all_heatmap_files) {
  sample_name <- tools::file_path_sans_ext(basename(file_path)) # Assumes clean names
  cat(paste0("\r    - Integrating ROIs for sample: ", sample_name, "      ")); flush.console()
  
  data_matrix_df <- read_csv(file_path, show_col_types = FALSE)
  drift_times <- data_matrix_df[[1]]
  retention_times <- as.numeric(colnames(data_matrix_df)[-1])
  
  current_sample_values <- c()
  for (feature_name in names(feature_rois)) {
    roi_df <- feature_rois[[feature_name]]
    row_indices <- which(drift_times >= roi_df$DT_min & drift_times <= roi_df$DT_max)
    col_indices <- which(retention_times >= roi_df$RT_min & retention_times <= roi_df$RT_max)
    
    if (length(row_indices) > 0 && length(col_indices) > 0) {
      integrated_value <- sum(as.matrix(data_matrix_df[row_indices, col_indices + 1, drop = FALSE]), na.rm = TRUE)
    } else {
      integrated_value <- 0
    }
    current_sample_values[feature_name] <- integrated_value
  }
  all_sample_features_roi[[sample_name]] <- current_sample_values
  
  # --- Generate and save the ROI validation heatmap ---
  id_col_name <- colnames(data_matrix_df)[1]
  data_long <- reshape2::melt(data_matrix_df, id.vars = id_col_name, variable.name = "RetentionTime", value.name = "Intensity") %>%
    mutate(RetentionTime = as.numeric(as.character(RetentionTime)))
  colnames(data_long)[1] <- "DriftTime"
  p_heatmap <- ggplot(data_long, aes(x = DriftTime, y = RetentionTime, fill = Intensity)) + geom_tile() + scale_fill_viridis_c(option = "plasma")
  p_validation <- p_heatmap + geom_rect(data = roi_plot_df, aes(xmin = DT_min, xmax = DT_max, ymin = RT_min, ymax = RT_max), color = "red", fill = NA, size = 0.5, inherit.aes = FALSE) +
    labs(title = paste("ROI Validation for Sample:", sample_name)) + theme_minimal()
  ggsave(file.path(roi_heatmap_dir, paste0(sample_name, "_ROI_map.png")), p_validation, width = 10, height = 7, dpi = 150)
}
cat("\n")

feature_matrix_roi_data_only <- bind_rows(all_sample_features_roi, .id = "SampleName")
feature_matrix_roi <- right_join(master_metadata, feature_matrix_roi_data_only, by = "SampleName") %>%
  mutate(across(starts_with("Cluster_"), ~replace_na(., 0))) %>%
  dplyr::select(SampleName, GroupName, everything())
write_csv(feature_matrix_roi, file.path(feature_matrix_dir, "feature_matrix_ROI_INTEGRATION.csv"))

cat("All four feature matrices generated and saved.\n")
cat("\n --- MASTER DATA PROCESSING WORKFLOW COMPLETE! ---\n")



