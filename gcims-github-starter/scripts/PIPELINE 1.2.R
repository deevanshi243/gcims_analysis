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
# Make sure you have these installed.
# install.packages(c("reticulate", "dplyr", "readr", "tools", "dbscan", "ggplot2", "reshape2"))
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
# =============================================================================

# --- A1: Python and Script Configuration ---
python_executable_path <- "C:/Users/deevanshi.walia/Desktop/new/Scripts/python.exe"

# Define the full path to your Python script.
python_script_path <- "C:/Users/deevanshi.walia/Desktop/Try 2.0/Onset of Symptoms/peak_picker.py" # Assuming it's here now

# --- A quick check to make sure both files actually exist ---
if (!file.exists(python_executable_path)) stop(paste("PYTHON NOT FOUND at:", python_executable_path))
if (!file.exists(python_script_path)) stop(paste("PYTHON SCRIPT NOT FOUND at:", python_script_path))


# --- A2: Directory Configuration ---
# Parent directory containing your group folders (e.g., "PCR_Positive", "PCR_Negative")
# These folders should contain your FINAL preprocessed heatmap CSVs.
base_data_dir <- "C:/Users/deevanshi.walia/Desktop/Try 2.0/PCR/Preprocessed_AA_All_Samples"
group_folders <- c("Mild_Low_Load (CT greater than 27)", "Severe_High_Load (CT_less_than_27)")

# Parent directory for all outputs
output_base_dir <- "C:/Users/deevanshi.walia/Desktop/Try 2.0/Organized_by_CT_Value_Severity_27CT/Final_Analysis_Output_Peaks_Clustering_9_0.05"

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
DBSCAN_EPS <- 0.08   # Epsilon: search radius (on SCALED data)
DBSCAN_MIN_PTS <- 9     # Minimum points to form a cluster

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

# Aggregate all generated peak lists
all_peak_files <- list.files(peaks_output_dir, recursive = TRUE, pattern = "^peaks_.*\\.csv$", full.names = TRUE)
all_peaks <- all_peak_files %>%
  lapply(function(f) {
    # Extract metadata from the file path
    group_name <- basename(dirname(f))
    sample_name <- tools::file_path_sans_ext(basename(f)) %>% gsub("peaks_", "", .)
    read_csv(f, show_col_types = FALSE) %>%
      mutate(SampleName = sample_name, GroupName = group_name)
  }) %>%
  bind_rows()

# Prepare data for clustering (scale coordinates)
peaks_for_clustering <- all_peaks %>%
  select(RetentionTime, DriftTime) %>%
  scale()

# Optional: Generate k-NN plot to validate 'eps'
kNNdistplot(peaks_for_clustering, k = DBSCAN_MIN_PTS)
abline(h = DBSCAN_EPS, lty = 2, col = "red")

# Run DBSCAN
dbscan_result <- dbscan(peaks_for_clustering, eps = DBSCAN_EPS, minPts = DBSCAN_MIN_PTS)

# Add cluster assignments back to the main peak data frame
all_peaks$Cluster <- dbscan_result$cluster

cat("✅  DBSCAN complete. Found", max(all_peaks$Cluster), "clusters.\n")
# Save the final clustered peak list
write_csv(all_peaks, file.path(clustering_output_dir, "master_peak_list_with_clusters.csv"))

# =============================================================================
# PART C.2: VISUALIZE THE DBSCAN CLUSTERING RESULTS (CORRECTED AXES)
# =============================================================================

cat("\n--- Generating DBSCAN Cluster Visualization Plot ---\n")

# --- Step 1: Prepare the data for plotting ---
plot_data <- all_peaks %>%
  mutate(Cluster_Factor = as.factor(Cluster))
noise_points <- plot_data %>% filter(Cluster == 0)
clustered_points <- plot_data %>% filter(Cluster != 0)

# --- NEW: Explicitly define the full range of your data ---
# This prevents ggplot from automatically zooming in and hiding outliers.
min_dt <- floor(min(all_peaks$DriftTime, na.rm = TRUE))
max_dt <- ceiling(max(all_peaks$DriftTime, na.rm = TRUE))
min_rt <- floor(min(all_peaks$RetentionTime, na.rm = TRUE))
max_rt <- ceiling(max(all_peaks$RetentionTime, na.rm = TRUE))

# --- Step 2: Create the plot using ggplot2 ---
dbscan_plot <- ggplot() +
  geom_point(data = noise_points, aes(x = DriftTime, y = RetentionTime),
             color = "grey80", size = 1.0, alpha = 0.5) +
  geom_point(data = clustered_points, aes(x = DriftTime, y = RetentionTime, color = Cluster_Factor),
             alpha = 0.8, size = 1.5) +
  
  # --- Formatting and Labels ---
  labs(
    title = "DBSCAN Clustering Results",
    subtitle = paste("eps =", DBSCAN_EPS, ", MinPts =", DBSCAN_MIN_PTS),
    x = "Drift Time (ms)",
    y = "Retention Time (s)",
    color = "Cluster"
  ) +
  
  # --- THIS IS THE KEY FIX ---
  # We manually set the limits for the axes to ensure all data is shown.
  scale_x_continuous(limits = c(min_dt, max_dt), breaks = scales::pretty_breaks(n = 8)) +
  scale_y_continuous(limits = c(min_rt, max_rt), breaks = scales::pretty_breaks(n = 10)) +
  
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 3)) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

# --- Step 3: Print the plot and save it ---
print(dbscan_plot)
plot_output_path <- file.path(clustering_output_dir, "DBSCAN_cluster_map.png")
ggsave(plot_output_path, plot = dbscan_plot, width = 12, height = 8, dpi = 300, bg = "white")
cat("✅  DBSCAN cluster map saved to:", plot_output_path, "\n")

# =============================================================================
# PART D: GENERATE ALL FOUR FEATURE MATRICES
# =============================================================================
cat("\n--- PART D: GENERATING FEATURE MATRICES ---\n")

# Filter out noise points for feature matrix creation
peaks_no_noise <- all_peaks %>% filter(Cluster != 0)

# --- 1. Binary (Presence/Absence) Matrix ---
cat("  - Generating Binary Matrix...\n")
feature_matrix_binary <- peaks_no_noise %>%
  distinct(SampleName, GroupName, Cluster) %>%
  mutate(Present = 1) %>%
  pivot_wider(
    id_cols = c(SampleName, GroupName),
    names_from = Cluster,
    values_from = Present,
    values_fill = 0,
    names_prefix = "Cluster_"
  )
write_csv(feature_matrix_binary, file.path(feature_matrix_dir, "feature_matrix_BINARY.csv"))

# --- 2. Maximum Intensity Matrix ---
cat("  - Generating Max Intensity Matrix...\n")
feature_matrix_max <- peaks_no_noise %>%
  group_by(SampleName, GroupName, Cluster) %>%
  summarise(Intensity = max(Intensity), .groups = 'drop') %>%
  pivot_wider(
    id_cols = c(SampleName, GroupName),
    names_from = Cluster,
    values_from = Intensity,
    values_fill = 0,
    names_prefix = "Cluster_"
  )
write_csv(feature_matrix_max, file.path(feature_matrix_dir, "feature_matrix_MAX_INTENSITY.csv"))

# --- 3. Median Intensity Matrix ---
cat("  - Generating Median Intensity Matrix...\n")
feature_matrix_median <- peaks_no_noise %>%
  group_by(SampleName, GroupName, Cluster) %>%
  summarise(Intensity = median(Intensity), .groups = 'drop') %>%
  pivot_wider(
    id_cols = c(SampleName, GroupName),
    names_from = Cluster,
    values_from = Intensity,
    values_fill = 0,
    names_prefix = "Cluster_"
  )
write_csv(feature_matrix_median, file.path(feature_matrix_dir, "feature_matrix_MEDIAN_INTENSITY.csv"))


# --- 4. ROI Integration Matrix ---
cat("  - Generating ROI Integration Matrix...\n")
# --- Step 4a: Define ROIs automatically from cluster boundaries ---
roi_definitions <- peaks_no_noise %>%
  group_by(Cluster) %>%
  summarise(
    RT_min = min(RetentionTime) - ROI_RT_PADDING,
    RT_max = max(RetentionTime) + ROI_RT_PADDING,
    DT_min = min(DriftTime) - ROI_DT_PADDING,
    DT_max = max(DriftTime) + ROI_DT_PADDING,
    .groups = 'drop'
  )
feature_rois <- setNames(split(roi_definitions, seq(nrow(roi_definitions))), paste0("Cluster_", roi_definitions$Cluster))
roi_plot_df <- bind_rows(feature_rois, .id = "FeatureName")

# --- Step 4b: Loop through samples to integrate signal and create plots ---
all_sample_features_roi <- list()
all_sample_files <- list.files(
  path = base_data_dir,
  recursive = TRUE,
  pattern = "_processed\\.csv$", # <-- More specific pattern
  full.names = TRUE
)

for (file_path in all_sample_files) {
  sample_name <- tools::file_path_sans_ext(basename(file_path))
  cat(paste0("    - Integrating ROIs for sample: ", sample_name, "\n"))
  
  data_matrix_df <- read_csv(file_path, show_col_types = FALSE)
  drift_times <- data_matrix_df[[1]]
  retention_times <- as.numeric(colnames(data_matrix_df)[-1])
  
  current_sample_values <- c()
  for (feature_name in names(feature_rois)) {
    roi_df <- feature_rois[[feature_name]]
    row_indices <- which(drift_times >= roi_df$DT_min & drift_times <= roi_df$DT_max)
    col_indices <- which(retention_times >= roi_df$RT_min & retention_times <= roi_df$RT_max)
    col_indices_in_df <- col_indices + 1
    
    if (length(row_indices) > 0 && length(col_indices_in_df) > 0) {
      integrated_value <- sum(as.matrix(data_matrix_df[row_indices, col_indices_in_df, drop = FALSE]), na.rm = TRUE)
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

feature_matrix_roi <- bind_rows(all_sample_features_roi, .id = "SampleName")
# Join with metadata to get GroupName
metadata_simple <- all_peaks %>% distinct(SampleName, GroupName)
feature_matrix_roi <- left_join(metadata_simple, feature_matrix_roi, by = "SampleName")
write_csv(feature_matrix_roi, file.path(feature_matrix_dir, "feature_matrix_ROI_INTEGRATION.csv"))

cat("✅  All four feature matrices generated and saved.\n")
cat("\n🎉 --- MASTER WORKFLOW COMPLETE! --- 🎉\n")






# =============================================================================
#
#   --- MASTER SCRIPT V14: BORUTA SELECTION & WILCOXON CHARACTERIZATION ---
#
# This script performs a complete and robust analysis pipeline:
#
#   1. FEATURE SELECTION: Uses the Boruta algorithm to identify a robust,
#      multivariate biomarker signature.
#   2. STATISTICAL CHARACTERIZATION: Performs a non-parametric Wilcoxon
#      rank-sum test on EACH of the Boruta-selected features to report
#      their individual p-values and significance.
#   3. MULTIVARIATE VALIDATION: Runs a 4-model gauntlet (LDA, RF, SVM, XGBoost)
#      on the complete Boruta-selected signature to assess its combined
#      predictive power with LOOCV.
#   4. REPORTING: Generates a final comparative performance summary, ROC
#      plots, and interaction plots for the final signature.
#
# =============================================================================

# --- 1. LOAD REQUIRED LIBRARIES ---
# install.packages(c("readr", "dplyr", "Boruta", "MASS", "caret", "ggplot2",
#                    "GGally", "vegan", "pROC", "e1071", "xgboost", "tibble"))

library(readr)
library(dplyr)
library(Boruta)
library(MASS)
library(caret)
library(ggplot2)
library(GGally)
library(vegan)
library(pROC)
library(e1071)
library(xgboost)
library(tibble)



# =============================================================================
# PART A: CONFIGURE AND LOAD DATA
# =============================================================================
cat("--- PART A: LOADING AND PREPARING DATA ---\n")
feature_matrix_path <-"C:/Users/deevanshi.walia/Desktop/Try 2.0/Onset of Symptoms/Final_Analysis_Output_Peaks_Clustering_9_0.05/Feature_Matrices/feature_matrix_MEDIAN_INTENSITY.csv"
data_full <- read_csv(feature_matrix_path, show_col_types = FALSE)
data_full$GroupName <- as.factor(data_full$GroupName)
if (any(is.na(data_full))) {
  warning("NA values detected. Removing rows with NAs.")
  data_full <- na.omit(data_full)
}
metadata <- data_full %>% dplyr:: select(SampleName, GroupName)
feature_data <- data_full %>% dplyr::select(starts_with("Cluster_"))
cat("✅  Data loaded successfully.\n")

# =============================================================================
# PART B: TWO-STAGE FEATURE SELECTION
# =============================================================================
cat("\n--- PART B: TWO-STAGE FEATURE SELECTION ---\n")

# --- STAGE 1: Broad Screen with Boruta ---
cat("\n--- B.1: Screening features with Boruta ---\n")
set.seed(456)
boruta_result <- Boruta(x = feature_data, y = metadata$GroupName, doTrace = 2, maxRuns = 500)
boruta_selected_features <- getSelectedAttributes(boruta_result, withTentative = FALSE)
cat("\n✅  Boruta identified the following candidate features:\n")
print(boruta_selected_features)

# --- STAGE 2: Statistical Confirmation with Wilcoxon Test ---
cat("\n--- B.2: Confirming features with Wilcoxon Test (p < 0.05) ---\n")
if (length(boruta_selected_features) > 0) {
  univariate_results_df <- lapply(boruta_selected_features, function(feature_name) {
    test_result <- wilcox.test(feature_data[[feature_name]] ~ metadata$GroupName)
    data.frame(Feature = feature_name, p.value = test_result$p.value)
  }) %>%
    bind_rows()
  
  SIGNIFICANCE_THRESHOLD <- 0.05
  selected_features <- univariate_results_df %>%
    filter(p.value < SIGNIFICANCE_THRESHOLD) %>%
    pull(Feature)
  
} else {
  selected_features <- c()
}

if (length(selected_features) == 0) {
  stop("The two-stage selection process resulted in zero features.")
}

cat("\n✅  The final, high-confidence signature (Boruta-confirmed AND p < 0.05) is:\n")
print(selected_features)
feature_data_best <- data_full %>% dplyr::select(all_of(selected_features))

# =============================================================================
# PART C: MODEL GAUNTLET - LEAVE-ONE-OUT CROSS-VALIDATION
# =============================================================================
cat("\n--- PART C: RUNNING MODEL GAUNTLET ON FINAL SIGNATURE ---\n")
num_samples <- nrow(feature_data_best)
y_response <- as.factor(metadata$GroupName)
y_response_numeric <- as.numeric(y_response) - 1
current_group_levels <- levels(y_response)
positive_class_name <- current_group_levels[1]
loocv_results <- data.frame(SampleName = metadata$SampleName, ActualGroup = y_response,
                            lda_scores=numeric(num_samples), rf_scores=numeric(num_samples),
                            svm_scores=numeric(num_samples), xgb_scores=numeric(num_samples))

cat("    - Running LOOCV for all 4 models...\n")
for (i in 1:num_samples) {
  cat(paste0("\r  - [", i, "/", num_samples, "]")); flush.console()
  x_train <- as.matrix(feature_data_best[-i, , drop=FALSE]); y_train <- y_response[-i]
  y_train_numeric <- y_response_numeric[-i]; x_test <- as.matrix(feature_data_best[i, , drop = FALSE])
  
  lda_model <- lda(x = x_train, grouping = y_train)
  loocv_results$lda_scores[i] <- predict(lda_model, newdata = data.frame(x_test))$posterior[, positive_class_name]
  
  min_class_size <- min(table(y_train))
  rf_model <- randomForest(x = x_train, y = y_train, ntree = 500, sampsize = c(min_class_size, min_class_size))
  loocv_results$rf_scores[i] <- predict(rf_model, newdata = x_test, type = "prob")[, positive_class_name]
  
  svm_model <- svm(x = x_train, y = y_train, probability = TRUE, kernel="radial")
  loocv_results$svm_scores[i] <- attr(predict(svm_model, newdata = data.frame(x_test), probability = TRUE), "probabilities")[, positive_class_name]
  
  dtrain <- xgb.DMatrix(data = x_train, label = y_train_numeric)
  dtest <- xgb.DMatrix(data = x_test)
  scale_pos_weight <- sum(y_train_numeric == 0) / sum(y_train_numeric == 1)
  params <- list(objective="binary:logistic", eval_metric="logloss", scale_pos_weight=scale_pos_weight)
  xgb_model <- xgboost(params = params, data = dtrain, nrounds = 50, verbose = 0)
  loocv_results$xgb_scores[i] <- predict(xgb_model, dtest)
}
cat("\n✅  LOOCV complete for all models.\n")

# =============================================================================
# PART D: PERFORMANCE COMPARISON AND VISUALIZATION
# =============================================================================
cat("\n--- PART D: COMPARING MODEL PERFORMANCES ---\n")
roc_lda <- roc(response = loocv_results$ActualGroup, predictor = loocv_results$lda_scores, levels=current_group_levels)
roc_rf  <- roc(response = loocv_results$ActualGroup, predictor = loocv_results$rf_scores, levels=current_group_levels)
roc_svm <- roc(response = loocv_results$ActualGroup, predictor = loocv_results$svm_scores, levels=current_group_levels)
roc_xgb <- roc(response = loocv_results$ActualGroup, predictor = loocv_results$xgb_scores, levels=current_group_levels)

model_performance_list <- list()
models_to_test <- list(LDA = roc_lda, `Random Forest` = roc_rf, SVM = roc_svm, XGBoost = roc_xgb)
cm_list <- list() 
for (model_name in names(models_to_test)) {
  current_roc <- models_to_test[[model_name]]
  best_coords <- coords(current_roc, "best", ret = c("threshold", "youden"), best.method="youden")
  best_thresh <- best_coords$threshold; youdens_j <- best_coords$youden
  predicted_classes <- ifelse(current_roc$predictor >= best_thresh, current_group_levels[2], current_group_levels[1])
  cm_stats <- confusionMatrix(data = factor(predicted_classes, levels = current_group_levels), reference = current_roc$response, positive = positive_class_name)
  cm_list[[model_name]] <- cm_stats
  model_performance_list[[model_name]] <- c(AUC=auc(current_roc), cm_stats$overall["Accuracy"], cm_stats$byClass[c("Sensitivity", "Specificity", "Balanced Accuracy")], "Youden_J"=youdens_j)
}
performance_summary <- as.data.frame(do.call(rbind, model_performance_list)) %>%
  tibble::rownames_to_column("Model") %>%
  mutate(across(where(is.numeric), ~ round(., 4))) %>%
  dplyr::select(Model, AUC, Accuracy, `Balanced Accuracy`, Youden_J, Sensitivity, Specificity)

cat("\n--- Performance Summary Table (at Optimal Youden's J Threshold) ---\n")
print(performance_summary)

cat("\n🎨  Generating combined ROC Curve Plot...\n")
roc_list_for_plot <- list(LDA = roc_lda, RF = roc_rf, SVM = roc_svm, XGBoost = roc_xgb)
roc_plot_combined <- ggroc(roc_list_for_plot, size = 1) +
  annotate("segment", x = 1, xend = 0, y = 0, yend = 1, color="grey", linetype="dashed") +
  labs(title = "Comparative ROC Curves for All Models", color = "Model") +
  theme_minimal(base_size = 14) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
print(roc_plot_combined)

best_model_name <- performance_summary$Model[which.max(performance_summary$`Balanced Accuracy`)]
cat("\n--- Detailed Confusion Matrix for the Best Performing Model:", best_model_name, "---\n")
print(cm_list[[best_model_name]])

# =============================================================================
# PART B: PAIRS PLOT OF ALL SELECTED FEATURES
# =============================================================================

cat("\n--- PART B: Generating pairs plot for all selected features ---\n")

# A pairs plot is most useful if you have between 2 and ~6 features.
if (length(selected_features) >= 2) {
  
  # 1. Select the data for the plot
  # This includes the GroupName and all the selected cluster intensities
  data_for_pairs_plot <- cbind(
    GroupName = metadata$GroupName,
    feature_data[, selected_features]
  )
  
  # 2. Create the pairs plot using the powerful ggpairs function
  pairs_plot <- ggpairs(
    data_for_pairs_plot,
    # The 'mapping' tells it to color the points by GroupName
    mapping = aes(color = GroupName, alpha = 0.7),
    # Tell it which columns to use for the grid
    columns = 2:ncol(data_for_pairs_plot),
    # Define what to show in each part of the grid
    lower = list(continuous = wrap("points", size = 2)), # Scatter plots in the lower triangle
    diag = list(continuous = wrap("densityDiag", alpha = 0.5)), # Density plots on the diagonal
    upper = list(continuous = "blank") # Keep the upper triangle clean
  ) +
    labs(title = paste("Pairwise Interactions of", length(selected_features), "Selected Features")) +
    theme(plot.title = element_text(face="bold", size=16, hjust=0.5))
  
  print(pairs_plot)
  cat("✅  Pairs plot generated.\n")
  
} else {
  cat("   - Skipping pairs plot: Fewer than 2 features were selected.\n")
  # If only one feature is selected, a boxplot is the most appropriate visualization
  boxplot_data <- cbind(metadata, feature_data[, selected_features, drop=FALSE])
  single_feature_name <- selected_features[1]
  
  single_feature_boxplot <- ggplot(boxplot_data, aes(x=GroupName, y=.data[[single_feature_name]], fill=GroupName)) +
    geom_boxplot(alpha=0.7) +
    geom_jitter(width=0.2, alpha=0.5) +
    labs(
      title = paste("Distribution of the Top Selected Feature:", single_feature_name),
      x = "Group",
      y = "Integrated Intensity"
    ) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "none")
  
  print(single_feature_boxplot)
}

cat("\n🎉 --- Visualization Complete! --- 🎉\n")

#   --- SCRIPT: BOXPLOT VISUALIZATION OF BORUTA-SELECTED FEATURES ---
#                  WITH STATISTICAL ANNOTATION
#
# Description:
# This script takes the final list of features identified by the Boruta
# algorithm. For each feature, it performs a Wilcoxon rank-sum test and
# generates a boxplot comparing the two groups, with the resulting
# p-value clearly annotated on the plot.
#
# =============================================================================

# --- 1. Load Required Libraries ---
# install.packages(c("dplyr", "readr", "ggplot2", "tidyr", "ggpubr"))
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(ggpubr) # For the stat_compare_means() function

# --- PREREQUISITE: Assumes you have the following objects from your master script ---
#   - `data_full`: Your complete data frame with 'GroupName' and all features.
#   - `selected_features`: The character vector of names for your best features from Boruta.
#
#   Example placeholder data (uncomment to run this script standalone):
#   selected_features <- c("Cluster_5", "Cluster_15", "Cluster_17")
#   data_full <- data.frame(
#     GroupName = factor(rep(c("Group_A", "Group_B"), each = 50)),
#     Cluster_5 = rnorm(100, mean = ifelse(rep(c("A","B"), each=50) == "A", 10, 12)),
#     Cluster_15 = rnorm(100, mean = ifelse(rep(c("A","B"), each=50) == "A", 20, 25)),
#     Cluster_17 = rnorm(100, mean = ifelse(rep(c("A","B"), each=50) == "A", 5, 5.5))
#   )


# =============================================================================
# PART A: PREPARE DATA FOR PLOTTING
# =============================================================================
cat("--- PART A: PREPARING DATA FOR VISUALIZATION ---\n")

# Check if we have features to plot
if (length(selected_features) == 0) {
  stop("The 'selected_features' vector is empty. Cannot create plots.")
}

# Reshape the data from "wide" format to "long" format, keeping only the selected features.
# This is the ideal format for creating faceted plots in ggplot2.
plot_data_long <- data_full %>%
  dplyr::select(GroupName, all_of(selected_features)) %>%
  pivot_longer(
    cols = -GroupName,           # Pivot all columns except GroupName
    names_to = "Feature",        # The new column for the feature names
    values_to = "Intensity"      # The new column for the intensity values
  ) %>%
  # Ensure the order of features in the plot is the same as your 'selected_features' vector
  mutate(Feature = factor(Feature, levels = selected_features))

cat("✅  Data prepared for plotting.\n")


# =============================================================================
# PART B: GENERATE THE ANNOTATED BOXPLOT PANEL
# =============================================================================
cat("\n--- PART B: GENERATING BOXPLOTS WITH P-VALUE ANNOTATIONS ---\n")

# Create the plot
boxplots_with_pvals <- ggplot(
  data = plot_data_long,
  aes(x = GroupName, y = Intensity, fill = GroupName)
) +
  # Add the boxplots
  geom_boxplot(
    alpha = 0.7,
    outlier.shape = NA # Hide default outliers so jitter can show them
  ) +
  
  # Overlay the individual data points to show the distribution
  geom_jitter(
    width = 0.2, # How much to jitter horizontally
    alpha = 0.4,
    size = 1.5
  ) +
  
  # --- THIS IS THE KEY FUNCTION FOR ADDING P-VALUES ---
  # `stat_compare_means` from the 'ggpubr' package automatically runs the
  # specified test and adds the p-value to the plot.
  stat_compare_means(
    method = "wilcox.test", # Specify the Wilcoxon test
    label = "p.format",     # Display the p-value in a nice format (e.g., p = 0.02)
    label.x.npc = 0.5,      # Position the label in the horizontal center
    label.y.npc = 0.9,      # Position the label 90% of the way up the plot
    size = 4
  ) +
  
  # Create a separate panel for each feature
  # `scales = "free_y"` allows each feature to have its own y-axis range
  facet_wrap(~ Feature, scales = "free_y", ncol = 3) + # Arrange in 3 columns
  
  # Apply custom colors and labels
  scale_fill_manual(
    name = "Group",
    values = c("#0073C2FF", "#EFC000FF") # Example colors: Blue and Yellow
  ) +
  
  labs(
    title = "Intensity Distribution of Boruta-Selected Features",
    subtitle = "P-values from Wilcoxon Rank-Sum Test",
    x = "Group",
    y = "Integrated Intensity"
  ) +
  
  # Use a clean theme
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    strip.text = element_text(face = "bold"), # Make the panel titles bold
    legend.position = "none", # Hide legend as colors are clear from x-axis
    axis.text.x = element_text(angle = 45, hjust = 1) # Angle x-axis labels
  )

# Print the plot to the RStudio "Plots" pane
print(boxplots_with_pvals)

# Save the plot to a file
output_dir <- dirname(feature_matrix_path) # Save near your input data
plot_output_path <- file.path(output_dir, "Boxplots_of_Selected_Features_with_PValues.png")
ggsave(
  plot_output_path,
  plot = boxplots_with_pvals,
  width = 12,
  height = 5 * ceiling(length(selected_features) / 3), # Auto-adjust height
  dpi = 300,
  bg = "white"
)

cat("\n✅  Annotated boxplot panel saved to:", plot_output_path, "\n")
cat("\n🎉 --- Visualization Complete --- 🎉\n")


