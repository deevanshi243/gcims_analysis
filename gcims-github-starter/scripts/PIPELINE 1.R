


#=============================================================================
#
#   --- MASTER ANALYSIS PIPELINE ---
#
# ===========================================================================

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
# PART A: LOAD AND PREPARE DATA
# =============================================================================

cat("--- PART A: LOADING AND PREPARING DATA ---\n")
feature_matrix_path <-  "C:/Users/deevanshi.walia/Desktop/Try 2.0/Onset of Symptoms/Final_Analysis_Output_Peaks_Clustering_9_0.05/Feature_Matrices/feature_matrix_ROI_INTEGRATION.csv"

data_full <- read_csv(feature_matrix_path, show_col_types = FALSE)
metadata <- data_full %>% dplyr::select(SampleName, GroupName)
feature_data <- data_full %>% dplyr::select(starts_with("Cluster_"))

if (any(is.na(data_full))) {
  warning("NA values detected. Removing rows with NAs.")
  data_full <- na.omit(data_full)
  metadata <- data_full %>% dplyr::select(SampleName, GroupName)
  feature_data <- data_full %>% dplyr::select(starts_with("Cluster_"))
}
cat("✅  Data loaded successfully with", nrow(feature_data), "samples and", ncol(feature_data), "features.\n")


# =============================================================================
# PART B: ROBUST FEATURE SELECTION WITH BORUTA
# =============================================================================

cat("\n--- PART B: PERFORMING ROBUST FEATURE SELECTION WITH BORUTA ---\n")
set.seed(321)
boruta_result <- Boruta(
  x = feature_data,
  y = as.factor(metadata$GroupName),
  doTrace = 2,
  maxRuns = 700
)

cat("\n✅ Boruta analysis complete.\n")
print(boruta_result)
save_result_path<-file.path(output_base_dir, "boruta_results.png")
png(save_result_path, width = 8, height = 6, units = "in", res = 300)
par(mar = c(8, 5, 4, 2))

plot(boruta_result,
     las = 2,
     cex.axis = 0.7,
     xlab = "",        
     ylab = ""        
)

# 4. Add the axis titles back manually with mtext()
#    - side=1 is bottom, side=2 is left.
#    - 'line' controls the distance from the plot. Adjust as needed.
#    - 'cex' controls the font size of the title.
mtext("Attributes", side = 1, line = 5, cex = 1)
mtext("Importance (Z-score)", side = 2, line = 2.5, cex = 1)


# 5. Close the device to save the file
dev.off()

selected_features <- getSelectedAttributes(boruta_result, withTentative = FALSE)
if (length(selected_features) < 1) stop("Boruta did not confirm any features.")

cat("\n✅  Boruta identified the definitive feature signature:\n")
print(selected_features)
feature_data_best <- feature_data[, selected_features, drop = FALSE]


# =============================================================================
# PART C: STATISTICAL SIGNIFICANCE TESTING OF THE SELECTED SIGNATURE
# =============================================================================

cat("\n--- PART C: TESTING SIGNIFICANCE WITH NON-PARAMETRIC PERMANOVA ---\n")
set.seed(987)
permanova_result <- adonis2(feature_data_best ~ GroupName, data = metadata, permutations = 9999, method="euclidean")
final_p_value <- permanova_result$`Pr(>F)`[1]

cat("✅  Significance test complete.\n")
cat("    - PERMANOVA P-value for the ", length(selected_features), "-feature signature: ", format.pval(final_p_value, digits = 4), "\n", sep="")
if (final_p_value < 0.05) {
  cat("    - CONCLUSION: The separation between groups is STATISTICALLY SIGNIFICANT.\n")
} else {
  cat("    - CONCLUSION: The separation is NOT statistically significant.\n")
}
# --- D2: Univariate Test (Wilcoxon Rank-Sum Test) on EACH selected feature ---
#cat("\n--- PART D.2: TESTING UNIVARIATE SIGNIFICANCE (WILCOXON TEST) ---\n")

# Create an empty list to store the results
#univariate_results_list <- list()

# Loop through ONLY the features that Boruta selected
#for (feature_name in selected_features) {
  
  # Create a formula for the test, e.g., Cluster_15 ~ GroupName
  #test_formula <- as.formula(paste("`", feature_name, "` ~ GroupName", sep=""))
  
  
  # Perform the Wilcoxon Rank-Sum test for two independent groups
 # test_result <- wilcox.test(
   # formula = test_formula,
   # data = cbind(feature_data_best, metadata)
 # )
  
  # Store the results
 # univariate_results_list[[feature_name]] <- data.frame(
  #  Feature = feature_name,
  #  W_statistic = test_result$statistic,
  #  p.value = test_result$p.value
#  )
#}

# Combine the results into a single, clean data frame
#univariate_results_df <- bind_rows(univariate_results_list)

# Apply a Benjamini-Hochberg correction to the p-values
# This is still best practice, even for a small number of tests.
#univariate_results_df$p.adj <- p.adjust(univariate_results_df$p.value, method = "BH")

# Sort the results by the adjusted p-value
#univariate_results_df <- univariate_results_df %>% arrange(p.adj)

#cat("--- Wilcoxon Rank-Sum Test Results for Boruta-Selected Features ---\n")
#print(univariate_results_df)
#######################################cleaner features
#high_confidence_features <- univariate_results_df %>%
 # filter(p.adj < 0.05) %>%
 # pull(Feature)

#cat("--- High-Confidence Feature Signature ---\n")
#print(high_confidence_features)


# =============================================================================
#
#   --- PART D: BALANCED MODEL GAUNTLET WITH LEAVE-ONE-OUT CV (DEFINITIVE) ---
#
# Description:
# This part performs a rigorous, comparative validation of four distinct
# classification models. It uses Leave-One-Out Cross-Validation (LOOCV) to
# get an unbiased performance estimate. Crucially, it applies class
# balancing techniques to Random Forest, SVM, and XGBoost to ensure a fair
# comparison on the imbalanced dataset.
#
# =============================================================================

cat("\n--- PART D: RUNNING BALANCED MODEL GAUNTLET WITH LOOCV ---\n")

# --- Initialization ---
# Assumes 'data_full' and 'selected_features' exist from previous parts.
metadata <- data_full %>% dplyr::select(SampleName, GroupName)
feature_data_best <- data_full %>% dplyr::select(all_of(selected_features))

num_samples <- nrow(feature_data_best)
y_response <- as.factor(metadata$GroupName)
y_response_numeric <- as.numeric(y_response) - 1
current_group_levels <- levels(y_response)
positive_class_name <- current_group_levels[1]

cat("    - Class of interest set to:", positive_class_name, "\n")
cat("    - Performing LOOCV for", num_samples, "samples...\n")

# Data frame to store the raw predictive scores from all models
loocv_results <- data.frame(
  SampleName = metadata$SampleName,
  ActualGroup = y_response,
  lda_scores = numeric(num_samples),
  rf_scores = numeric(num_samples),
  svm_scores = numeric(num_samples),
  xgb_scores = numeric(num_samples)
)

# --- The Grand LOOCV Loop ---
for (i in 1:num_samples) {
  # Progress bar
  cat(paste0("\r  - [", i, "/", num_samples, "]")); flush.console()
  
  # --- 1. Split data for this fold ---
  x_train <- as.matrix(feature_data_best[-i, , drop = FALSE])
  y_train <- y_response[-i]
  y_train_numeric <- y_response_numeric[-i]
  x_test <- as.matrix(feature_data_best[i, , drop = FALSE])
  
  # --- 2. Train and Predict with LDA (Unbalanced Baseline) ---
  # Standard LDA has no direct balancing parameter.
  lda_model <- lda(x = x_train, grouping = y_train)
  lda_pred <- predict(lda_model, newdata = data.frame(x_test))
  loocv_results$lda_scores[i] <- lda_pred$posterior[, positive_class_name]
  
  # --- 3. Train and Predict with Random Forest (BALANCED) ---
  # Use 'sampsize' for balancing
  min_class_size <- min(table(y_train))
  rf_model <- randomForest(
    x = x_train, y = y_train, ntree = 500,
    sampsize = c(min_class_size, min_class_size)
  )
  rf_pred_prob <- predict(rf_model, newdata = x_test, type = "prob")
  loocv_results$rf_scores[i] <- predict(rf_model, newdata = x_test, type = "prob")[, 2]
  
  # --- 4. Train and Predict with SVM (BALANCED) ---
  # Use 'class.weights' for balancing
  class_weights_svm <- 1 / table(y_train)
  svm_model <- svm(
    x = x_train, y = y_train, probability = TRUE, kernel = "radial",
    class.weights = class_weights_svm
  )
  svm_pred <- predict(svm_model, newdata = data.frame(x_test), probability = TRUE)
  loocv_results$svm_scores[i] <- attr(svm_pred, "probabilities")[, positive_class_name]
  
  # --- 5. Train and Predict with XGBoost (BALANCED) ---
  # Use 'scale_pos_weight' for balancing
  dtrain <- xgb.DMatrix(data = x_train, label = y_train_numeric)
  dtest <- xgb.DMatrix(data = x_test)
  scale_pos_weight <- sum(y_train_numeric == 0) / sum(y_train_numeric == 1)
  params <- list(objective = "binary:logistic", eval_metric = "logloss", scale_pos_weight = scale_pos_weight)
  xgb_model <- xgboost(params = params, data = dtrain, nrounds = 50, verbose = 0)
  loocv_results$xgb_scores[i] <- predict(xgb_model, dtest)
}
cat("\n✅  LOOCV complete for all models.\n")

# =============================================================================
# PART E: PERFORMANCE COMPARISON WITH YOUDEN'S J INDEX
# =============================================================================
cat("\n--- PART E: COMPARING MODEL PERFORMANCES ---\n")

# --- Calculate ROC curves ---
roc_lda <- roc(response = loocv_results$ActualGroup, predictor = loocv_results$lda_scores, levels=current_group_levels)
roc_rf  <- roc(response = loocv_results$ActualGroup, predictor = loocv_results$rf_scores, levels=current_group_levels)
roc_svm <- roc(response = loocv_results$ActualGroup, predictor = loocv_results$svm_scores, levels=current_group_levels)
roc_xgb <- roc(response = loocv_results$ActualGroup, predictor = loocv_results$xgb_scores, levels=current_group_levels)

# --- Calculate Performance Metrics for Each Model ---
model_performance_list <- list()
models_to_test <- list(LDA = roc_lda, `Random Forest` = roc_rf, SVM = roc_svm, XGBoost = roc_xgb)
cm_list <- list() 

for (model_name in names(models_to_test)) {
  current_roc <- models_to_test[[model_name]]
  
  # --- Find the optimal threshold using YOUDEN'S J INDEX ---
  best_coords <- coords(current_roc, "best", ret = c("threshold", "youden"), best.method="youden")
  best_thresh <- best_coords$threshold
  youdens_j <- best_coords$youden
  
  predicted_classes <- ifelse(current_roc$predictor >= best_thresh, current_group_levels[2], current_group_levels[1])
  cm_stats <- confusionMatrix(data = factor(predicted_classes, levels = current_group_levels), reference = current_roc$response, positive = positive_class_name)
  cm_list[[model_name]] <- cm_stats
  
  # --- Add Youden's J to the summary ---
  model_performance_list[[model_name]] <- c(
    AUC = auc(current_roc),
    cm_stats$overall["Accuracy"],
    cm_stats$byClass[c("Sensitivity", "Specificity", "Balanced Accuracy", "Pos Pred Value", "Neg Pred Value")],
    "Youden_J" = youdens_j
  )
}

# --- Create a clean summary table ---
performance_summary <- as.data.frame(do.call(rbind, model_performance_list)) %>%
  tibble::rownames_to_column("Model") %>%
  mutate(across(where(is.numeric), ~ round(., 4))) %>%
  # Reorder for clarity, putting Youden's J next to Balanced Accuracy
 dplyr:: select(Model, AUC, Accuracy, `Balanced Accuracy`, Youden_J, Sensitivity, Specificity,`Neg Pred Value`,`Pos Pred Value`)

cat("\n--- Performance Summary Table (at Optimal Youden's J Threshold) ---\n")
print(performance_summary)
write.table(performance_summary, 
            file = "C:/Users/deevanshi.walia/Desktop/performance_summary.csv", 
            sep = ";", 
            row.names = FALSE)


# --- Create the combined ROC Curve Plot ---
cat("\n🎨  Generating combined ROC Curve Plot...\n")
roc_list_for_plot <- list(LDA = roc_lda, RF = roc_rf, SVM = roc_svm, XGBoost = roc_xgb)
roc_plot_combined <- ggroc(roc_list_for_plot, size = 1) +
  annotate("segment", x = 1, xend = 0, y = 0, yend = 1, color="grey", linetype="dashed") +
  labs(title = "Comparative ROC Curves for All Models", color = "Model") +
  theme_minimal(base_size = 14) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
print(roc_plot_combined)


# --- Print the best model's detailed confusion matrix ---
best_model_name <- performance_summary$Model[which.max(performance_summary$`Balanced Accuracy`)]
cat("\n--- Detailed Confusion Matrix for the Best Performing Model:", best_model_name, "---\n")
print(cm_list[[best_model_name]])

# Extract the overall statistics
cm_object <- cm_list[[best_model_name]]
overall_stats <- as.data.frame(t(cm_object$overall))

# Extract class-specific metrics (like sensitivity, specificity)
by_class_stats <- as.data.frame(t(cm_object$byClass))

# Combine both if you want
full_stats <- cbind(overall_stats, by_class_stats)
###################################################################################FIX
write.table(full_stats,
            file = "C:/Users/deevanshi.walia/Desktop/best_model_metrics.csv",
            sep = ";",
            row.names = FALSE)

cat("\n🎉 --- Master Analysis Workflow Complete! --- 🎉\n")
# =============================================================================
# PART B: PAIRS PLOT OF ALL SELECTED FEATURES
# =============================================================================

cat("\n--- PART B: Generating pairs plot for all selected features ---\n")
num_features <- length(selected_features)

# Initialize a variable to hold our final plot and its filename
final_plot <- NULL
plot_filename <- ""

if (num_features >= 2) {
  cat("   - Found", num_features, "features. Generating ggpairs plot...\n")
  
  data_for_pairs_plot <- cbind(
    GroupName = metadata$GroupName,
    feature_data[, selected_features]
  )
  
  # Assign the created plot to our 'final_plot' variable
  final_plot <- GGally::ggpairs(
    data_for_pairs_plot,
    mapping = aes(color = GroupName, alpha = 0.7),
    columns = 2:ncol(data_for_pairs_plot),
    lower = list(continuous = wrap("points", size = 2)),
    diag = list(continuous = wrap("densityDiag", alpha = 0.5)),
    upper = list(continuous = "blank")
  ) +
    labs(title = paste("Pairwise Interactions of", num_features, "Selected Features")) +
    theme(plot.title = element_text(face="bold", size=16, hjust=0.5))
  
  # Set the desired filename
  plot_filename <- "selected_features_pairs_plot.png"
  
} else if (num_features == 1) {
  cat("   - Found 1 feature. Generating a boxplot...\n")
  
  boxplot_data <- cbind(metadata, feature_data[, selected_features, drop=FALSE])
  single_feature_name <- selected_features[1]
  
  # Assign the created plot to our 'final_plot' variable
  final_plot <- ggplot(boxplot_data, aes(x=GroupName, y=.data[[single_feature_name]], fill=GroupName)) +
    geom_boxplot(alpha=0.7) +
    geom_jitter(width=0.2, alpha=0.5) +
    labs(
      title = paste("Distribution of the Top Selected Feature:", single_feature_name),
      x = "Group",
      y = "Integrated Intensity"
    ) +
    theme_minimal(base_size = 16, face= "bold") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size=14), legend.position = "none")
  
  # Set the desired filename
  plot_filename <- "single_feature_boxplot.png"
  
} else {
  cat("   - Skipping plot: 0 features were selected.\n")
}

# --- Step 2: Print and Save the final plot (if one was created) ---

if (!is.null(final_plot)) {
  
  # Print the plot to the RStudio viewer
  print(final_plot)
  cat("✅ Plot generated.\n")
  
  # Save the plot to a file using the filename defined earlier
  # Using ggsave() is the modern and easier way for ggplot objects.
  ggsave(
    plot_filename,
    plot = final_plot,
    width = 12,
    height = 8,
    dpi = 600,
    units = "in"
  )
  dev.off()
}
########################################################
######################################################boxplots with inten and p values wilcoxins
# =============================================================================
#
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
    values = c("darkgreen", "#EFC000FF") # Example colors: Blue and Yellow
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
    strip.text = element_text(face = "bold",  size = 14), # Make the panel titles bold
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
  width = 15,
  height = 5 * ceiling(length(selected_features) / 4), # Auto-adjust height
  dpi = 300,
  bg = "white"
)
dev.off()
cat("\n✅  Annotated boxplot panel saved to:", plot_output_path, "\n")
cat("\n🎉 --- Visualization Complete --- 🎉\n")

