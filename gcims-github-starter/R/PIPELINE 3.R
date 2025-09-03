#   --- MASTER SCRIPT V15.1: UNIVARIATE SELECTION & 4-MODEL GAUNTLET ---
#
# This definitive script uses a classical statistical approach for feature
# selection and then validates the resulting signature with a comprehensive
# 4-model gauntlet and final visualizations.
#
# =============================================================================

# --- 1. LOAD REQUIRED LIBRARIES ---
library(dplyr); library(readr); library(caret); library(tidyr); library(MASS)
library(pROC); library(ggplot2); library(randomForest); library(e1071)
library(xgboost); library(GGally); library(tibble); library(vegan)

# =============================================================================
# PART A: CONFIGURE YOUR SCRIPT
# =============================================================================
cat("--- PART A: CONFIGURATION ---\n")
feature_matrix_path <- "C:/Users/deevanshi.walia/Desktop/Try 2.0/Grouping acc to Jessy/Final_Analysis_Output_Peaks_Clustering_9_0.05/Feature_Matrices/feature_matrix_ROI_INTEGRATION.csv"
SIGNIFICANCE_THRESHOLD <- 0.05

# =============================================================================
# PART B: UNIVARIATE SCREENING (WILCOXON ON ALL FEATURES)
# =============================================================================
cat("\n--- PART B: UNIVARIATE SCREENING TO IDENTIFY CANDIDATE BIOMARKERS ---\n")
data_full <- read_csv(feature_matrix_path, show_col_types = FALSE)
data_full$GroupName <- as.factor(data_full$GroupName)
if (any(is.na(data_full))) {
  warning("NA values detected. Removing rows with NAs.")
  data_full <- na.omit(data_full)
}
metadata <- data_full %>% dplyr:: select(SampleName, GroupName)
feature_data <- data_full %>% dplyr::select(starts_with("Cluster_"))

stats_results_all <- lapply(names(feature_data), function(cluster_name) {
  test_result <- wilcox.test(feature_data[[cluster_name]] ~ metadata$GroupName)
  return(data.frame(Feature = cluster_name, p.value = test_result$p.value))
})
stats_df_all <- bind_rows(stats_results_all)
stats_df_all$p.adj <- p.adjust(stats_df_all$p.value, method = "BH")
stats_df_all <- stats_df_all %>% arrange(p.adj)

cat("--- Top Candidate Features from Univariate Screening ---\n")
print(head(stats_df_all))

selected_features <- stats_df_all %>%
  filter(p.adj < SIGNIFICANCE_THRESHOLD) %>%
  pull(Feature)

if (length(selected_features) == 0) {
  stop("No features were found to be statistically significant after BH correction. Cannot proceed.")
}
cat("\n--- Final Biomarker Signature (p.adj <", SIGNIFICANCE_THRESHOLD, ") ---\n")
print(selected_features)

# =============================================================================
# PART C: MULTIVARIATE CONFIRMATION (PERMANOVA ON SELECTED FEATURES)
# =============================================================================
cat("\n--- PART C: TESTING MULTIVARIATE SIGNIFICANCE OF THE SIGNATURE ---\n")
feature_data_best <- feature_data[, selected_features, drop = FALSE] # This is our final feature set
set.seed(123)
permanova_result <- adonis2(feature_data_best ~ GroupName, data = metadata, permutations = 9999, method = "euclidean")
final_p_value <- permanova_result$`Pr(>F)`[1]
cat("✅  PERMANOVA P-value for the ", length(selected_features), "-feature signature: ", format.pval(final_p_value, digits = 4), "\n", sep="")
if (final_p_value < 0.05) cat("    - CONCLUSION: The overall signature is STATISTICALLY SIGNIFICANT.\n")

# =============================================================================
# PART D: MODEL GAUNTLET - LEAVE-ONE-OUT CROSS-VALIDATION (CORRECTED)
# =============================================================================
cat("\n--- PART D: RUNNING 4-MODEL GAUNTLET WITH LOOCV ---\n")

# --- Initialization ---
# ❗ KEY FIX: Use `feature_data_best` which was defined in the previous step.
num_samples <- nrow(feature_data_best)
y_response <- as.factor(metadata$GroupName)
y_response_numeric <- as.numeric(y_response) - 1
current_group_levels <- levels(y_response)
positive_class_name <- current_group_levels[1]
loocv_results <- data.frame(SampleName=metadata$SampleName, ActualGroup=y_response,
                            lda_scores=numeric(num_samples), rf_scores=numeric(num_samples),
                            svm_scores=numeric(num_samples), xgb_scores=numeric(num_samples))

cat("    - Running LOOCV for LDA, Random Forest, SVM, and XGBoost...\n")
for (i in 1:num_samples) {
  cat(paste0("\r  - [", i, "/", num_samples, "]")); flush.console()
  
  # ❗ KEY FIX: Use `feature_data_best` in all these lines
  x_train <- as.matrix(feature_data_best[-i, , drop=FALSE])
  y_train <- y_response[-i]
  y_train_numeric <- y_response_numeric[-i]
  x_test <- as.matrix(feature_data_best[i, , drop = FALSE])
  
  # 1. LDA
  lda_model <- lda(x = x_train, grouping = y_train)
  lda_pred <- predict(lda_model, newdata = data.frame(x_test))
  loocv_results$lda_scores[i] <- lda_pred$posterior[, positive_class_name]
  
  # 2. Random Forest (balanced)
  min_class_size <- min(table(y_train))
  rf_model <- randomForest(x = x_train, y = y_train, ntree = 500, sampsize = c(min_class_size, min_class_size))
  rf_pred_prob <- predict(rf_model, newdata = x_test, type = "prob")
  loocv_results$rf_scores[i] <- rf_pred_prob[, positive_class_name]
  
  # 3. SVM
  svm_model <- svm(x = x_train, y = y_train, probability = TRUE, kernel="radial")
  svm_pred <- predict(svm_model, newdata = data.frame(x_test), probability = TRUE)
  loocv_results$svm_scores[i] <- attr(svm_pred, "probabilities")[, positive_class_name]
  
  # 4. XGBoost (balanced)
  dtrain <- xgb.DMatrix(data = x_train, label = y_train_numeric)
  dtest <- xgb.DMatrix(data = x_test)
  scale_pos_weight <- sum(y_train_numeric == 0) / sum(y_train_numeric == 1)
  params <- list(objective="binary:logistic", eval_metric="logloss", scale_pos_weight=scale_pos_weight)
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
  dplyr::select(Model, AUC, Accuracy, `Balanced Accuracy`, Youden_J, Sensitivity, Specificity,`Neg Pred Value`,`Pos Pred Value`)

cat("\n--- Performance Summary Table (at Optimal Youden's J Threshold) ---\n")
print(performance_summary)


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

cat("\n🎉 --- Master Analysis Workflow Complete! --- 🎉\n")
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

#######
# =============================================================================
# --- 3. Intensity Boxplots of All Selected Features ---
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
    title = "Intensity Distribution of Selected Features",
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


