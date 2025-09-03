# =============================================================================
#
#   --- GRAND MASTER SCRIPT V3: RFE-CV, 4-MODEL GAUNTLET, & FINAL SUMMARY ---
#
# This script performs a complete, publication-quality analysis pipeline from
# feature selection to final model validation and visualization.
#
# =============================================================================

# --- 1. LOAD REQUIRED LIBRARIES ---
# install.packages(c("dplyr", "readr", "caret", "tidyr", "MASS", "pROC",
#                    "ggplot2", "randomForest", "e1071", "xgboost", "GGally", "tibble"))
library(dplyr)
library(readr)
library(caret)
library(tidyr)
library(MASS)
library(pROC)
library(ggplot2)
library(randomForest)
library(e1071)
library(xgboost)
library(GGally)
library(tibble)

# =============================================================================
# PART A: CONFIGURE YOUR SCRIPT
# =============================================================================
cat("--- PART A: CONFIGURATION ---\n")
feature_matrix_path <-"C:/Users/deevanshi.walia/Desktop/Try 2.0/Onset of Symptoms/Final_Analysis_Output_Peaks_Clustering_9_0.05/Feature_Matrices/feature_matrix_MEDIAN_INTENSITY.csv"
NUM_FOLDS <- 10
SELECTION_STABILITY_THRESHOLD <- 0.8 # e.g., 0.8 means selected in at least 8/10 folds

# =============================================================================
# PART B: LOAD AND PREPARE DATA
# =============================================================================
cat("\n--- PART B: LOADING DATA ---\n")
data_full <- read_csv(feature_matrix_path, show_col_types = FALSE)
data_full$GroupName <- as.factor(data_full$GroupName)
if (any(is.na(data_full))) {
  warning("NA values detected. Removing rows with NAs.")
  data_full <- na.omit(data_full)
}
cat("✅  Data loaded successfully with", nrow(data_full), "samples.\n")

# =============================================================================
# PART C: CROSS-VALIDATED FEATURE SELECTION (RFE-CV)
# =============================================================================
cat(paste("\n--- PART C: RUNNING", NUM_FOLDS, "-FOLD CROSS-VALIDATED FEATURE SELECTION ---\n"))
set.seed(123)
folds <- createFolds(data_full$GroupName, k = NUM_FOLDS, list = TRUE, returnTrain = FALSE)
list_of_fold_results <- list()

for (i in 1:NUM_FOLDS) {
  train_indices <- setdiff(1:nrow(data_full), folds[[i]])
  train_data <- data_full[train_indices, ]
  
  metadata_train <- train_data %>% dplyr::select(SampleName, GroupName)
  feature_data_train <- train_data %>% dplyr::select(starts_with("Cluster_"))
  
  stats_results_train <- lapply(names(feature_data_train), function(cn) {
    p_val <- tryCatch(wilcox.test(feature_data_train[[cn]] ~ metadata_train$GroupName)$p.value, error=function(e) 1.0)
    return(data.frame(Feature = cn, p.value = p_val))
  })
  stats_df_train <- bind_rows(stats_results_train)
  list_of_fold_results[[i]] <- stats_df_train
}
all_fold_p_values <- bind_rows(list_of_fold_results, .id="Fold")

p_value_summary <- all_fold_p_values %>%
  group_by(Feature) %>%
  summarise(Selection_Count = sum(p.value < 0.05, na.rm = TRUE)) %>%
  arrange(desc(Selection_Count))

cat("--- Feature Selection Stability Summary ---\n")
print.data.frame(p_value_summary)

selected_features <- p_value_summary %>%
  filter(Selection_Count >= (SELECTION_STABILITY_THRESHOLD * NUM_FOLDS)) %>%
  pull(Feature)

if (length(selected_features) == 0) {
  stop("No features were selected consistently enough to meet the threshold. Consider lowering SELECTION_STABILITY_THRESHOLD.")
}
cat("\n--- Final Robust Feature Signature ---\n")
print(selected_features)

# =============================================================================
# PART D: MODEL GAUNTLET - LEAVE-ONE-OUT CROSS-VALIDATION
# =============================================================================
cat("\n--- PART D: RUNNING MODEL GAUNTLET WITH LOOCV ---\n")

# --- Prepare data for the gauntlet ---
metadata <- data_full %>% dplyr::select(SampleName, GroupName)
feature_data_best <- data_full %>% dplyr::select(all_of(selected_features))

num_samples <- nrow(feature_data_best)
y_response <- as.factor(metadata$GroupName)
y_response_numeric <- as.numeric(y_response) - 1
group_levels <- levels(y_response)
positive_class_name <- group_levels[1]

loocv_results <- data.frame(
  SampleName = metadata$SampleName, ActualGroup = y_response,
  lda_scores=numeric(num_samples), rf_scores=numeric(num_samples),
  svm_scores=numeric(num_samples), xgb_scores=numeric(num_samples)
)

cat("    - Running LOOCV for all 4 models...\n")
for (i in 1:num_samples) {
  cat(paste0("\r  - [", i, "/", num_samples, "]")); flush.console()
  
  x_train <- as.matrix(feature_data_best[-i, , drop=FALSE])
  y_train <- y_response[-i]
  y_train_numeric <- y_response_numeric[-i]
  x_test <- as.matrix(feature_data_best[i, , drop = FALSE])
  
  # 1. LDA
  lda_model <- lda(x = x_train, grouping = y_train)
  loocv_results$lda_scores[i] <- predict(lda_model, newdata = data.frame(x_test))$posterior[, 2]
  
  # 2. Random Forest (balanced)
  min_class_size <- min(table(y_train))
  rf_model <- randomForest(x = x_train, y = y_train, ntree = 500, sampsize = c(min_class_size, min_class_size))
  loocv_results$rf_scores[i] <- predict(rf_model, newdata = x_test, type = "prob")[, 2]
  
  # 3. SVM
  svm_model <- svm(x = x_train, y = y_train, probability = TRUE, kernel="radial")
  svm_pred <- predict(svm_model, newdata = data.frame(x_test), probability = TRUE)
  loocv_results$svm_scores[i] <- attr(svm_pred, "probabilities")[, 2]
  
  # 4. XGBoost (balanced)
  dtrain <- xgb.DMatrix(data = x_train, label = y_train_numeric)
  dtest <- xgb.DMatrix(data = x_test)
  scale_pos_weight <- sum(y_train_numeric == 0) / sum(y_train_numeric == 1)
  params <- list(objective = "binary:logistic", eval_metric = "logloss", scale_pos_weight = scale_pos_weight, eta = 0.1, max_depth = 3)
  xgb_model <- xgboost(params = params, data = dtrain, nrounds = 50, verbose = 0)
  loocv_results$xgb_scores[i] <- predict(xgb_model, dtest)
}
cat("\n✅  LOOCV complete for all models.\n")

# =============================================================================
# PART E: PERFORMANCE COMPARISON
# =============================================================================
cat("\n--- PART E: COMPARING MODEL PERFORMANCES ---\n")

roc_lda <- roc(response = loocv_results$ActualGroup, predictor = loocv_results$lda_scores, levels=group_levels)
roc_rf <- roc(response = loocv_results$ActualGroup, predictor = loocv_results$rf_scores, levels=group_levels)
roc_svm <- roc(response = loocv_results$ActualGroup, predictor = loocv_results$svm_scores, levels=group_levels)
roc_xgb <- roc(response = loocv_results$ActualGroup, predictor = loocv_results$xgb_scores, levels=group_levels)

model_performance_list <- list()
models_to_test <- list(LDA = roc_lda, `Random Forest` = roc_rf, SVM = roc_svm, XGBoost = roc_xgb)
cm_list <- list()

for (model_name in names(models_to_test)) {
  current_roc <- models_to_test[[model_name]]
  best_thresh <- coords(current_roc, "best", ret = "threshold")$threshold
  predicted_classes <- ifelse(current_roc$predictor >= best_thresh, group_levels[2], group_levels[1])
  cm_stats <- confusionMatrix(data = factor(predicted_classes, levels = group_levels), reference = current_roc$response, positive = positive_class_name)
  cm_list[[model_name]] <- cm_stats
  model_performance_list[[model_name]] <- c(AUC = auc(current_roc), cm_stats$overall["Accuracy"], cm_stats$byClass[c("Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value", "Balanced Accuracy")])
}
performance_summary <- as.data.frame(do.call(rbind, model_performance_list)) %>%
  tibble::rownames_to_column("Model") %>%
  mutate(across(where(is.numeric), ~ round(., 3)))

cat("\n--- Performance Summary Table (at Optimal Thresholds) ---\n")
print(performance_summary)

cat("\n🎨  Generating combined ROC Curve Plot...\n")
roc_list <- list(LDA = roc_lda, `Random Forest` = roc_rf, SVM = roc_svm, XGBoost = roc_xgb)
roc_plot_combined <- ggroc(roc_list, size = 1) +
  annotate("segment", x = 1, xend = 0, y = 0, yend = 1, color="grey", linetype="dashed") +
  labs(title = "Comparative ROC Curves for Different Models", color = "Model") +
  theme_minimal(base_size = 14) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
print(roc_plot_combined)

best_model_name <- performance_summary$Model[which.max(performance_summary$`Balanced Accuracy`)]
cat("\n--- Detailed Confusion Matrix for the Best Performing Model:", best_model_name, "---\n")
print(cm_list[[best_model_name]])

# =============================================================================
# PART F: FINAL SUMMARY, VISUALIZATIONS, AND FEATURE IMPORTANCE
# =============================================================================
cat("\n--- PART F: FINAL FEATURE CHARACTERIZATION & VISUALIZATION ---\n")

# --- 1. Create a Final Summary Table of the Selected Features ---
cat("\n--- Final Summary Table of Selected Biomarkers ---\n")

# --- a. Get univariate stats (p-value and fold change) for the selected features ---
# This part assumes 'feature_data', 'metadata', and 'selected_features' exist
univariate_stats_df <- bind_rows(
  lapply(
    selected_features,
    function(feature_name, fd, md) {
      group_levels_func <- levels(as.factor(md$GroupName))
      test_result <- wilcox.test(fd[[feature_name]] ~ md$GroupName)
      mean_g1 <- mean(fd[[feature_name]][md$GroupName == group_levels_func[1]]) + 1e-9
      mean_g2 <- mean(fd[[feature_name]][md$GroupName == group_levels_func[2]]) + 1e-9
      log2_fc <- log2(mean_g2 / mean_g1)
      return(data.frame(Feature = feature_name, p.value = test_result$p.value, log2FoldChange = log2_fc))
    },
    fd = feature_data_best,
    md = metadata
  )
)

# --- b. Get importance scores from the BEST performing model ---
# This part assumes 'best_model_name' exists from PART E
cat("    - Extracting importance scores from the best model:", best_model_name, "\n")
feature_data_best <- data_full %>% dplyr::select(all_of(selected_features))
y_response <- as.factor(metadata$GroupName)
y_response_numeric <- as.numeric(y_response) - 1

if (best_model_name == "Random Forest") {
  final_model <- randomForest(x=as.matrix(feature_data_best), y=y_response, importance=TRUE)
  importance_df <- as.data.frame(importance(final_model)) %>%
    tibble::rownames_to_column("Feature") %>% dplyr::select(Feature, Importance = MeanDecreaseGini)
} else if (best_model_name == "XGBoost") {
  dtrain_final <- xgb.DMatrix(data = as.matrix(feature_data_best), label = y_response_numeric)
  final_model <- xgboost(params = list(objective="binary:logistic"), data = dtrain_final, nrounds = 50, verbose=0)
  importance_df <- xgb.importance(model = final_model) %>%
    as.data.frame() %>% dplyr::select(Feature, Importance = Gain)
} else { # LDA is the default
  final_model <- lda(x = feature_data_best, grouping = y_response)
  importance_df <- as.data.frame(final_model$scaling) %>%
    tibble::rownames_to_column("Feature") %>% dplyr::mutate(Importance = abs(LD1)) %>% dplyr::select(Feature, Importance)
}

# --- c. Join everything together into one final table ---
final_summary_table <- univariate_stats_df %>%
  left_join(importance_df, by = "Feature") %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  dplyr::select(Feature, log2FoldChange, p.value, p.adj, Importance) %>%
  arrange(desc(Importance))

print.data.frame(final_summary_table)


# =============================================================================
# --- 2. NEW: GENERATE FEATURE IMPORTANCE BAR PLOT ---
# =============================================================================

cat("\n🎨  Generating feature importance bar plot...\n")

# Ensure the 'Importance' column is numeric for plotting
final_summary_table$Importance <- as.numeric(final_summary_table$Importance)

importance_barplot <- ggplot(
  data = final_summary_table,
  # We use reorder() to automatically sort the bars from most to least important
  aes(x = Importance, y = reorder(Feature, Importance), fill = Importance)
) +
  geom_col() + # Use geom_col for bar charts
  scale_fill_viridis_c(guide = "none") + # Use a nice color scale, remove the legend
  labs(
    title = paste("Feature Importance from", best_model_name, "Model"),
    subtitle = "Based on RFE-CV selected features",
    x = "Importance Score (e.g., Mean Decrease Gini or LD1 coefficient)",
    y = "Feature"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(importance_barplot)

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
# =============================================================================
# --- 4. Pairs Plot for Interactions ---
# =============================================================================

cat("\n--- PART B: Generating pairs plot for all selected features ---\n")

# A pairs plot is most useful if you have between 2 and ~6 features.
if (length(selected_features) >= 2) {
  
  # 1. Select the data for the plot
  # This includes the GroupName and all the selected cluster intensities
  data_for_pairs_plot <- cbind(
    GroupName = metadata$GroupName,
    feature_data_best[, selected_features]
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

####

