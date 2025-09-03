# =============================================================================
#
#   --- SCRIPT: CALCULATE INTENSITY-WEIGHTED CLUSTER CENTROIDS ---
#
# Description:
# This script loads an aggregated peak list with cluster assignments (from DBSCAN)
# and calculates the intensity-weighted centroid for each cluster. This provides a
# more robust and representative coordinate for each feature, as it is pulled

# towards the highest-intensity peaks. The results are saved and visualized.
#
# =============================================================================

# --- 1. Load Required Libraries ---
# install.packages(c("dplyr", "readr", "ggplot2", "ggrepel"))
library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel) # For non-overlapping labels in the plot

# =============================================================================
# PART A: CONFIGURE YOUR SCRIPT
# =============================================================================

cat("--- PART A: CONFIGURATION ---\n")

# --- Path to your aggregated peak list ---
# This MUST be the file that contains the final cluster assignments.
# It should have columns like 'RetentionTime', 'DriftTime', 'Intensity', 'Cluster'.
aggregated_peaks_path <- "C:/Users/deevanshi.walia/Desktop/Try 2.0/CT-Value_25CT/Final_Analysis_Output_Peaks_Clustering_9_0.05/Clustering_Results/master_peak_list_with_clusters.csv"

# --- Define the output directory ---
output_dir <- "C:/Users/deevanshi.walia/Desktop/Try 2.0/CT-Value_25CT/Final_Analysis_Output_Peaks_Clustering_9_0.05/Clustering_Results" # Or any other folder
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# PART B: SCRIPT EXECUTION
# =============================================================================

# --- Step 1: Load and Prepare Data ---
cat("\n--- PART B: LOADING AND PREPARING DATA ---\n")
tryCatch({
  all_peaks <- read_csv(aggregated_peaks_path, show_col_types = TRUE)
}, error = function(e) { stop("Could not read the aggregated peaks file. Error: ", e$message) })

# --- Filter out noise points (Cluster 0) ---
# We only want to calculate centroids for the real, numbered clusters.
peaks_no_noise <- all_peaks %>%
  filter(Cluster != 0)

# Check if there are any clusters to process
if (nrow(peaks_no_noise) == 0) {
  stop("No clustered peaks found (all points were classified as noise). Cannot calculate centroids.")
}

cat("✅  Data loaded. Found", length(unique(peaks_no_noise$Cluster)), "clusters to process.\n")


# --- Step 2: Calculate the Intensity-Weighted Centroids ---
cat("\n--- PART C: CALCULATING WEIGHTED CENTROIDS ---\n")

# This is the core calculation using dplyr
weighted_centroids <- peaks_no_noise %>%
  # Group the data by the 'Cluster' column
  group_by(Cluster) %>%
  # For each group, calculate the summary statistics
  summarise(
    # Weighted Retention Time: sum(RT * Intensity) / sum(Intensity)
    weighted_RT = sum(RetentionTime * Intensity, na.rm = TRUE) / sum(Intensity, na.rm = TRUE),
    
    # Weighted Drift Time: sum(DT * Intensity) / sum(Intensity)
    weighted_DT = sum(DriftTime * Intensity, na.rm = TRUE) / sum(Intensity, na.rm = TRUE),
    
    # Also calculate some other useful info
    Total_Intensity = sum(Intensity, na.rm = TRUE),
    Peak_Count = n() # Counts the number of peaks in the cluster
  ) %>%
  # Arrange by cluster number for a clean output
  arrange(Cluster)

cat("✅  Weighted centroids calculated successfully.\n")
cat("--- First few rows of the Centroid Table ---\n")
print(head(weighted_centroids))


# --- Step 3: Save the Centroid List ---
centroid_output_path <- file.path(output_dir, "weighted_centroids.csv")
write_csv(weighted_centroids, centroid_output_path)
cat("\n✅  Centroid list saved to:", centroid_output_path, "\n")


# --- Step 4: Create a Validation Plot ---
cat("\n--- PART D: GENERATING VALIDATION PLOT ---\n")

centroid_plot <- ggplot() +
  # Layer 1: Plot all the individual peaks, colored by their cluster
  geom_point(
    data = peaks_no_noise,
    aes(x = DriftTime, y = RetentionTime, color = as.factor(Cluster)),
    alpha = 0.5,
    size = 1.5
  ) +
  
  # Layer 2: Overlay the calculated weighted centroids as black stars
  geom_point(
    data = weighted_centroids,
    aes(x = weighted_DT, y = weighted_RT),
    color = "black",
    shape = 8, # An asterisk shape
    size = 4,
    stroke = 1
  ) +
  
  # Layer 3: Add labels to the centroids
  geom_text_repel(
    data = weighted_centroids,
    aes(x = weighted_DT, y = weighted_RT, label = Cluster),
    color = "black",
    fontface = "bold",
    box.padding = 0.5
  ) +
  
  labs(
    title = "DBSCAN Clusters with Intensity-Weighted Centroids",
    x = "Drift Time (ms)",
    y = "Retention Time (s)"
  ) +
  
  theme_minimal() +
  guides(color = "none") # Hide the legend for the individual points to reduce clutter

# Print the plot and save it
print(centroid_plot)
plot_output_path <- file.path(output_dir, "cluster_centroids_visualization.png")
ggsave(plot_output_path, plot = centroid_plot, width = 12, height = 9, dpi = 300, bg = "white")

cat("\n✅  Centroid validation plot saved to:", plot_output_path, "\n")
cat("\n🎉 --- Workflow Complete! --- 🎉\n")

