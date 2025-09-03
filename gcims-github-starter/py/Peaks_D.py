# -----------------------------------------------------------------------------
# Script for 2D Peak Identification on Aligned and RIP-Removed GC-IMS Data
# -----------------------------------------------------------------------------
import pandas as pd
import matplotlib.pyplot as plt
from skimage.feature import peak_local_max
import numpy as np
import os

# --- 1. Define Input and Output Directories ---

# Directory with your ALIGNED and RIP-REMOVED sample CSVs from the R script.
input_dir = r"C:\Users\deevanshi.walia\Desktop\Try 2.0\PCR\Preprocessed_AA_All_Samples\Divided_Samples\PCR_Positive"

# A new directory to save the final peak lists and plots.
output_dir = r"C:\Users\deevanshi.walia\Desktop\Try 2.0\PCR_Pos_Peaks_9_0.05"

# --- 2. Create Output Directory ---
os.makedirs(output_dir, exist_ok=True)

# --- 3. Find all CSV files to process ---
try:
    csv_files = [f for f in os.listdir(input_dir) if f.endswith(".csv")]
    if not csv_files:
        raise FileNotFoundError
except FileNotFoundError:
    print(f"Error: No CSV files found in the specified input directory:\n{input_dir}")
    exit()

print(f" Found {len(csv_files)} processed matrices to process for peak picking.\n")


# --- 4. Loop Through Each Processed Matrix ---
for file_name in csv_files:
    print(f"  Processing: {file_name}")
    file_path = os.path.join(input_dir, file_name)

    # --- Read the data, using comma as the separator ---
    df = pd.read_csv(file_path, sep=',')

    # --- Add a safety check for empty or 1D data ---
    if df.shape[0] < 2 or df.shape[1] < 3:
        print(f"    -  WARNING: Not enough data in '{file_name}' to process. Skipping.")
        continue

    # --- Separate coordinates from the intensity matrix ---
    drift_times = df.iloc[:, 0].values
    intensity_matrix = df.iloc[:, 1:].values
    retention_times = df.columns[1:].astype(float)

    # --- Peak Picking using scikit-image ---
    # The image is transposed so that the output coordinates are (retention_time_idx, drift_time_idx)
    image_for_peaks = intensity_matrix.T
    
    peaks = peak_local_max(
        image_for_peaks, 
        min_distance=9,      # Tune this
        threshold_abs=0.05 # Tune this
    )

    print(f"    - Detected {len(peaks)} peaks.")
    
    if len(peaks) > 0:
        # --- Store peaks in a new DataFrame using the real coordinates ---
        peak_df = pd.DataFrame({
            'RetentionTime': retention_times[peaks[:, 0]],
            'DriftTime': drift_times[peaks[:, 1]],
            'Intensity': image_for_peaks[peaks[:, 0], peaks[:, 1]]
        })

        # --- Save the list of peaks to a new CSV file ---
        base_name = file_name.replace('_processed.csv', '')
        peak_output_path = os.path.join(output_dir, f"peaks_{base_name}.csv")
        peak_df.to_csv(peak_output_path, index=False, sep=',') # Use comma separator
        print(f"    - Saved peak list to: {os.path.basename(peak_output_path)}")

        # --- Create and Save a Publication-Ready Heatmap ---
        plt.figure(figsize=(12, 8))
        
        # Use imshow with the TRANSPOSED matrix and corrected 'extent' for your desired orientation
        # extent=[x_min, x_max, y_min, y_max] -> [DriftTime_min, DriftTime_max, RetentionTime_min, RetentionTime_max]
        plt.imshow(
            intensity_matrix.T,
            cmap='plasma', # Using 'plasma' to match your R plots
            interpolation='nearest', 
            aspect='auto',
            origin='lower',
            extent=[drift_times[0], drift_times[-1], retention_times[0], retention_times[-1]]
        )
        
        plt.colorbar(label='Intensity')

        # --- Plot the detected peaks using their real DriftTime and RetentionTime values ---
        plt.plot(
            peak_df['DriftTime'],      # X-coordinate
            peak_df['RetentionTime'],  # Y-coordinate
            'ro',                      # Red circles
            markersize=4, 
            markerfacecolor='none', 
            markeredgewidth=1
        )

        plt.title(f"Aligned Heatmap with Peaks - {base_name}")
        plt.xlabel("Drift Time (ms)")
        plt.ylabel("Retention Time (s)")
        plt.tight_layout()
        
        # Save the plot
        plot_output_path = os.path.join(output_dir, f"heatmap_{base_name}.png")
        plt.savefig(plot_output_path, dpi=300)
        plt.close() # Close the figure to free up memory
        print(f"    - Saved heatmap plot to: {os.path.basename(plot_output_path)}\n")

    else:
        print("    - No peaks found with the current settings.\n")

print("🎉 --- Peak picking complete for all files. ---")