import argparse
import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def analyze_convergence(fragment_metrics_csv, output_dir):
    """
    Analyzes convergence metrics from fragment metrics data.

    Args:
        fragment_metrics_csv (str): Path to the fragment metrics CSV file.
        output_dir (str): Directory to save output plots.
    """
    print(f"Analyzing convergence from {fragment_metrics_csv}")

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Load data
    try:
        df = pd.read_csv(fragment_metrics_csv)
    except FileNotFoundError:
        print(f"Error: Fragment metrics file not found at {fragment_metrics_csv}")
        return
    except Exception as e:
        print(f"Error loading fragment metrics CSV: {e}")
        return

    print(f"Debug: Columns in df after loading CSV: {df.columns}") # Added debug print

    # Implement unique-facet coverage calculation and plotting
    print("Calculating unique-facet coverage...")
    unique_facets = []
    all_facets = df['nearest_pt_class'].dropna().unique() # Use nearest_pt_class and drop NaN
    frames = df['frame'].unique()

    for frame in sorted(frames):
        facets_up_to_frame = df[df['frame'] <= frame]['nearest_pt_class'].dropna().unique() # Use nearest_pt_class and drop NaN
        unique_facets.append(len(facets_up_to_frame))

    plt.figure(figsize=(10, 6))
    plt.plot(sorted(frames), unique_facets)
    plt.xlabel("Frame")
    plt.ylabel("Number of Unique Pt Regions Covered") # Updated label
    plt.title("Unique Pt Region Coverage Over Time") # Updated title
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, "unique_pt_region_coverage.png")) # Updated filename
    plt.close()
    print("Unique Pt Region coverage plot saved.") # Updated print message

    # Implement average Rg over time plotting
    print("Calculating average Radius of Gyration over time...") # Updated print message
    average_rg = df.groupby('frame')['radius_gyration'].mean() # Use 'radius_gyration'

    plt.figure(figsize=(10, 6))
    plt.plot(average_rg.index, average_rg.values)
    plt.xlabel("Frame")
    plt.ylabel("Average Radius of Gyration (Å)") # Updated label
    plt.title("Average Radius of Gyration Over Time") # Updated title
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, "average_radius_gyration_over_time.png")) # Updated filename
    plt.close()
    print("Average Radius of Gyration over time plot saved.") # Updated print message

    # Implement MSD calculation and plotting
    print("Calculating MSD over time...")
    msd_data = []
    for fragment_id in df['fragment_id'].unique():
        fragment_df = df[df['fragment_id'] == fragment_id].sort_values('frame')
        if not fragment_df.empty:
            # Access the first row and then select the columns
            first_row = fragment_df.iloc[0]
            initial_com = first_row[['com_x', 'com_y', 'com_z']].values # Use lowercase column names
            displacements = fragment_df[['com_x', 'com_y', 'com_z']].values - initial_com # Use lowercase column names
            squared_displacements = np.sum(displacements**2, axis=1)
            msd_data.append(pd.DataFrame({'frame': fragment_df['frame'], 'msd': squared_displacements}))

    if msd_data:
        all_msd = pd.concat(msd_data)
        average_msd = all_msd.groupby('frame')['msd'].mean()

        plt.figure(figsize=(10, 6))
        plt.plot(average_msd.index, average_msd.values)
        plt.xlabel("Frame")
        plt.ylabel("Mean Squared Displacement (MSD)")
        plt.title("Average MSD Over Time")
        plt.grid(True)
        plt.savefig(os.path.join(output_dir, "average_msd_over_time.png"))
        plt.close()
        print("Average MSD over time plot saved.")
    else:
        print("No fragment data available to calculate MSD.")


    # Implement event-count growth calculation and plotting
    print("Calculating event count growth over time...")
    event_counts = []
    frames = df['frame'].unique()
    cutoff = 3.5 # Assuming a default cutoff, this should ideally come from input parameters

    for frame in sorted(frames):
        df_up_to_frame = df[df['frame'] <= frame]
        current_events = 0
        # Iterate through fragments and their data up to the current frame
        for fragment_id in df_up_to_frame['fragment_id'].unique():
            fragment_df = df_up_to_frame[df_up_to_frame['fragment_id'] == fragment_id].sort_values('frame')
            if fragment_df.empty:
                continue

            # Identify continuous binding events (distance < cutoff) by iterating through rows
            is_bound_series = fragment_df['min_distance'] < cutoff # Access min_distance by label here
            event_starts = (is_bound_series) & (~is_bound_series.shift(1, fill_value=False))
            current_events += event_starts.sum()

        event_counts.append(current_events)

    plt.figure(figsize=(10, 6))
    plt.plot(sorted(frames), event_counts)
    plt.xlabel("Frame")
    plt.ylabel("Cumulative Binding Event Count")
    plt.title("Binding Event Count Growth Over Time")
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, "event_count_growth.png"))
    plt.close()
    print("Event count growth plot saved.")

    print("Convergence analysis complete.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze convergence metrics from fragment metrics data.")
    parser.add_argument("fragment_metrics_csv", type=str, help="Path to the fragment metrics CSV file.")
    parser.add_argument("output_dir", type=str, help="Directory to save output plots.")

    args = parser.parse_args()

    analyze_convergence(args.fragment_metrics_csv, args.output_dir)