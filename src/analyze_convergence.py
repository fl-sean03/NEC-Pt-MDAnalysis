import json
import argparse
import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def analyze_convergence(fragment_metrics_csv, output_dir, run_id, output_json=None):
    """
    Analyzes convergence metrics from fragment metrics data.

    Args:
        fragment_metrics_csv (str): Path to the fragment metrics CSV file.
        output_dir (str): Directory to save output plots.
        run_id (str): Run ID for plot filenames.
        output_json (str, optional): Path to the output JSON file for the summary results.
    """
    print(f"Analyzing convergence from {fragment_metrics_csv}")

    # If output-json is provided, ensure its directory exists
    if output_json:
        output_json_dir = os.path.dirname(output_json)
        os.makedirs(output_json_dir, exist_ok=True)

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Load data
    try:
        df = pd.read_csv(fragment_metrics_csv)
    except FileNotFoundError:
        print(f"Error: Fragment metrics file not found at {fragment_metrics_csv}")
        # If output_json is provided, write an error summary
        if output_json:
            error_summary = {"stage": "Convergence Analysis", "error": f"Fragment metrics file not found at {fragment_metrics_csv}"}
            with open(output_json, 'w') as f:
                json.dump(error_summary, f, indent=4)
        sys.exit(1) # Exit with error code
    except Exception as e:
        print(f"Error loading fragment metrics CSV: {e}")
        # If output_json is provided, write an error summary
        if output_json:
            error_summary = {"stage": "Convergence Analysis", "error": f"Error loading fragment metrics CSV: {e}"}
            with open(output_json, 'w') as f:
                json.dump(error_summary, f, indent=4)
        sys.exit(1) # Exit with error code

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
    plt.title(f"{run_id} Unique Pt Region Coverage Over Time") # Updated title
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, f"{run_id}_unique_pt_region_coverage.png")) # Updated filename
    plt.close()
    print("Unique Pt Region coverage plot saved.") # Updated print message

    # Implement average Rg over time plotting
    print("Calculating average Radius of Gyration over time...") # Updated print message
    average_rg = df.groupby('frame')['radius_gyration'].mean() # Use 'radius_gyration'

    plt.figure(figsize=(10, 6))
    plt.plot(average_rg.index, average_rg.values)
    plt.xlabel("Frame")
    plt.ylabel("Average Radius of Gyration (Å)") # Updated label
    plt.title(f"{run_id} Average Radius of Gyration Over Time") # Updated title
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, f"{run_id}_average_radius_gyration_over_time.png")) # Updated filename
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
        plt.title(f"{run_id} Average MSD Over Time")
        plt.grid(True)
        plt.savefig(os.path.join(output_dir, f"{run_id}_average_msd_over_time.png"))
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
            is_bound_series = fragment_df['min_dist'] < cutoff # Access min_dist by label here
            event_starts = (is_bound_series) & (~is_bound_series.shift(1, fill_value=False))
            current_events += event_starts.sum()

        event_counts.append(current_events)

    plt.figure(figsize=(10, 6))
    plt.plot(sorted(frames), event_counts)
    plt.xlabel("Frame")
    plt.ylabel("Cumulative Binding Event Count")
    plt.title(f"{run_id} Binding Event Count Growth Over Time")
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, f"{run_id}_event_count_growth.png"))
    plt.close()
    print("Event count growth plot saved.")

    # Prepare results for the comprehensive summary
    summary_results = {
        "stage": "Convergence Analysis",
        "convergence_checks_performed": [
            "Unique Pt Region Coverage Over Time",
            "Average Radius of Gyration Over Time",
            "Average MSD Over Time",
            "Binding Event Count Growth Over Time"
        ],
        "generated_plots": [
            os.path.join(output_dir, "unique_pt_region_coverage.png"),
            os.path.join(output_dir, "average_radius_gyration_over_time.png"),
            os.path.join(output_dir, "average_msd_over_time.png"),
            os.path.join(output_dir, "event_count_growth.png")
        ]
    }

    # Save summary results to JSON file if output_json is provided
    if output_json:
        try:
            with open(output_json, 'w') as f:
                json.dump(summary_results, f, indent=4)
            print(f"Exported summary: {output_json}")
        except Exception as e:
            print(f"Error saving summary JSON to {output_json}: {e}")
            # Continue without exiting, as the analysis itself was successful

    # Print summary results as JSON to stdout (still keep for orchestrator)
    print(json.dumps(summary_results))

    print("Convergence analysis complete.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze convergence metrics from fragment metrics data.")
    parser.add_argument("--csv", required=True, type=str, help="Path to the fragment metrics CSV file.")
    parser.add_argument("--outdir", required=True, type=str, help="Directory to save output plots.")
    parser.add_argument("--run-id", required=True, type=str, help="Run ID for plot filenames.") # Added run_id argument
    parser.add_argument(
        '--output-json', required=False, type=str, # Make output-json optional
        help='Path to the output JSON file for the summary results.'
    )

    args = parser.parse_args()

    # Update function call to pass run_id, output_json and use correct argument names
    analyze_convergence(args.csv, args.outdir, args.run_id, args.output_json)