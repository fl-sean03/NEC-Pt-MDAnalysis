import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse

def calculate_unique_facet_coverage(df):
    """
    Calculate the fraction of total unique facets contacted over time.
    """
    # Assuming 'frame' and 'facet_id' columns exist
    # Need to determine the total number of unique facets from the data or context
    # For now, let's assume facet_id is a categorical identifier
    unique_facets_per_frame = df.groupby('frame')['facet_id'].nunique()
    cumulative_unique_facets = unique_facets_per_frame.cumsum()
    total_unique_facets = df['facet_id'].nunique() # This assumes all possible facets appear at some point

    # Calculate fraction
    fraction_coverage = cumulative_unique_facets / total_unique_facets
    return fraction_coverage

def plot_unique_facet_coverage(fraction_coverage, output_path):
    """
    Plot unique-facet coverage vs. time.
    """
    plt.figure(figsize=(10, 6))
    fraction_coverage.plot()
    plt.xlabel('Frame')
    plt.ylabel('Fraction of Unique Facets Contacted')
    plt.title('Unique Facet Coverage Over Time')
    plt.grid(True)
    plt.savefig(output_path)
    plt.close()

def plot_radius_of_gyration_over_time(df, output_path):
    """
    Plot the average Rg of fragments over time.
    Assumes 'frame' and 'radius_of_gyration' columns exist.
    """
    avg_rg_per_frame = df.groupby('frame')['radius_of_gyration'].mean()

    plt.figure(figsize=(10, 6))
    avg_rg_per_frame.plot()
    plt.xlabel('Frame')
    plt.ylabel('Average Radius of Gyration (Rg)')
    plt.title('Average Radius of Gyration Over Time')
    plt.grid(True)
    plt.savefig(output_path)
    plt.close()

def calculate_msd(df):
    """
    Calculate Mean Squared Displacement (MSD) of fragment centers of mass.
    Assumes 'frame', 'fragment_id', 'com_x', 'com_y', 'com_z' columns exist.
    """
    # This is a simplified MSD calculation. A proper one would involve
    # averaging over multiple time origins.
    # For simplicity, let's calculate displacement from the first frame's position.

    msd_data = []
    fragment_groups = df.groupby('fragment_id')
    initial_positions = fragment_groups.head(1).set_index('fragment_id')[['com_x', 'com_y', 'com_z']]

    for frame, frame_df in df.groupby('frame'):
        if frame == 0:
            continue # Skip the first frame as it's the origin

        # Merge current frame positions with initial positions
        merged_df = pd.merge(frame_df, initial_positions, on='fragment_id', suffixes=('_current', '_initial'))

        # Calculate displacement squared for each fragment
        displacement_squared = (
            (merged_df['com_x_current'] - merged_df['com_x_initial'])**2 +
            (merged_df['com_y_current'] - merged_df['com_y_initial'])**2 +
            (merged_df['com_z_current'] - merged_df['com_z_initial'])**2
        )

        # Calculate mean squared displacement for the frame
        msd_frame = displacement_squared.mean()
        msd_data.append({'frame': frame, 'msd': msd_frame})

    msd_df = pd.DataFrame(msd_data)
    return msd_df.set_index('frame')['msd']


def plot_msd(msd_series, output_path):
    """
    Plot Mean Squared Displacement (MSD) over time.
    """
    plt.figure(figsize=(10, 6))
    msd_series.plot()
    plt.xlabel('Frame')
    plt.ylabel('Mean Squared Displacement (MSD)')
    plt.title('Mean Squared Displacement Over Time')
    plt.grid(True)
    plt.savefig(output_path)
    plt.close()

def calculate_event_count_growth(df):
    """
    Track and plot the cumulative event count per facet over time.
    This requires identifying 'events'. An event could be defined as
    a fragment contacting a facet. We need to track the first time
    each fragment contacts each facet.
    Assumes 'frame', 'fragment_id', and 'facet_id' columns exist.
    """
    # Identify unique fragment-facet contacts
    fragment_facet_contacts = df[['frame', 'fragment_id', 'facet_id']].drop_duplicates()

    # Sort by frame to get the first contact time
    fragment_facet_contacts = fragment_facet_contacts.sort_values(by='frame')

    # Get the first frame each fragment-facet pair is observed
    first_contact_frames = fragment_facet_contacts.groupby(['fragment_id', 'facet_id'])['frame'].min().reset_index()

    # Count cumulative unique events (fragment-facet pairs) over time
    # An 'event' is a unique fragment-facet contact occurring at a specific frame
    cumulative_events = first_contact_frames.groupby('frame').size().cumsum()

    # Reindex to include all frames, filling missing frames with the last known cumulative count
    all_frames = df['frame'].unique()
    cumulative_events = cumulative_events.reindex(all_frames, method='ffill').fillna(0)

    return cumulative_events


def plot_event_count_growth(cumulative_events, output_path):
    """
    Plot cumulative event count growth over time.
    """
    plt.figure(figsize=(10, 6))
    cumulative_events.plot()
    plt.xlabel('Frame')
    plt.ylabel('Cumulative Event Count (Fragment-Facet Contacts)')
    plt.title('Cumulative Event Count Growth Over Time')
    plt.grid(True)
    plt.savefig(output_path)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Analyze molecular dynamics simulation data for convergence.')
    parser.add_argument('input_csv', help='Path to the _fragment_metrics.csv file')
    parser.add_argument('output_dir', help='Directory to save output plots')
    args = parser.parse_args()

    input_csv_path = args.input_csv
    output_base_dir = args.output_dir

    # Extract run_id from the input file path
    # Assuming the path is like outputs/<run_id>/<run_id>_fragment_metrics.csv
    try:
        run_id = input_csv_path.split(os.sep)[-2]
    except IndexError:
        print(f"Could not extract run_id from path: {input_csv_path}")
        run_id = "default_run" # Fallback run_id

    output_convergence_dir = os.path.join(output_base_dir, run_id, 'convergence_plots')
    os.makedirs(output_convergence_dir, exist_ok=True)

    print(f"Reading data from {input_csv_path}")
    df = pd.read_csv(input_csv_path)

    # Ensure required columns exist
    required_cols = ['frame', 'fragment_id', 'facet_id', 'radius_of_gyration', 'com_x', 'com_y', 'com_z']
    if not all(col in df.columns for col in required_cols):
        missing = [col for col in required_cols if col not in df.columns]
        print(f"Error: Input CSV is missing required columns: {missing}")
        return

    print("Calculating and plotting Unique-facet coverage...")
    fraction_coverage = calculate_unique_facet_coverage(df)
    plot_unique_facet_coverage(fraction_coverage, os.path.join(output_convergence_dir, f'{run_id}_unique_facet_coverage.png'))
    print("Unique-facet coverage plot saved.")

    print("Calculating and plotting Radius of gyration over time...")
    plot_radius_of_gyration_over_time(df, os.path.join(output_convergence_dir, f'{run_id}_avg_radius_of_gyration.png'))
    print("Radius of gyration plot saved.")

    print("Calculating and plotting MSD...")
    msd_series = calculate_msd(df)
    plot_msd(msd_series, os.path.join(output_convergence_dir, f'{run_id}_msd.png'))
    print("MSD plot saved.")

    print("Calculating and plotting Event-count growth...")
    cumulative_events = calculate_event_count_growth(df)
    plot_event_count_growth(cumulative_events, os.path.join(output_convergence_dir, f'{run_id}_event_count_growth.png'))
    print("Event-count growth plot saved.")

    print("Convergence analysis complete.")

if __name__ == "__main__":
    main()