#!/usr/bin/env python3
"""
Regenerated script to generate 3D trajectory traces for selected molecules, colored by Pt region.
"""
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import sys
from scipy.spatial.distance import cdist

def calculate_spatial_excursion(df_fragment: pd.DataFrame) -> float:
    """
    Calculate the maximum spatial excursion (max pairwise distance) for a fragment.

    Args:
        df_fragment (pd.DataFrame): DataFrame containing 'com_x', 'com_y', 'com_z' for a single fragment.

    Returns:
        float: Maximum pairwise distance between COM positions. Returns 0.0 if less than 2 points.
    """
    if len(df_fragment) < 2:
        return 0.0
    coords = df_fragment[['com_x', 'com_y', 'com_z']].values
    # Calculate pairwise distances and return the maximum
    return cdist(coords, coords).max()

def calculate_residence_events(df_fragment: pd.DataFrame, min_duration_frames: int) -> int:
    """
    Calculate the number of residence events for a fragment based on contiguous regions.

    Args:
        df_fragment (pd.DataFrame): DataFrame containing 'nearest_pt_class' for a single fragment.
        min_duration_frames (int): Minimum duration of a residence event in frames.

    Returns:
        int: The number of residence events.
    """
    if df_fragment.empty:
        return 0

    # Identify consecutive identical regions (including NaN)
    regions = df_fragment['nearest_pt_class'].fillna('NaN').tolist()
    event_count = 0
    current_event_duration = 0
    current_region = None

    for region in regions:
        if region == current_region:
            current_event_duration += 1
        else:
            # Check if the previous block was a valid residence event (not NaN and long enough)
            if current_region != 'NaN' and current_event_duration >= min_duration_frames:
                event_count += 1
            current_region = region
            current_event_duration = 1 # Start new block duration

    # Check the last block
    if current_region != 'NaN' and current_event_duration >= min_duration_frames:
        event_count += 1

    return event_count


def calculate_unique_regions(df_fragment: pd.DataFrame) -> int:
    """
    Calculate the number of unique Pt regions a fragment visited.

    Args:
        df_fragment (pd.DataFrame): DataFrame containing 'nearest_pt_class' for a single fragment.

    Returns:
        int: The number of unique non-NaN Pt regions.
    """
    return df_fragment['nearest_pt_class'].dropna().nunique()

def plot_single_trace(df_single_fragment: pd.DataFrame, df_pt_atoms: pd.DataFrame, output_path: str, title: str, region_color_map: dict, regions: list, distance_cutoff: float):
    """
    Generate a 3D trajectory trace for a single fragment and plot Pt atoms, colored by Pt region.

    Args:
        df_single_fragment (pd.DataFrame): DataFrame containing data for a single fragment.
        df_pt_atoms (pd.DataFrame): DataFrame containing Pt atom coordinates and regions
                                    with 'x', 'y', 'z', and 'region' columns.
        output_path (str): Full path to save the output plot.
        title (str): Title for the plot.
        region_color_map (dict): Mapping from region name to color.
        regions (list): Sorted list of unique Pt regions.
    """
    if df_single_fragment.empty:
        print(f"Warning: No data to plot for {title}. Skipping.")
        return

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Plot segments colored by facet ID
    for i in range(len(df_single_fragment) - 1):
        p1 = df_single_fragment.iloc[i]
        p2 = df_single_fragment.iloc[i+1]
        # Determine color based on distance and region
        if p1['min_dist'] <= distance_cutoff:
            region = p1['nearest_pt_class']
            color = region_color_map.get(region, 'gray') # Use gray for unknown/NaN regions within cutoff
            print(f"Frame: {p1['frame']}, Fragment: {p1['fragment_id']}, MinDist: {p1['min_dist']:.2f}, NearestRegion: {region}, AssignedColor: {color}") # Debug print
        else:
            color = 'gray' # Use gray if far from Pt surface

        ax.plot([p1['com_x'], p2['com_x']], [p1['com_y'], p2['com_y']], [p1['com_z'], p2['com_z']], color=color, alpha=0.7, linewidth=1)

    # Plot Pt atoms
    if not df_pt_atoms.empty:
        for region, color in region_color_map.items():
            df_region = df_pt_atoms[df_pt_atoms['region'] == region]
            if not df_region.empty:
                ax.scatter(df_region['x'], df_region['y'], df_region['z'], color=color, s=10, label=f'Pt ({region})', alpha=0.5)
        # Plot NaN/Unknown Pt atoms if any
        df_unknown_pt = df_pt_atoms[df_pt_atoms['region'].isnull()]
        if not df_unknown_pt.empty:
             ax.scatter(df_unknown_pt['x'], df_unknown_pt['y'], df_unknown_pt['z'], color='gray', s=10, label='Pt (Unknown/NaN)', alpha=0.5)


    # Calculate total distance traveled
    total_distance = 0.0
    for i in range(len(df_single_fragment) - 1):
        p1 = df_single_fragment.iloc[i]
        p2 = df_single_fragment.iloc[i+1]
        distance = np.sqrt((p2['com_x'] - p1['com_x'])**2 + (p2['com_y'] - p1['com_y'])**2 + (p2['com_z'] - p1['com_z'])**2)
        total_distance += distance

    # Calculate number of unique regions contacted within cutoff
    contacted_regions = set()
    for i in range(len(df_single_fragment)):
        p1 = df_single_fragment.iloc[i]
        if p1['min_dist'] <= distance_cutoff and pd.notna(p1['nearest_pt_class']):
            contacted_regions.add(p1['nearest_pt_class'])
    num_contacted_regions = len(contacted_regions)

    # Mark start and end points
    start_point = df_single_fragment.iloc[0]
    end_point = df_single_fragment.iloc[-1]
    ax.scatter(start_point['com_x'], start_point['com_y'], start_point['com_z'], color='green', s=50, label='Start', marker='o')
    ax.scatter(end_point['com_x'], end_point['com_y'], end_point['com_z'], color='red', s=50, label='End', marker='s')

    # Add text annotations for metrics
    # Position the text relative to the axes for now
    # Add text annotations for metrics
    # Position the text relative to the axes
    ax.text2D(1.08, 0.85, f'Total Distance: {total_distance:.2f} Å', transform=ax.transAxes, fontsize=10)
    ax.text2D(1.08, 0.80, f'Contacted Regions: {num_contacted_regions}', transform=ax.transAxes, fontsize=10)


    ax.set_xlabel('COM X (Å)')
    ax.set_ylabel('COM Y (Å)')
    ax.set_zlabel('COM Z (Å)')
    ax.set_title(title)

    # Create legend for fragment trace colors
    fragment_legend_elements = [plt.Line2D([0], [0], color=region_color_map[region], lw=4, label=region) for region in regions if region != 'NaN']
    if 'NaN' in regions or (df_single_fragment['min_dist'] > distance_cutoff).any():
         fragment_legend_elements.append(plt.Line2D([0], [0], color='gray', lw=4, label='Fragment (Far/Unknown)'))

    # Create legend for Pt atom colors/markers
    pt_legend_elements = [plt.Line2D([0], [0], marker='o', color='w', label=f'Pt ({region})',
                                     markerfacecolor=region_color_map.get(region, 'gray'), markersize=8, linestyle='None')
                          for region in regions if region != 'NaN']
    if 'NaN' in regions:
         pt_legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', label='Pt (Unknown/NaN)',
                                              markerfacecolor='gray', markersize=8, linestyle='None'))

    # Create legend for Start/End markers
    start_end_legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', label='Start', markerfacecolor='green', markersize=8, linestyle='None'),
        plt.Line2D([0], [0], marker='s', color='w', label='End', markerfacecolor='red', markersize=8, linestyle='None')
    ]

    # Position the legends
    # Fragment legend
    legend1 = ax.legend(handles=fragment_legend_elements, title="Fragment Region", loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0.)
    ax.add_artist(legend1) # Add the first legend manually to the axes

    # Pt atom legend
    legend2 = ax.legend(handles=pt_legend_elements, title="Pt Region", loc='upper left', bbox_to_anchor=(1.05, 0.7), borderaxespad=0.)
    ax.add_artist(legend2) # Add the second legend manually

    # Start/End legend
    legend3 = ax.legend(handles=start_end_legend_elements, title="Trace Markers", loc='upper left', bbox_to_anchor=(1.05, 0.4), borderaxespad=0.)
    ax.add_artist(legend3) # Add the third legend manually


    plt.tight_layout(rect=[0, 0, 0.85, 1]) # Adjust layout to make space for the legends

    try:
        plt.savefig(output_path)
        print(f'Saved residence time traces plot: {output_path}')
    except Exception as e:
        print(f"Error saving plot to {output_path}: {e}")
    finally:
        plt.close(fig) # Close the figure to free up memory


import MDAnalysis as mda
import json

def main():
    parser = argparse.ArgumentParser(
        description='Generate 3D trajectory traces for selected molecules, colored by Pt region, and visualize Pt atoms.'
    )
    parser.add_argument(
        '--csv', required=True,
        help='Path to fragment metrics CSV (e.g. outputs/run_id/fragment_metrics.csv)'
    )
    parser.add_argument(
        '--psf', required=True,
        help='Path to the PSF file.'
    )
    parser.add_argument(
        '--dcd', required=True,
        help='Path to the DCD trajectory file.'
    )
    parser.add_argument(
        '--pt-json', required=True,
        help='Path to the Pt classification JSON file (e.g. outputs/run_id/pt_classification.json)'
    )
    parser.add_argument(
        '--outdir', default='visualization_plots',
        help='Directory to save output plots within the run_id directory structure'
    )
    parser.add_argument(
        '--min-residence-duration-ns', type=float, default=0.1,
        help='Minimum duration in nanoseconds for a continuous residence event (default: 0.1)'
    )
    parser.add_argument(
        '--frame-time-step-ns', type=float, required=True,
        help='Time step between frames in nanoseconds (required to convert duration to frames)'
    )
    parser.add_argument(
        '--distance-cutoff', type=float, default=3.0,
        help='Maximum distance from Pt surface for fragment trace to be colored by region (default: 3.0 Å)'
    )
    args = parser.parse_args()

    # Load data
    try:
        df = pd.read_csv(args.csv)
    except FileNotFoundError:
        print(f"Error: Input fragment metrics CSV file not found at {args.csv}")
        sys.exit(1)
    except Exception as e:
        print(f"Error loading fragment metrics CSV file: {e}")
        sys.exit(1)

    print("\n--- DataFrame Info ---")
    print("Columns:", df.columns.tolist())
    print("----------------------\n")

    # Check min_dist values for selected fragments
    print("\n--- MinDist for Selected Fragments ---")
    for frag_id, selection_type in selected_fragments.items():
        if frag_id in df['fragment_id'].unique():
            df_frag = df[df['fragment_id'] == frag_id]
            print(f"Fragment ID {frag_id} ({selection_type}):")
            print(df_frag[['frame', 'min_dist', 'nearest_pt_class']].to_string())
            print("-" * 20)
        else:
            print(f"Fragment ID {frag_id} ({selection_type}) not found in DataFrame.")
    print("------------------------------------\n")


    # Load universe to get Pt atom coordinates
    try:
        u = mda.Universe(args.psf, args.dcd)
        pt_atoms = u.select_atoms("name PT*")
        if not pt_atoms:
            print("Error: No Pt atoms found in the universe.")
            sys.exit(1)
        pt_coords = pt_atoms.positions
    except FileNotFoundError as e:
        print(f"Error: Input trajectory file not found: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error loading universe or selecting Pt atoms: {e}")
        sys.exit(1)

    # Load Pt classification
    try:
        with open(args.pt_json, 'r') as f:
            pt_classification = json.load(f)
        # Ensure keys are integers
        pt_classification = {int(k): v for k, v in pt_classification.items()}
    except FileNotFoundError:
        print(f"Error: Pt classification JSON file not found at {args.pt_json}")
        sys.exit(1)
    except json.JSONDecodeError:
        print(f"Error: Could not decode Pt classification JSON file at {args.pt_json}")
        sys.exit(1)
    except Exception as e:
        print(f"Error loading Pt classification JSON file: {e}")
        sys.exit(1)

    # Create DataFrame for Pt atoms
    pt_data = []
    for i, atom in enumerate(pt_atoms):
        region = pt_classification.get(atom.index, None) # Use atom.index for lookup
        pt_data.append({'x': atom.position[0], 'y': atom.position[1], 'z': atom.position[2], 'region': region})
    df_pt_atoms = pd.DataFrame(pt_data)


    # Infer run_id from the directory structure
    # Assuming CSV path is like outputs/<run_id>/fragment_metrics.csv
    run_id = os.path.basename(os.path.dirname(args.csv))
    if not run_id or run_id == 'outputs':
         # Fallback if the structure is not as expected
         run_id = df['run_id'].iloc[0] if 'run_id' in df.columns else os.path.splitext(os.path.basename(args.csv))[0]
         print(f"Inferred run_id: {run_id}")

    # Construct the output directory based on run_id and the specified outdir
    # Assuming output structure like outputs/<run_id>/visualization_plots
    base_output_dir = os.path.dirname(os.path.dirname(args.csv)) # Go up two levels from CSV path
    if os.path.basename(base_output_dir) != "outputs":
         # Adjust if CSV path is not in expected outputs/<run_id>/ structure
         base_output_dir = "outputs" # Default to 'outputs' if structure is unexpected
         run_output_dir = os.path.join(base_output_dir, run_id)
    else:
        run_output_dir = os.path.join(base_output_dir, run_id)

    plot_output_dir = os.path.join(run_output_dir, args.outdir)
    os.makedirs(plot_output_dir, exist_ok=True)

    # Calculate metrics for each fragment
    fragment_metrics = {}
    min_duration_frames = int(args.min_residence_duration_ns / args.frame_time_step_ns)

    for frag_id, df_fragment in df.groupby('fragment_id'):
        spatial_excursion = calculate_spatial_excursion(df_fragment)
        residence_events = calculate_residence_events(df_fragment, min_duration_frames)
        unique_regions = calculate_unique_regions(df_fragment)
        fragment_metrics[frag_id] = {
            'spatial_excursion': spatial_excursion,
            'residence_events': residence_events,
            'unique_regions': unique_regions
        }

    metrics_df = pd.DataFrame.from_dict(fragment_metrics, orient='index')
    print("\n--- Fragment Metrics Summary ---")
    print(metrics_df.to_string()) # Use to_string() to ensure full DataFrame is printed
    print("------------------------------\n")

    # Identify top molecules
    most_moved_frag_id = metrics_df['spatial_excursion'].idxmax() if not metrics_df.empty else None
    most_events_frag_id = metrics_df['residence_events'].idxmax() if not metrics_df.empty else None
    most_unique_regions_frag_id = metrics_df['unique_regions'].idxmax() if not metrics_df.empty else None

    selected_fragments = {}
    if most_moved_frag_id is not None:
        selected_fragments[most_moved_frag_id] = 'most_moved'
    if most_events_frag_id is not None and most_events_frag_id not in selected_fragments:
         selected_fragments[most_events_frag_id] = 'most_events'
    if most_unique_regions_frag_id is not None and most_unique_regions_frag_id not in selected_fragments:
         selected_fragments[most_unique_regions_frag_id] = 'most_unique_regions'

    if not selected_fragments:
        print("No fragments found or metrics could not be calculated.")
        sys.exit(1)
    # Include NaN for consistent coloring
    all_regions = df['nearest_pt_class'].unique()
    # Sort non-NaN regions, add NaN at the end if present
    regions = sorted([r for r in all_regions if pd.notna(r)])
    if pd.isna(all_regions).any():
        regions.append('NaN') # Represent NaN as a string key

    num_regions = len(regions)
    cmap = matplotlib.colormaps.get_cmap('tab10')
    colors = [cmap(i) for i in range(num_regions)]
    # Create color map, explicitly mapping 'NaN' to gray
    region_color_map = {region: colors[i] if region != 'NaN' else 'gray' for i, region in enumerate(regions)}


    # Plot traces for selected fragments
    for frag_id, selection_type in selected_fragments.items():
        df_single_fragment = df[df['fragment_id'] == frag_id].sort_values('frame').reset_index(drop=True)
        output_fname = f"{run_id}_fragment_{frag_id}_{selection_type}_trace_with_pt.png"
        output_path = os.path.join(plot_output_dir, output_fname)
        title = f'{run_id} Fragment {frag_id} Trace ({selection_type.replace("_", " ").title()}) with Pt Atoms'
        plot_single_trace(df_single_fragment, df_pt_atoms, output_path, title, region_color_map, regions, args.distance_cutoff)


    print("Residence time traces plotting complete for selected fragments with Pt atoms.")

if __name__ == '__main__':
    main()