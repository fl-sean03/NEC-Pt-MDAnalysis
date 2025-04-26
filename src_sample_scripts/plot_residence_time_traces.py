#!/usr/bin/env python3
"""
plot_residence_time_traces.py: Generate 3D trajectory traces for selected molecules, colored by facet ID.
"""
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm

def plot_traces(df, outdir, run_id, fragment_ids=None):
    """
    Generate 3D trajectory traces for specified fragment IDs, colored by facet ID.
    If fragment_ids is None, plots traces for all fragments.
    """
    os.makedirs(outdir, exist_ok=True)

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    if fragment_ids is None:
        fragment_ids = df['fragment_id'].unique()

    # Get unique facet regions for coloring
    regions = sorted(df['nearest_pt_class'].dropna().unique())
    num_regions = len(regions)
    colors = cm.viridis(np.linspace(0, 1, num_regions))
    region_color_map = {region: colors[i] for i, region in enumerate(regions)}

    for frag_id in fragment_ids:
        frag_df = df[df['fragment_id'] == frag_id].sort_values('time_ps').reset_index(drop=True)

        # Plot segments colored by facet ID
        for i in range(len(frag_df) - 1):
            p1 = frag_df.iloc[i]
            p2 = frag_df.iloc[i+1]
            region = p1['nearest_pt_class'] # Color segment based on the starting point's region

            color = region_color_map.get(region, 'gray') # Default to gray if region is NaN or unknown

            ax.plot([p1['com_x'], p2['com_x']], [p1['com_y'], p2['com_y']], [p1['com_z'], p2['com_z']], color=color, alpha=0.7)

    ax.set_xlabel('COM X (Å)')
    ax.set_ylabel('COM Y (Å)')
    ax.set_zlabel('COM Z (Å)')
    ax.set_title(f'{run_id} Residence Time Traces (Colored by Facet ID)')

    # Create a legend for facet colors
    legend_elements = [plt.Line2D([0], [0], color=region_color_map[region], lw=4, label=region) for region in regions]
    ax.legend(handles=legend_elements, title="Facet ID")


    fname = f"{run_id}_residence_time_traces_facet_colored.png"
    outpath = os.path.join(outdir, fname)
    plt.savefig(outpath)
    plt.close()
    print(f'Saved residence time traces plot: {outpath}')


def main():
    parser = argparse.ArgumentParser(
        description='Generate 3D trajectory traces for selected molecules, colored by facet ID.'
    )
    parser.add_argument(
        '--csv', required=True,
        help='Path to fragment metrics CSV (e.g. 0HPt_fragment_metrics.csv)'
    )
    parser.add_argument(
        '--outdir', default='visualization_plots',
        help='Directory to save output plots'
    )
    parser.add_argument(
        '--fragment-ids', nargs='+', type=int,
        help='List of fragment IDs to plot (e.g., 1 5 10). Plots all if not specified.'
    )
    args = parser.parse_args()

    # Load data
    df = pd.read_csv(args.csv)

    # Infer run_id
    run_id = df['run_id'].iloc[0] if 'run_id' in df.columns else os.path.splitext(os.path.basename(args.csv))[0]

    # Construct the output directory based on run_id
    base_output_dir = "2-1250Run/outputs"
    run_output_dir = os.path.join(base_output_dir, run_id)
    plot_output_dir = os.path.join(run_output_dir, args.outdir)

    # Plot residence time traces
    plot_traces(df, plot_output_dir, run_id, args.fragment_ids)


if __name__ == '__main__':
    main()