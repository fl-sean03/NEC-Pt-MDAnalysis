#!/usr/bin/env python3
"""
Regenerated script to compute residence time histograms for fragment contacts per Pt region.
"""
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import itertools
import json # Added import for json
from residence_analysis_utils import find_residence_events_per_region # Import the new utility function

def compute_dt_ps(df):
    """
    Determine sampling interval (dt) in picoseconds from the time_ps column.
    """
    times = np.sort(df['time_ps'].unique())
    if len(times) < 2:
        raise ValueError('Not enough time points to determine dt')
    dt = np.median(np.diff(times))
    return float(dt)

# Removed the old compute_residence_times function

def plot_histograms(durations, outdir, run_id, max_bin=None, bins=50):
    """
    Plot and save a histogram of residence times for each Pt region.
    """
    os.makedirs(outdir, exist_ok=True)
    for region, times in durations.items():
        if not times:
            print(f'No contact/residence events to plot for region: {region}')
            continue # Corrected continue statement
        arr = np.array(times)
        mb = max_bin if max_bin is not None else arr.max()
        edges = np.linspace(0, mb, bins + 1)
        plt.figure()
        plt.hist(arr, bins=edges, color='C0', edgecolor='black')
        plt.xlabel('Residence time (ns)')
        plt.ylabel('Count')
        plt.title(f'{run_id} residence times\nregion: {region} (n={len(arr)})')
        plt.tight_layout()
        fname = f"{run_id}_{region}_residence_time_hist.png"
        outpath = os.path.join(outdir, fname)
        plt.savefig(outpath)
        plt.close()
        print(f'Saved {outpath}')
    return

def plot_overlay(durations, outdir, run_id, max_bin=None, bins=50):
    """
    Plot and save an overlaid histogram of residence times for all Pt regions.
    """
    os.makedirs(outdir, exist_ok=True)
    all_times = [t for times in durations.values() for t in times]
    if not all_times:
        print('No residence events to plot in overlay')
        return
    mb = max_bin if max_bin is not None else max(all_times)
    edges = np.linspace(0, mb, bins + 1)
    plt.figure(figsize=(8, 5))
    cmap = plt.get_cmap('tab10')
    regions = list(durations.keys())
    for idx, region in enumerate(regions):
        times = durations[region]
        if not times:
            continue
        arr = np.array(times)
        color = cmap(idx % cmap.N)
        plt.hist(arr, bins=edges, histtype='step', linewidth=2,
                 color=color, label=f"{region} (n={len(arr)})")
    plt.xlabel('Residence time (ns)')
    plt.ylabel('Count')
    plt.title(f'{run_id} residence times overlay')
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0)
    plt.tight_layout(rect=[0, 0, 0.75, 1])
    fname = f"{run_id}_residence_time_overlay.png"
    outpath = os.path.join(outdir, fname)
    plt.savefig(outpath)
    plt.close()
    print(f'Saved overlay {outpath}')

def plot_side_by_side(durations, outdir, run_id, max_bin=None, bins=50):
    """
    Plot four histograms side by side for each Pt region.
    """
    os.makedirs(outdir, exist_ok=True)
    regions = sorted(durations.keys())
    all_times = [t for times in durations.values() for t in times]
    if not all_times:
        print('No residence events to plot side by side')
        return
    mb = max_bin if max_bin is not None else max(all_times)
    edges = np.linspace(0, mb, bins + 1)
    counts_per_region = {}
    for region in regions:
        times = durations.get(region, [])
        if times:
            counts, _ = np.histogram(times, bins=edges)
        else:
            counts = np.zeros(len(edges) - 1, dtype=int)
        counts_per_region[region] = counts
    y_max = max(counts.max() for counts in counts_per_region.values())
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=False, sharey=False)
    cmap = plt.get_cmap('tab10')
    for idx, region in enumerate(regions):
        ax = axes.flatten()[idx]
        ax.set_xlim(0, mb)
        ax.set_ylim(0, y_max * 1.05)
        ax.tick_params(axis='x', labelbottom=True)
        ax.tick_params(axis='y', labelleft=True)
        times = durations.get(region, [])
        if times:
            ax.hist(times, bins=edges, color=cmap(idx), edgecolor='black', alpha=0.7)
        ax.set_title(f'{region} (n={len(times)})')
        ax.set_xlabel('Residence time (ns)')
        ax.set_ylabel('Count')
        ax.grid(True, linestyle='--', alpha=0.4)
    for extra_ax in axes.flatten()[len(regions):]:
        fig.delaxes(extra_ax)
    plt.suptitle(f'{run_id} residence times by region', fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fname = f"{run_id}_residence_time_side_by_side.png"
    outpath = os.path.join(outdir, fname)
    plt.savefig(outpath)
    plt.close()
    print(f'Saved side-by-side {outpath}')


def main():
    parser = argparse.ArgumentParser(
        description='Compute and plot residence time histograms by Pt region'
    )
    parser.add_argument(
        '--csv', required=True,
        help='Path to fragment metrics CSV (e.g. 0HPt_fragment_metrics.csv)'
    )
    parser.add_argument(
        '--outdir', default='residence_time_histograms',
        help='Directory to save histogram figures'
    )
    parser.add_argument(
        '--max-bin', type=float, default=None,
        help='Maximum residence time (ns) for histogram x-axis'
    )
    parser.add_argument(
        '--bins', type=int, default=50,
        help='Number of bins for histograms'
    )
    parser.add_argument(
        '--cutoff', type=float, default=2.5,
        help='Cutoff distance (Å) for surface residence'
    )
    parser.add_argument(
        '--min-duration', type=float, default=0.1,
        help='Minimum contiguous residence time (ns) to qualify'
    )
    parser.add_argument(
        '--overlay', action='store_true',
        help='Overlay histograms for all regions into one plot'
    )
    parser.add_argument(
        '--side-by-side', action='store_true',
        help='Plot four histograms side by side in one figure'
    )
    parser.add_argument(
        '--output-json', required=False, # Make output-json optional
        help='Path to the output JSON file for the summary results.'
    )
    args = parser.parse_args()

    # load data
    try:
        df = pd.read_csv(args.csv)
    except FileNotFoundError:
        print(f"Error: Input CSV file not found at {args.csv}")
        # If output-json is provided, write an error summary
        if args.output_json:
            error_summary = {"stage": "Residence Time Analysis", "error": f"Input CSV file not found at {args.csv}"}
            with open(args.output_json, 'w') as f:
                json.dump(error_summary, f, indent=4)
        sys.exit(1) # Exit with error code
    except Exception as e:
        print(f"Error loading CSV file: {e}")
        # If output-json is provided, write an error summary
        if args.output_json:
            error_summary = {"stage": "Residence Time Analysis", "error": f"Error loading CSV file: {e}"}
            with open(args.output_json, 'w') as f:
                json.dump(error_summary, f, indent=4)
        sys.exit(1) # Exit with error code


    # infer run_id
    run_id = df['run_id'].iloc[0] if 'run_id' in df.columns else os.path.splitext(os.path.basename(args.csv))[0]

    # compute sampling interval (time_ps in ps, convert to ns)
    try:
        dt_ps = compute_dt_ps(df)
        dt = dt_ps / 1000.0
        print(f'Determined sampling interval: {dt_ps:.3f} ps ({dt:.6f} ns)')
    except ValueError as e:
        print(f"Error determining sampling interval: {e}")
        # If output-json is provided, write an error summary
        if args.output_json:
            error_summary = {"stage": "Residence Time Analysis", "error": f"Error determining sampling interval: {e}"}
            with open(args.output_json, 'w') as f:
                json.dump(error_summary, f, indent=4)
        sys.exit(1) # Exit with error code
    except Exception as e:
        print(f"An unexpected error occurred while determining sampling interval: {e}")
        # If output-json is provided, write an error summary
        if args.output_json:
            error_summary = {"stage": "Residence Time Analysis", "error": f"An unexpected error occurred while determining sampling interval: {e}"}
            with open(args.output_json, 'w') as f:
                json.dump(error_summary, f, indent=4)
        sys.exit(1) # Exit with error code


    # Get unique non-NaN regions
    regions = sorted(df['nearest_pt_class'].dropna().unique())
    if not regions:
        print("No valid Pt regions found in the data.")
        # If output-json is provided, write a summary indicating no regions found
        if args.output_json:
            no_regions_summary = {"stage": "Residence Time Analysis", "message": "No valid Pt regions found in the data."}
            with open(args.output_json, 'w') as f:
                json.dump(no_regions_summary, f, indent=4)
        sys.exit(0) # Exit successfully as no analysis is possible


    # Compute residence durations (in ns) using the utility function
    all_durations = {reg: [] for reg in regions}
    try:
        for frag_id, grp in df.groupby('fragment_id'):
            frag_durations = find_residence_events_per_region(
                grp, dt, args.cutoff, args.min_duration, regions
            )
            for region, durations_list in frag_durations.items():
                all_durations[region].extend(durations_list)
    except Exception as e:
        print(f"Error computing residence times: {e}")
        # If output-json is provided, write an error summary
        if args.output_json:
            error_summary = {"stage": "Residence Time Analysis", "error": f"Error computing residence times: {e}"}
            with open(args.output_json, 'w') as f:
                json.dump(error_summary, f, indent=4)
        sys.exit(1) # Exit with error code


    # Construct the output directory based on run_id
    # Assuming output structure like outputs/<run_id>/residence_time_histograms
    base_output_dir = os.path.dirname(os.path.dirname(args.csv)) # Go up two levels from CSV path
    if os.path.basename(base_output_dir) != "outputs":
         # Fallback or adjust if CSV path is not in expected outputs/<run_id>/ structure
         base_output_dir = "outputs" # Default to 'outputs' if structure is unexpected
         run_output_dir = os.path.join(base_output_dir, run_id)
    else:
        run_output_dir = os.path.join(base_output_dir, run_id)


    plot_output_dir = os.path.join(run_output_dir, args.outdir)

    # plot and save histograms
    generated_plots = []
    try:
        if args.side_by_side:
            plot_side_by_side(all_durations, plot_output_dir, run_id, max_bin=args.max_bin, bins=args.bins)
            fname = f"{run_id}_residence_time_side_by_side.png"
            generated_plots.append(os.path.join(plot_output_dir, fname))
        elif args.overlay:
            plot_overlay(all_durations, plot_output_dir, run_id, max_bin=args.max_bin, bins=args.bins)
            fname = f"{run_id}_residence_time_overlay.png"
            generated_plots.append(os.path.join(plot_output_dir, fname))
        else:
            plot_histograms(all_durations, plot_output_dir, run_id, max_bin=args.max_bin, bins=args.bins)
            for region in regions:
                fname = f"{run_id}_{region}_residence_time_hist.png"
                generated_plots.append(os.path.join(plot_output_dir, fname))
    except Exception as e:
        print(f"Error generating plots: {e}")
        # If output-json is provided, write an error summary
        if args.output_json:
            error_summary = {"stage": "Residence Time Analysis", "error": f"Error generating plots: {e}"}
            with open(args.output_json, 'w') as f:
                json.dump(error_summary, f, indent=4)
        sys.exit(1) # Exit with error code


    # Prepare results for the comprehensive summary
    summary_results = {
        "stage": "Residence Time Analysis",
        "regions_analyzed": regions,
        "residence_event_summary": {},
        "generated_plots": generated_plots # Use the collected plot paths
    }

    for region, durations_list in all_durations.items():
        summary_results["residence_event_summary"][region] = {
            "event_count": len(durations_list),
            "average_duration_ns": float(np.mean(durations_list)) if durations_list else 0.0,
            "median_duration_ns": float(np.median(durations_list)) if durations_list else 0.0,
            "total_residence_time_ns": float(np.sum(durations_list)) if durations_list else 0.0
        }

    # Save summary results to JSON file if output-json is provided
    if args.output_json:
        try:
            with open(args.output_json, 'w') as f:
                json.dump(summary_results, f, indent=4)
            print(f"Exported summary: {args.output_json}")
        except Exception as e:
            print(f"Error saving summary JSON to {args.output_json}: {e}")
            # Continue without exiting, as the analysis itself was successful

    # Print summary results as JSON to stdout (still keep for orchestrator)
    print(json.dumps(summary_results))

    print("Residence time analysis complete.")

if __name__ == '__main__':
    main()