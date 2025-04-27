#!/usr/bin/env python3
"""
Regenerated script to plot the maximum contiguous residence time per fragment for each Pt region.
"""
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

def compute_dt_ns(df):
    """
    Determine sampling interval in nanoseconds from the time_ps column.
    """
    times = np.sort(df['time_ps'].unique())
    if len(times) < 2:
        raise ValueError('Not enough time points to determine dt')
    dt_ps = np.median(np.diff(times))
    return float(dt_ps) / 1000.0

def compute_max_residence_times(df, dt_ns, cutoff):
    """
    For each fragment and each region, find the longest contiguous run where
    min_dist <= cutoff, then record its duration (in ns).
    Returns a dict: region -> list of max durations (ns) per fragment.
    """
    regions = sorted(df['nearest_pt_class'].dropna().unique())
    durations = {reg: [] for reg in regions}
    for frag_id, grp in df.groupby('fragment_id'):
        grp = grp[['time_ps', 'nearest_pt_class', 'min_dist']]
        grp = grp.sort_values('time_ps').reset_index(drop=True)
        n = len(grp)
        for region in regions:
            mask = ((grp['nearest_pt_class'] == region) & (grp['min_dist'] <= cutoff)).values
            max_run = 0
            current = 0
            for m in mask:
                if m:
                    current += 1
                else:
                    if current > max_run:
                        max_run = current
                    current = 0
            if current > max_run:
                max_run = current
            if max_run > 0:
                durations[region].append(max_run * dt_ns)
    return durations

def plot_side_by_side(durations, outdir, run_id, max_bin=None, bins=50):
    """
    Plot four histograms side by side for each region, using max residence times.
    """
    os.makedirs(outdir, exist_ok=True)
    regions = sorted(durations.keys())
    # collect all durations for common x-range
    all_times = [t for times in durations.values() for t in times]
    if not all_times:
        print('No residence data to plot')
        return
    mb = max_bin if max_bin is not None else max(all_times)
    edges = np.linspace(0, mb, bins + 1)
    # precompute histogram counts to unify y-axis range
    counts_per_region = {}
    for region in regions:
        times = durations.get(region, [])
        if times:
            counts, _ = np.histogram(times, bins=edges)
        else:
            counts = np.zeros(len(edges) - 1, dtype=int)
        counts_per_region[region] = counts
    y_max = max(counts.max() for counts in counts_per_region.values())

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=False, sharey=False) # Set sharey to False to allow explicit setting
    cmap = plt.get_cmap('tab10')
    for idx, region in enumerate(regions):
        ax = axes.flatten()[idx]
        # enforce common x- and y-limits
        ax.set_xlim(0, mb)
        ax.set_ylim(0, y_max * 1.05) # Set y-limit
        # ensure ticks and labels on all subplots
        ax.tick_params(axis='x', labelbottom=True)
        ax.tick_params(axis='y', labelleft=True)
        times = durations.get(region, [])
        if times:
            ax.hist(times, bins=edges, color=cmap(idx), edgecolor='black', alpha=0.7)
        ax.set_title(f'{region} (n={len(times)})')
        ax.set_xlabel('Max residence time (ns)')
        ax.set_ylabel('Count')
        ax.grid(True, linestyle='--', alpha=0.4)
    # remove unused axes
    for extra in axes.flatten()[len(regions):]:
        fig.delaxes(extra)
    plt.suptitle(f'{run_id} max residence times by region', fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    fname = f"{run_id}_max_residence_time_side_by_side.png"
    outpath = os.path.join(outdir, fname)
    plt.savefig(outpath)
    plt.close()
    print(f'Saved side-by-side max-residence {outpath}')


def main():
    parser = argparse.ArgumentParser(
        description='Plot max contiguous residence times per fragment per region.'
    )
    parser.add_argument('--csv', required=True, help='Fragment metrics CSV')
    parser.add_argument('--outdir', default='max_residence_plots', help='Output directory')
    parser.add_argument('--cutoff', type=float, default=2.5, help='Distance cutoff (Å)')
    parser.add_argument('--bins', type=int, default=50, help='Number of histogram bins')
    parser.add_argument('--max-bin', type=float, default=None, help='Max x-axis (ns)')
    args = parser.parse_args()

    # load data
    try:
        df = pd.read_csv(args.csv)
    except FileNotFoundError:
        print(f"Error: Input CSV file not found at {args.csv}")
        return
    except Exception as e:
        print(f"Error loading CSV file: {e}")
        return

    # infer run_id
    run_id = df['run_id'].iloc[0] if 'run_id' in df.columns else os.path.splitext(os.path.basename(args.csv))[0]

    # compute sampling interval (time_ps in ps, convert to ns)
    try:
        dt_ns = compute_dt_ns(df)
        print(f'Sampling interval: {dt_ns:.3f} ns')
    except ValueError as e:
        print(f"Error determining sampling interval: {e}")
        return
    except Exception as e:
        print(f"An unexpected error occurred while determining sampling interval: {e}")
        return

    # compute max residence durations (in ns)
    try:
        durations = compute_max_residence_times(df, dt_ns, args.cutoff)
    except Exception as e:
        print(f"Error computing maximum residence times: {e}")
        return

    # Construct the output directory based on run_id
    # Assuming output structure like outputs/<run_id>/max_residence_plots
    base_output_dir = os.path.dirname(os.path.dirname(args.csv)) # Go up two levels from CSV path
    if os.path.basename(base_output_dir) != "outputs":
         # Fallback or adjust if CSV path is not in expected outputs/<run_id>/ structure
         base_output_dir = "outputs" # Default to 'outputs' if structure is unexpected
         run_output_dir = os.path.join(base_output_dir, run_id)
    else:
        run_output_dir = os.path.join(base_output_dir, run_id)

    plot_output_dir = os.path.join(run_output_dir, args.outdir)

    # plot and save histograms
    try:
        plot_side_by_side(durations, plot_output_dir, run_id, max_bin=args.max_bin, bins=args.bins)
    except Exception as e:
        print(f"Error generating plots: {e}")
        return

    # Prepare results for the comprehensive summary
    summary_results = {
        "stage": "Maximum Residence Time Analysis",
        "regions_analyzed": sorted(durations.keys()),
        "max_residence_time_summary": {},
        "generated_plots": []
    }

    for region, durations_list in durations.items():
        summary_results["max_residence_time_summary"][region] = {
            "fragment_count_with_max_residence": len(durations_list),
            "average_max_duration_ns": float(np.mean(durations_list)) if durations_list else 0.0,
            "median_max_duration_ns": float(np.median(durations_list)) if durations_list else 0.0,
        }

    # Collect generated plot path
    # Assuming plot is saved with run_id in filename
    fname = f"{run_id}_max_residence_time_side_by_side.png"
    summary_results["generated_plots"].append(os.path.join(plot_output_dir, fname))

    # Print summary results as JSON to stdout
    print(json.dumps(summary_results))

    print("Maximum residence time analysis complete.")

if __name__ == '__main__':
    main()