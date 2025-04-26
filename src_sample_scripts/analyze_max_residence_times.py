#!/usr/bin/env python3
"""
analyze_max_residence_times.py: Plot the maximum contiguous residence time per fragment for each Pt region.
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
    parser.add_argument('--outdir', default='max_residence_plots', help='Output directory') # This default will be overridden
    parser.add_argument('--cutoff', type=float, default=2.5, help='Distance cutoff (Å)')
    parser.add_argument('--bins', type=int, default=50, help='Number of histogram bins')
    parser.add_argument('--max-bin', type=float, default=None, help='Max x-axis (ns)')
    args = parser.parse_args()
    df = pd.read_csv(args.csv)
    run_id = df['run_id'].iloc[0] if 'run_id' in df.columns else os.path.splitext(os.path.basename(args.csv))[0]
    dt_ns = compute_dt_ns(df)
    print(f'Sampling interval: {dt_ns:.3f} ns')
    durations = compute_max_residence_times(df, dt_ns, args.cutoff)

    # Construct the output directory based on run_id
    base_output_dir = "2-1250Run/outputs"
    run_output_dir = os.path.join(base_output_dir, run_id)
    plot_output_dir = os.path.join(run_output_dir, args.outdir) # Use the original default or provided --outdir

    plot_side_by_side(durations, plot_output_dir, run_id, max_bin=args.max_bin, bins=args.bins)

if __name__ == '__main__':
    main()