#!/usr/bin/env python3
"""
analyze_residence_times.py: Compute residence time histograms for fragment contacts per Pt region.
"""
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import itertools

def compute_dt_ps(df):
    """
    Determine sampling interval (dt) in picoseconds from the time_ps column.
    """
    times = np.sort(df['time_ps'].unique())
    if len(times) < 2:
        raise ValueError('Not enough time points to determine dt')
    dt = np.median(np.diff(times))
    return float(dt)

def compute_residence_times(df, dt, cutoff, min_frames):
    """
    For each Pt region and each fragment, find contiguous residence events
    where the nearest atom is within cutoff (Å) for at least min_frames samples.
    Returns a dict: region -> list of durations (in ns).
    """
    regions = sorted(df['nearest_pt_class'].dropna().unique())
    durations = {reg: [] for reg in regions}
    # group by fragment
    for frag_id, grp in df.groupby('fragment_id'):
        grp = grp[['time_ps', 'nearest_pt_class', 'min_dist']]
        grp = grp.sort_values('time_ps').reset_index(drop=True)
        n = len(grp)
        # scan per region
        for region in regions:
            mask = ((grp['nearest_pt_class'] == region) &
                    (grp['min_dist'] <= cutoff)).values
            i = 0
            while i < n:
                if mask[i]:
                    # start of residence event
                    j = i + 1
                    while j < n and mask[j]:
                        j += 1
                    run_length = j - i
                    if run_length >= min_frames:
                        durations[region].append(run_length * dt)
                    i = j
                else:
                    i += 1
    return durations

def plot_histograms(durations, outdir, run_id, max_bin=None, bins=50):
    """
    Plot and save a histogram of residence times for each Pt region.
    """
    os.makedirs(outdir, exist_ok=True)
    for region, times in durations.items():
        if not times:
            print(f'No contact/residence events to plot for region: {region}')
            continue
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
    # optionally produce overlay plot? handled in main
    return

def plot_overlay(durations, outdir, run_id, max_bin=None, bins=50):
    """
    Plot and save an overlaid histogram of residence times for all Pt regions.
    """
    os.makedirs(outdir, exist_ok=True)
    # flatten all durations to determine x-range
    # collect all durations to set consistent bin edges
    all_times = [t for times in durations.values() for t in times]
    if not all_times:
        print('No residence events to plot in overlay')
        return
    mb = max_bin if max_bin is not None else max(all_times)
    edges = np.linspace(0, mb, bins + 1)
    plt.figure(figsize=(8, 5))
    # choose qualitative colors
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
    # legend outside
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
    # determine common x-axis limits and histogram bins
    all_times = [t for times in durations.values() for t in times]
    if not all_times:
        print('No residence events to plot side by side')
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
    # 2x2 grid; no shared axes so each subplot shows its ticks
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=False, sharey=False)
    cmap = plt.get_cmap('tab10')
    for idx, region in enumerate(regions):
        ax = axes.flatten()[idx]
        # enforce common x- and y-limits
        ax.set_xlim(0, mb)
        ax.set_ylim(0, y_max * 1.05)
        # ensure ticks and labels on all subplots
        ax.tick_params(axis='x', labelbottom=True)
        ax.tick_params(axis='y', labelleft=True)
        # plot histogram
        times = durations.get(region, [])
        if times:
            ax.hist(times, bins=edges, color=cmap(idx), edgecolor='black', alpha=0.7)
        ax.set_title(f'{region} (n={len(times)})')
        ax.set_xlabel('Residence time (ns)')
        ax.set_ylabel('Count')
        ax.grid(True, linestyle='--', alpha=0.4)
    # remove extra axes if any
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
        '--outdir', default='residence_time_histograms', # This default will be overridden
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
    args = parser.parse_args()
    # load data
    df = pd.read_csv(args.csv)
    # infer run_id
    run_id = df['run_id'].iloc[0] if 'run_id' in df.columns else os.path.splitext(os.path.basename(args.csv))[0]
    # compute sampling interval (time_ps in ps, convert to ns)
    dt_ps = compute_dt_ps(df)
    dt = dt_ps / 1000.0
    print(f'Determined sampling interval: {dt_ps:.3f} ps ({dt:.6f} ns)')
    # compute minimum frames from min-duration
    min_frames = int(math.ceil(args.min_duration / dt))
    print(f'Filtering for residence >= {args.min_duration} ns = {min_frames} frames (cutoff {args.cutoff} Å)')
    # compute residence durations (in ns)
    durations = compute_residence_times(df, dt, args.cutoff, min_frames)
    # Construct the output directory based on run_id
    base_output_dir = "2-1250Run/outputs"
    run_output_dir = os.path.join(base_output_dir, run_id)
    plot_output_dir = os.path.join(run_output_dir, args.outdir) # Use the original default or provided --outdir

    # plot and save histograms
    if args.side_by_side:
        # side-by-side panel of all regions
        plot_side_by_side(durations, plot_output_dir, run_id, max_bin=args.max_bin, bins=args.bins)
    elif args.overlay:
        # combined overlay of all regions
        plot_overlay(durations, plot_output_dir, run_id, max_bin=args.max_bin, bins=args.bins)
    else:
        plot_histograms(durations, plot_output_dir, run_id, max_bin=args.max_bin, bins=args.bins)

if __name__ == '__main__':
    main()