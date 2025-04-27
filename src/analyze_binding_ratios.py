#!/usr/bin/env python3
"""
Regenerated script to compute and plot dissociation constants per molecule and per Pt region.
"""
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from residence_analysis_utils import find_residence_events_per_region, count_residence_events_per_region # Import utility functions

def compute_dt_ns(df):
    """
    Determine sampling interval in nanoseconds from the time_ps column.
    """
    times = np.sort(df['time_ps'].unique())
    if len(times) < 2:
        raise ValueError('Not enough time points to determine dt')
    dt_ps = np.median(np.diff(times))
    return float(dt_ps) / 1000.0

def calculate_dissociation_constants(df, dt_ns, cutoff):
    """
    Calculate overall and per-facet dissociation constants for each fragment.
    Returns a dictionary: fragment_id -> { 'overall_constant': float, 'facet_constants': { region: float }, 'contacted_facets': [region] }
    """
    dissociation_data = {}
    regions = sorted(df['nearest_pt_class'].dropna().unique())
    all_fragments = df['fragment_id'].unique()

    for frag_id in all_fragments:
        frag_df = df[df['fragment_id'] == frag_id].sort_values('time_ps').reset_index(drop=True)
        total_frames = len(frag_df)
        total_sim_time = total_frames * dt_ns # in ns

        # Overall binding ratio
        on_surface_frames = frag_df[frag_df['min_dist'] <= cutoff]
        time_on_surface = len(on_surface_frames) * dt_ns
        time_off_surface = total_sim_time - time_on_surface

        # Corrected K_D calculation: total time off / total time on
        overall_constant = float('inf') if time_on_surface == 0 and time_off_surface > 0 else (time_off_surface / time_on_surface if time_on_surface > 0 else 0.0)

        # Per-facet dissociation constants and contacted facets
        facet_constants = {}
        contacted_facets = set()
        for region in regions:
            on_facet_frames = frag_df[(frag_df['nearest_pt_class'] == region) & (frag_df['min_dist'] <= cutoff)]
            time_on_facet = len(on_facet_frames) * dt_ns
            time_anywhere_else = total_sim_time - time_on_facet

            # Corrected K_D calculation: total time off / total time on for the facet
            facet_constant = float('inf') if time_on_facet == 0 and time_anywhere_else > 0 else (time_anywhere_else / time_on_facet if time_on_facet > 0 else 0.0)
            facet_constants[region] = facet_constant

            if time_on_facet > 0:
                contacted_facets.add(region)

        dissociation_data[frag_id] = {
            'overall_constant': overall_constant,
            'facet_constants': facet_constants,
            'contacted_facets': list(contacted_facets)
        }

    return dissociation_data, regions

# Removed the old count_residence_events function

def calculate_delta_g_d(k_d, temp_k=453.0):
    """
    Calculate Free Energy of Dissociation (ΔG_D) from K_D.
    R = 8.314 J/(mol·K) = 0.008314 kJ/(mol·K)
    ΔG_D = R * T * ln(K_D)
    """
    R = 0.008314 # kJ/(mol·K)
    # Handle K_D = 0 or inf
    if k_d <= 0:
        return float('inf') # Or some large negative number depending on convention
    if k_d == float('inf'):
        return float('inf') # Dissociation is infinitely favorable
    return R * temp_k * math.log(k_d)

def calculate_average_facet_metrics(df, dt_ns, cutoff, dissociation_data, regions, event_counts, min_event_threshold, temp_k):
    """
    Calculate the average dissociation constant, Delta G_D, mean residence time,
    and std dev of residence time per facet, filtering by minimum event count.
    Returns a dictionary: region -> { 'avg_k_d': float, 'avg_delta_g_d': float, 'mean_tau': float, 'std_tau': float, 'n_events': int }
    """
    avg_facet_metrics = {}
    min_tau_duration_ns = 0.1 # Use 0.1 ns minimum duration for tau calculation

    for region in regions:
        if region not in event_counts or event_counts[region] < min_event_threshold:
            print(f"Excluding region {region} due to insufficient events ({event_counts.get(region, 0)} < {min_event_threshold})")
            continue

        constants_for_facet = []
        all_durations_for_facet = [] # in ns

        for frag_id, data in dissociation_data.items():
            # Use the utility function to get residence durations for this fragment and region
            frag_df = df[df['fragment_id'] == frag_id]
            frag_durations = find_residence_events_per_region(
                frag_df, dt_ns, cutoff, min_tau_duration_ns, [region] # Pass only the current region
            )
            all_durations_for_facet.extend(frag_durations.get(region, []))


            if region in data['contacted_facets']:
                 if data['facet_constants'][region] != float('inf'):
                    constants_for_facet.append(data['facet_constants'][region])

        avg_k_d = np.mean(constants_for_facet) if constants_for_facet else 0.0
        avg_delta_g_d = calculate_delta_g_d(avg_k_d, temp_k)

        mean_tau = np.mean(all_durations_for_facet) if all_durations_for_facet else 0.0
        std_tau = np.std(all_durations_for_facet) if all_durations_for_facet else 0.0

        avg_facet_metrics[region] = {
            'avg_k_d': avg_k_d,
            'avg_delta_g_d': avg_delta_g_d,
            'mean_tau': mean_tau,
            'std_tau': std_tau,
            'n_events': event_counts.get(region, 0)
        }

    return avg_facet_metrics

# ----- Plotting functions -----
def plot_average_k_d(avg_facet_metrics, outdir, run_id):
    """
    Plot average K_D per facet.
    """
    os.makedirs(outdir, exist_ok=True)
    regions = list(avg_facet_metrics.keys())
    avg_k_d_values = [avg_facet_metrics[region]['avg_k_d'] for region in regions]

    plt.figure(figsize=(8, 5))
    plt.bar(regions, avg_k_d_values, color='skyblue')
    plt.xlabel('Pt Region')
    plt.ylabel('Average K_D')
    plt.title(f'{run_id} Average Dissociation Constant (K_D) per Pt Region')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    fname = f"{run_id}_avg_k_d_per_region.png"
    outpath = os.path.join(outdir, fname)
    plt.savefig(outpath)
    plt.close()
    print(f'Saved average K_D plot: {outpath}')

def plot_average_delta_g_d(avg_facet_metrics, outdir, run_id):
    """
    Plot average ΔG_D per facet.
    """
    os.makedirs(outdir, exist_ok=True)
    regions = list(avg_facet_metrics.keys())
    avg_delta_g_d_values = [avg_facet_metrics[region]['avg_delta_g_d'] for region in regions]

    plt.figure(figsize=(8, 5))
    plt.bar(regions, avg_delta_g_d_values, color='lightcoral')
    plt.xlabel('Pt Region')
    plt.ylabel('Average ΔG_D (kJ/mol)')
    plt.title(f'{run_id} Average Free Energy of Dissociation (ΔG_D) per Pt Region')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    fname = f"{run_id}_avg_delta_g_d_per_region.png"
    outpath = os.path.join(outdir, fname)
    plt.savefig(outpath)
    plt.close()
    print(f'Saved average ΔG_D plot: {outpath}')

def plot_mean_tau(avg_facet_metrics, outdir, run_id):
    """
    Plot mean residence time (τ) per facet.
    """
    os.makedirs(outdir, exist_ok=True)
    regions = list(avg_facet_metrics.keys())
    mean_tau_values = [avg_facet_metrics[region]['mean_tau'] for region in regions]
    std_tau_values = [avg_facet_metrics[region]['std_tau'] for region in regions]

    plt.figure(figsize=(8, 5))
    plt.bar(regions, mean_tau_values, yerr=std_tau_values, capsize=5, color='lightgreen')
    plt.xlabel('Pt Region')
    plt.ylabel('Mean Residence Time (τ) (ns)')
    plt.title(f'{run_id} Mean Residence Time (τ) per Pt Region')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    fname = f"{run_id}_mean_tau_per_region.png"
    outpath = os.path.join(outdir, fname)
    plt.savefig(outpath)
    plt.close()
    print(f'Saved mean τ plot: {outpath}')


def main():
    parser = argparse.ArgumentParser(
        description='Compute and plot dissociation constants per molecule and per Pt region.'
    )
    parser.add_argument('--csv', required=True, help='Fragment metrics CSV')
    parser.add_argument('--outdir', default='binding_metrics_plots', help='Output directory')
    parser.add_argument('--cutoff', type=float, default=2.5, help='Distance cutoff (Å)')
    parser.add_argument('--temperature', type=float, default=453.0, help='Temperature in Kelvin for ΔG_D calculation')
    parser.add_argument('--min-event-duration', type=float, default=0.1, help='Minimum contiguous residence time (ns) for event counting')
    parser.add_argument('--min-event-threshold', type=int, default=10, help='Minimum event count for facet inclusion in averages')
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

    # Get unique non-NaN regions
    regions = sorted(df['nearest_pt_class'].dropna().unique())
    if not regions:
        print("No valid Pt regions found in the data.")
        return

    # calculate dissociation constants per fragment
    try:
        dissociation_data, regions = calculate_dissociation_constants(df, dt_ns, args.cutoff)
    except Exception as e:
        print(f"Error calculating dissociation constants: {e}")
        return

    # count residence events per region using the utility function
    all_event_counts = {reg: 0 for reg in regions}
    try:
        for frag_id, grp in df.groupby('fragment_id'):
             frag_event_counts = count_residence_events_per_region(
                 grp, dt_ns, args.cutoff, args.min_event_duration, regions
             )
             for region, count in frag_event_counts.items():
                 all_event_counts[region] += count
    except Exception as e:
        print(f"Error counting residence events: {e}")
        return


    # calculate average facet metrics
    try:
        avg_facet_metrics = calculate_average_facet_metrics(df, dt_ns, args.cutoff, dissociation_data, regions, all_event_counts, args.min_event_threshold, args.temperature)
    except Exception as e:
        print(f"Error calculating average facet metrics: {e}")
        return

    # Construct the output directory based on run_id
    # Assuming output structure like outputs/<run_id>/binding_metrics_plots
    base_output_dir = os.path.dirname(os.path.dirname(args.csv)) # Go up two levels from CSV path
    if os.path.basename(base_output_dir) != "outputs":
         # Fallback or adjust if CSV path is not in expected outputs/<run_id>/ structure
         base_output_dir = "outputs" # Default to 'outputs' if structure is unexpected
         run_output_dir = os.path.join(base_output_dir, run_id)
    else:
        run_output_dir = os.path.join(base_output_dir, run_id)

    plot_output_dir = os.path.join(run_output_dir, args.outdir)
    os.makedirs(plot_output_dir, exist_ok=True)

    # plot and save metrics
    if avg_facet_metrics:
        try:
            plot_average_k_d(avg_facet_metrics, plot_output_dir, run_id)
            plot_average_delta_g_d(avg_facet_metrics, plot_output_dir, run_id)
            plot_mean_tau(avg_facet_metrics, plot_output_dir, run_id)
        except Exception as e:
            print(f"Error generating plots: {e}")
            return
    else:
        print("No average facet metrics to plot.")

    # Prepare results for the comprehensive summary
    summary_results = {
        "stage": "Binding Metrics Analysis",
        "average_facet_metrics": avg_facet_metrics,
        "generated_plots": []
    }

    # Collect generated plot paths
    if avg_facet_metrics:
        summary_results["generated_plots"].append(os.path.join(plot_output_dir, f"{run_id}_avg_k_d_per_region.png"))
        summary_results["generated_plots"].append(os.path.join(plot_output_dir, f"{run_id}_avg_delta_g_d_per_region.png"))
        summary_results["generated_plots"].append(os.path.join(plot_output_dir, f"{run_id}_mean_tau_per_region.png"))

    # Print summary results as JSON to stdout
    print(json.dumps(summary_results))

    print("Binding metrics analysis complete.")

if __name__ == '__main__':
    main()