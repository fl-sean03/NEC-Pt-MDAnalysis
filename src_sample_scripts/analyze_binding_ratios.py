#!/usr/bin/env python3
"""
analyze_binding_ratios.py: Compute and plot dissociation constants per molecule and per Pt region.
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
        # If only one time point, dt is undefined in this context, assume 0 or handle appropriately
        # For binding ratio, if only one frame, time on/off surface is trivial.
        # Let's return a small non-zero value or handle the single frame case later.
        # For now, raise error as in other scripts for consistency.
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

def count_residence_events(df, dt_ns, cutoff, min_duration_ns):
    """
    For each Pt region, count contiguous residence events where the nearest atom
    is within cutoff (Å) for at least min_duration_ns.
    Returns a dict: region -> event count.
    """
    regions = sorted(df['nearest_pt_class'].dropna().unique())
    event_counts = {reg: 0 for reg in regions}
    min_frames = int(math.ceil(min_duration_ns / dt_ns))

    for frag_id, grp in df.groupby('fragment_id'):
        grp = grp[['time_ps', 'nearest_pt_class', 'min_dist']]
        grp = grp.sort_values('time_ps').reset_index(drop=True)
        n = len(grp)

        for region in regions:
            mask = ((grp['nearest_pt_class'] == region) &
                    (grp['min_dist'] <= cutoff)).values
            i = 0
            while i < n:
                if mask[i]:
                    j = i + 1
                    while j < n and mask[j]:
                        j += 1
                    run_length = j - i
                    if run_length >= min_frames:
                        event_counts[region] += 1
                    i = j
                else:
                    i += 1
    return event_counts

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
    min_frames = int(math.ceil(0.1 / dt_ns)) # Use 0.1 ns minimum duration for tau calculation

    for region in regions:
        if region not in event_counts or event_counts[region] < min_event_threshold:
            print(f"Excluding region {region} due to insufficient events ({event_counts.get(region, 0)} < {min_event_threshold})")
            continue

        constants_for_facet = []
        durations_for_facet = [] # in ns

        for frag_id, data in dissociation_data.items():
            # Re-process fragment data to get residence durations for this region
            frag_df = df[df['fragment_id'] == frag_id].sort_values('time_ps').reset_index(drop=True)
            n = len(frag_df)
            mask = ((frag_df['nearest_pt_class'] == region) &
                    (frag_df['min_dist'] <= cutoff)).values
            i = 0
            while i < n:
                if mask[i]:
                    j = i + 1
                    while j < n and mask[j]:
                        j += 1
                    run_length = j - i
                    if run_length >= min_frames:
                        durations_for_facet.append(run_length * dt_ns)
                    i = j
                else:
                    i += 1

            if region in data['contacted_facets']:
                 if data['facet_constants'][region] != float('inf'):
                    constants_for_facet.append(data['facet_constants'][region])

        avg_k_d = np.mean(constants_for_facet) if constants_for_facet else 0.0
        avg_delta_g_d = calculate_delta_g_d(avg_k_d, temp_k)

        mean_tau = np.mean(durations_for_facet) if durations_for_facet else 0.0
        std_tau = np.std(durations_for_facet) if durations_for_facet else 0.0

        avg_facet_metrics[region] = {
            'avg_k_d': avg_k_d,
            'avg_delta_g_d': avg_delta_g_d,
            'mean_tau': mean_tau,
            'std_tau': std_tau,
            'n_events': event_counts.get(region, 0)
        }

    return avg_facet_metrics