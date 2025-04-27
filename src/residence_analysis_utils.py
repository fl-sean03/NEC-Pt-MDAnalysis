import numpy as np
import pandas as pd
import math

def _find_contiguous_blocks(mask, min_length):
    """
    Helper function to find contiguous blocks of True values in a boolean mask
    that are at least min_length long.
    Returns a list of (start_index, end_index) tuples for each block.
    """
    blocks = []
    i = 0
    n = len(mask)
    while i < n:
        if mask[i]:
            # Start of a potential block
            j = i
            while j < n and mask[j]:
                j += 1
            # Block found from i to j-1
            run_length = j - i
            if run_length >= min_length:
                blocks.append((i, j - 1))
            i = j # Continue search after this block
        else:
            i += 1 # Continue search after this False value
    return blocks

def find_residence_events_per_region(df_fragment: pd.DataFrame, dt_ns: float, cutoff: float, min_duration_ns: float, regions: list):
    """
    For a single fragment, find contiguous residence events per Pt region
    where the fragment is within the distance cutoff AND is nearest to a specific
    non-NaN Pt region, for at least a minimum duration.

    Args:
        df_fragment (pd.DataFrame): DataFrame containing data for a single fragment,
                                    must include 'time_ps', 'nearest_pt_class', 'min_dist'.
        dt_ns (float): Sampling interval in nanoseconds.
        cutoff (float): Distance cutoff in Angstroms for surface residence.
        min_duration_ns (float): Minimum contiguous residence time in nanoseconds.
        regions (list): List of unique non-NaN Pt regions.

    Returns:
        dict: region -> list of durations (in ns) for residence events in that region.
    """
    durations = {reg: [] for reg in regions}
    min_frames = int(math.ceil(min_duration_ns / dt_ns))

    df_fragment = df_fragment.sort_values('time_ps').reset_index(drop=True)
    n = len(df_fragment)

    for region in regions:
        # Mask for frames within cutoff AND nearest to the specific region
        mask = ((df_fragment['nearest_pt_class'] == region) &
                (df_fragment['min_dist'] <= cutoff)).values

        # Find contiguous blocks meeting the minimum duration
        blocks = _find_contiguous_blocks(mask, min_frames)

        # Calculate durations for each block
        for start_idx, end_idx in blocks:
            run_length = end_idx - start_idx + 1
            durations[region].append(run_length * dt_ns)

    return durations

def count_residence_events_per_region(df_fragment: pd.DataFrame, dt_ns: float, cutoff: float, min_duration_ns: float, regions: list):
    """
    For a single fragment, count contiguous residence events per Pt region
    where the fragment is within the distance cutoff AND is nearest to a specific
    non-NaN Pt region, for at least a minimum duration.

    Args:
        df_fragment (pd.DataFrame): DataFrame containing data for a single fragment,
                                    must include 'time_ps', 'nearest_pt_class', 'min_dist'.
        dt_ns (float): Sampling interval in nanoseconds.
        cutoff (float): Distance cutoff in Angstroms for surface residence.
        min_duration_ns (float): Minimum contiguous residence time in nanoseconds.
        regions (list): List of unique non-NaN Pt regions.

    Returns:
        dict: region -> event count for residence events in that region.
    """
    event_counts = {reg: 0 for reg in regions}
    min_frames = int(math.ceil(min_duration_ns / dt_ns))

    df_fragment = df_fragment.sort_values('time_ps').reset_index(drop=True)
    n = len(df_fragment)

    for region in regions:
        # Mask for frames within cutoff AND nearest to the specific region
        mask = ((df_fragment['nearest_pt_class'] == region) &
                (df_fragment['min_dist'] <= cutoff)).values

        # Find contiguous blocks meeting the minimum duration
        blocks = _find_contiguous_blocks(mask, min_frames)

        event_counts[region] = len(blocks)

    return event_counts


def count_total_residence_events(df_fragment: pd.DataFrame, dt_ns: float, cutoff: float, min_duration_ns: float):
    """
    For a single fragment, count contiguous residence events where the fragment
    is within the distance cutoff AND is nearest to *any* non-NaN Pt region,
    for at least a minimum duration.

    Args:
        df_fragment (pd.DataFrame): DataFrame containing data for a single fragment,
                                    must include 'time_ps', 'nearest_pt_class', 'min_dist'.
        dt_ns (float): Sampling interval in nanoseconds.
        cutoff (float): Distance cutoff in Angstroms for surface residence.
        min_duration_ns (float): Minimum contiguous residence time in nanoseconds.

    Returns:
        int: The total number of residence events for the fragment.
    """
    min_frames = int(math.ceil(min_duration_ns / dt_ns))

    df_fragment = df_fragment.sort_values('time_ps').reset_index(drop=True)
    n = len(df_fragment)

    # Mask for frames within cutoff AND nearest to ANY non-NaN region
    mask = ((df_fragment['min_dist'] <= cutoff) &
            (df_fragment['nearest_pt_class'].notna())).values

    # Find contiguous blocks meeting the minimum duration
    blocks = _find_contiguous_blocks(mask, min_frames)

    return len(blocks)