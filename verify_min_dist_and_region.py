import MDAnalysis as mda
import pandas as pd
import json
import numpy as np
from MDAnalysis.lib import distances
from scipy.spatial import cKDTree # Need cKDTree for classify_pt_atoms
import sys
import traceback

# ----- Pt classification (copied from nec_pt_fragment_metrics.py) -----
def classify_pt_atoms(u, pt_cutoff):
    """
    Classify Pt atoms into regions (bulk, 100, 111, edge, vertex)
    based on coordination number using cKDTree with cutoff pt_cutoff.

    Returns:
      pt_class (dict): atom_index -> region string
      pts_by_region (dict): region -> AtomGroup of Pt atoms
    """
    pt_all = u.select_atoms("name PT*")
    coords = pt_all.positions
    tree = cKDTree(coords)
    # compute coordination numbers
    cn = []
    for pos in coords:
        neighbors = tree.query_ball_point(pos, r=pt_cutoff)
        cn.append(len(neighbors) - 1)
    cn = np.array(cn)
    # classify
    def _classify(c):
        if c >= 11:
            return 'bulk'
        elif c == 9:
            return '111'
        elif c == 8:
            return '100'
        elif 6 <= c < 8:
            return 'edge'
        else:
            return 'vertex'
    regions = np.array([_classify(c) for c in cn])
    # build mapping
    pt_class = {idx: reg for idx, reg in zip(pt_all.indices, regions)}
    # group atoms by region
    pts_by_region = {}
    for reg in ['bulk','100','111','edge','vertex']:
        sel = [i for i, r in pt_class.items() if r == reg]
        pts_by_region[reg] = u.atoms[sel]
    return pt_class, pts_by_region

# ----- Main verification logic -----
psf_path = 'data/0H/0HPt.psf'
dcd_path = 'data/0H/out_eq.dcd'
csv_path = 'outputs/com_metrics_run/com_fragment_metrics.csv' # Updated to use the new CSV
json_path = 'outputs/com_metrics_run/com_pt_classification.json' # Updated to use the new JSON
fragment_id_to_check = 20 # Fragment ID for "most unique regions"
sample_interval = 500 # Updated to match the CSV generation sample interval
distance_tolerance = 0.05 # Increased tolerance for comparing distances
pt_classification_cutoff = 3.0 # Cutoff used in nec_pt_fragment_metrics.py

try:
    # Load universe
    u = mda.Universe(psf_path, dcd_path)

    # Classify Pt atoms using the same method as nec_pt_fragment_metrics.py
    # Classify on the first frame (frame 0) for consistency
    u.trajectory[0]
    pt_classification, pts_by_region = classify_pt_atoms(u, pt_classification_cutoff)

    # Select surface Pt atoms based on the classification
    surface_regs = ['100','111','edge','vertex']
    all_surface_idxs = []
    for r in surface_regs:
        if r in pts_by_region:
            all_surface_idxs.extend(pts_by_region[r].indices.tolist())
    pts_surface = u.atoms[sorted(all_surface_idxs)]
    surface_pt_coords = pts_surface.positions


    # Load fragment metrics CSV
    df = pd.read_csv(csv_path)
    # Filter for the specific fragment ID and sort by frame
    df_fragment = df[df['fragment_id'] == fragment_id_to_check].sort_values('frame').reset_index(drop=True)

    if df_fragment.empty:
        print(f"Error: Fragment ID {fragment_id_to_check} not found in CSV.", file=sys.stderr)
        sys.exit(1)

    print(f"Verifying min_dist and nearest_pt_class for Fragment ID {fragment_id_to_check} every {sample_interval} frames:")

    # Use a counter to match frames with the filtered dataframe
    df_index = 0

    # Iterate through the trajectory frame by frame
    for ts in u.trajectory:
        frame = int(ts.frame)

        # Check if the current frame is one of the frames in the filtered dataframe
        if df_index < len(df_fragment) and int(df_fragment.iloc[df_index]['frame']) == frame:
            row = df_fragment.iloc[df_index]
            csv_min_dist = row['min_dist']
            csv_nearest_pt_class = row['nearest_pt_class']

            # Get fragment COM at this frame (current timestep)
            # Need to select the fragment correctly - assuming "not name PT*" selects the NEC fragments
            # and fragment_id is 1-based
            nec_fragments = u.select_atoms("not name PT*").fragments
            if fragment_id_to_check - 1 < len(nec_fragments):
                 fragment_to_check = nec_fragments[fragment_id_to_check - 1]
                 com = fragment_to_check.center_of_mass()
            else:
                 print(f"Error: Fragment ID {fragment_id_to_check} out of bounds.", file=sys.stderr)
                 df_index += 1 # Move to the next row to avoid infinite loop if fragment ID is bad
                 continue # Skip this frame/fragment


            # Calculate distance from fragment COM to surface Pt atoms
            if len(pts_surface) > 0:
                # Pass com as a (1, 3) numpy array
                dists_to_surface_pt = distances.distance_array(np.array([com]), pts_surface.positions, box=u.dimensions)[0]

                # Find minimum distance and nearest surface Pt atom index
                calculated_min_dist = np.min(dists_to_surface_pt)
                nearest_pt_atom_index_in_surface_selection = np.argmin(dists_to_surface_pt)
                nearest_pt_atom_index_in_universe = pts_surface.indices[nearest_pt_atom_index_in_surface_selection]

                # Get calculated nearest Pt atom classification
                calculated_nearest_pt_class = pt_classification.get(nearest_pt_atom_index_in_universe, 'Unknown')

                # Debug print statements for Frame 500 and Fragment ID 20
            else:
                # Handle case where there are no surface Pt atoms
                calculated_min_dist = np.inf
                nearest_pt_atom_index_in_universe = -1
                calculated_nearest_pt_class = 'none'


            # Compare with CSV values
            distance_match = abs(calculated_min_dist - csv_min_dist) < distance_tolerance
            region_match = str(calculated_nearest_pt_class) == str(csv_nearest_pt_class) # Convert to string for comparison with potential NaN in CSV

            status = "Match" if distance_match and region_match else "Mismatch"
            discrepancy_details = []
            if not distance_match:
                discrepancy_details.append(f"Distance Mismatch: Calculated={calculated_min_dist:.3f}, CSV={csv_min_dist:.3f}")
            if not region_match:
                 discrepancy_details.append(f"Region Mismatch: Calculated={calculated_nearest_pt_class}, CSV={csv_nearest_pt_class}")

            print(f"Frame {frame}: Status = {status}. COM = ({com[0]:.3f}, {com[1]:.3f}, {com[2]:.3f}) Å, Calculated MinDist = {calculated_min_dist:.3f}, Calculated NearestRegion = {calculated_nearest_pt_class}. CSV MinDist = {csv_min_dist:.3f}, CSV NearestRegion = {csv_nearest_pt_class}. {'(' + ', '.join(discrepancy_details) + ')' if discrepancy_details else ''}")

            df_index += 1 # Move to the next row in the filtered dataframe
        # If the current frame is not in the filtered dataframe, just continue to the next frame in the trajectory
        # No warning needed here as we are only checking frames present in the CSV

except FileNotFoundError as e:
   print(f"Error: File not found - {e}", file=sys.stderr)
   sys.exit(1)
except Exception as e:
   print(f"An error occurred: {e}", file=sys.stderr)
   traceback.print_exc(file=sys.stderr)
   sys.exit(1)