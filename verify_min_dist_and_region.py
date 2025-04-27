import MDAnalysis as mda
import pandas as pd
import json
import numpy as np
from MDAnalysis.lib import distances
import sys
import traceback

psf_path = 'data/0H/0HPt.psf'
dcd_path = 'data/0H/out_eq.dcd'
csv_path = 'outputs/run_100frames/fragment_metrics.csv'
json_path = 'outputs/run_100frames/pt_classification.json'
fragment_id_to_check = 20 # Fragment ID for "most unique regions"
sample_interval = 100
distance_tolerance = 0.01 # Tolerance for comparing distances

try:
    # Load universe
    u = mda.Universe(psf_path, dcd_path)
    pt_atoms = u.select_atoms("name PT*")
    pt_coords = pt_atoms.positions

    # Load Pt classification
    with open(json_path, 'r') as f:
        pt_classification = json.load(f)
    pt_classification = {int(k): v for k, v in pt_classification.items()}

    # Load fragment metrics CSV
    df = pd.read_csv(csv_path)
    df_fragment = df[df['fragment_id'] == fragment_id_to_check].sort_values('frame').reset_index(drop=True)

    if df_fragment.empty:
        print(f"Error: Fragment ID {fragment_id_to_check} not found in CSV.", file=sys.stderr)
        sys.exit(1)

    print(f"Verifying min_dist and nearest_pt_class for Fragment ID {fragment_id_to_check} every {sample_interval} frames:")

    for index, row in df_fragment.iterrows():
        frame = int(row['frame'])
        csv_min_dist = row['min_dist']
        csv_nearest_pt_class = row['nearest_pt_class']

        # Set the trajectory to the current frame
        u.trajectory[frame]

        # Get fragment COM at this frame
        fragment_to_check = u.select_atoms("not name PT*").fragments[fragment_id_to_check - 1]
        com = fragment_to_check.center_of_mass()

        # Calculate distance from fragment COM to all Pt atoms
        dists_to_pt = distances.distance_array(com.reshape(1, 3), pt_coords, box=u.dimensions)[0]

        # Find minimum distance and nearest Pt atom index
        calculated_min_dist = np.min(dists_to_pt)
        nearest_pt_atom_index_in_selection = np.argmin(dists_to_pt)
        nearest_pt_atom_index_in_universe = pt_atoms.indices[nearest_pt_atom_index_in_selection]

        # Get calculated nearest Pt atom classification
        calculated_nearest_pt_class = pt_classification.get(nearest_pt_atom_index_in_universe, 'Unknown')

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


except FileNotFoundError as e:
    print(f"Error: File not found - {e}", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"An error occurred: {e}", file=sys.stderr)
    traceback.print_exc(file=sys.stderr)
    sys.exit(1)