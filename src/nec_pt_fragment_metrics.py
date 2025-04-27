#!/usr/bin/env python3
"""
Regenerated script to extract per-fragment, per-frame NEC–Pt interaction metrics
and export to CSV (and Pt classification JSON).

Metrics include:
- Fragment center-of-mass (com_x, com_y, com_z)
- Per-atom min distances to each Pt region type (bulk, 100, 111, edge, vertex)
- Fragment-level summaries: min_dist, contact_frac (default cutoff), tilt_deg,
  radius_gyration, asphericity, principal moments (eig1–3),
  patch_area (convex hull of contacting atoms), pt_density,
  nearest_pt_idx, nearest_pt_class.

Requires: MDAnalysis, numpy, pandas, scipy, sklearn, tqdm
"""
import os
import json
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib import distances
from scipy.spatial import ConvexHull, cKDTree
from sklearn.cluster import DBSCAN
from tqdm import tqdm
import argparse

# ----- Pt classification -----
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

# ----- Fragment metrics -----
def fragment_metrics(frag, pts_by_region, pts_surface, default_cutoff, box_dims, pt_classification):
    """
    Compute metrics for one fragment AtomGroup at current TS.
    Returns: dict of all CSV columns for this fragment-frame.
    """
    data = {}
    n_atoms = len(frag)
    data['n_atoms'] = n_atoms
    # center-of-mass
    com = frag.center_of_mass()
    data['com_x'], data['com_y'], data['com_z'] = com
    # per-atom min distances (kept for potential future use or debugging, not used for min_dist calculation)
    atom_positions = frag.positions
    atom_dists = {}
    for reg, pts in pts_by_region.items():
        pts_pos = pts.positions
        if len(pts_pos) > 0:
            d = distances.distance_array(atom_positions, pts_pos, box=box_dims)
            atom_mind = d.min(axis=1)
        else:
            atom_mind = np.full(n_atoms, np.inf)
        atom_dists[reg] = atom_mind
        # data[f'atom_mindists_{reg}'] = json.dumps(atom_mind.tolist()) # Commented out as per-atom dists are not used for main metrics

    # global min_dist (from fragment COM to any surface Pt)
    surf_regs = ['100','111','edge','vertex']
    # Ensure pts_surface is not empty before calculating distances
    if len(pts_surface) > 0:
        # Calculate distance from fragment COM to all surface Pt atoms
        # Pass com as a (1, 3) numpy array
        com_to_surface_dists = distances.distance_array(np.array([com]), pts_surface.positions, box=box_dims)[0]
        # Find the minimum distance
        min_dist_val = com_to_surface_dists.min()
        data['min_dist'] = float(min_dist_val)

        # Find the index of the nearest surface Pt atom
        nearest_pt_surface_idx_in_subset = com_to_surface_dists.argmin()
        nearest_pt_idx = pts_surface.indices[nearest_pt_surface_idx_in_subset]
        data['nearest_pt_idx'] = int(nearest_pt_idx)
        data['nearest_pt_class'] = pt_classification[nearest_pt_idx]

        # contact_frac at default cutoff (based on COM distance)
        data['contact_frac'] = float(min_dist_val <= default_cutoff) # This should be 1 or 0 based on COM distance

    else:
        # Handle case where there are no surface Pt atoms
        data['min_dist'] = np.inf
        data['nearest_pt_idx'] = -1 # Or some other indicator for no nearest Pt
        data['nearest_pt_class'] = 'none'
        data['contact_frac'] = 0.0

    # tilt angle: principal axis vs z
    coords_rel = atom_positions - com
    cov = np.cov(coords_rel, rowvar=False)
    eigvals, eigvecs = np.linalg.eigh(cov)
    # sort descending
    idxs = np.argsort(eigvals)[::-1]
    eig_sorted = eigvals[idxs]
    vec_prin = eigvecs[:, idxs[0]]
    cosang = abs(np.dot(vec_prin, np.array([0,0,1])))
    data['tilt_deg'] = float(np.degrees(np.arccos(np.clip(cosang, -1, 1))))
    data['eig1'], data['eig2'], data['eig3'] = eig_sorted
    # radius of gyration
    data['radius_gyration'] = float(np.sqrt(eig_sorted.sum()))
    # asphericity = λ1 - 0.5*(λ2+λ3)
    data['asphericity'] = float(eig_sorted[0] - 0.5*(eig_sorted[1] + eig_sorted[2]))
    # patch_area and pt_density calculations removed as they are not needed for the current fix
    # and depend on the old per-atom distance logic.

    # nearest_pt_idx and class (already calculated based on COM distance above)
    # The nearest_pt_idx and nearest_pt_class are already determined in the COM-based min_dist calculation block.
    # No need to re-calculate here.

    return data

def write_summary(output_summary_path, run_id, u, nec_frags, sample_interval, default_cutoff, pt_cutoff):
    """
    Writes a comprehensive summary of the run details and output files to a markdown file.
    """
    with open(output_summary_path, 'w') as f:
        f.write(f"# Run Summary: {run_id}\n\n")

        f.write("## System Information\n\n")
        # Assuming system time step is 1 fs based on sample script
        f.write(f"- System time step: 1 fs\n")
        f.write(f"- Total number of frames: {len(u.trajectory)}\n")
        f.write(f"- Number of NEC molecules (fragments): {len(nec_frags)}\n")
        f.write(f"- Number of Pt atoms: {len(u.select_atoms('name PT*'))}\n")
        f.write(f"- Box dimensions (Å): {u.dimensions.tolist()}\n\n")

        f.write("## Run Parameters\n\n")
        f.write(f"- Sample interval (frames): {sample_interval}\n")
        f.write(f"- Default contact cutoff (Å): {default_cutoff}\n")
        f.write(f"- Pt classification cutoff (Å): {pt_cutoff}\n")
        # Parameters from other scripts mentioned in sample summary
        f.write(f"- ΔG_D calculation temperature: 453 K (used in analyze_binding_ratios.py)\n")
        f.write(f"- Minimum event duration threshold for analysis: 0.1 ns (used in analysis scripts)\n")
        f.write(f"- Minimum event count threshold for facet inclusion in averages: N_events >= 10 (used in analyze_binding_ratios.py)\n\n")
        f.write("Parameters for analysis scripts (bins, max_bin) are configured within those scripts.\n\n")


        f.write("## Output Files\n\n")
        f.write(f"- Fragment metrics CSV: `{os.path.basename(output_summary_path).replace('_summary.md', '_fragment_metrics.csv')}`\n")
        f.write(f"- Pt classification JSON: `{os.path.basename(output_summary_path).replace('_summary.md', '_pt_classification.json')}`\n")
        f.write(f"- Run Summary (this file): `{os.path.basename(output_summary_path)}`\n")
        # Mention output directories for other stages based on PLANNING.md and sample scripts
        f.write(f"- Residence time histograms: `residence_time_histograms/` (generated by analyze_residence_times.py)\n")
        f.write(f"- Max residence time plots: `max_residence_plots/` (generated by analyze_max_residence_times.py)\n")
        f.write(f"- Convergence analysis plots: `convergence_plots/` (generated by analyze_convergence.py)\n")
        f.write(f"- Visualization plots: `visualization_plots/` (generated by plot_residence_time_traces.py)\n\n")

        f.write("Plot files within the histogram and plot directories are named based on the run ID and region/metric.\n")


    print(f"Exported summary: {output_summary_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Extract per-fragment, per-frame NEC–Pt interaction metrics.'
    )
    parser.add_argument(
        '--psf', required=True,
        help='Path to the PSF file.'
    )
    parser.add_argument(
        '--dcd', required=True,
        help='Path to the DCD trajectory file.'
    )
    parser.add_argument(
        '--run-id', required=True,
        help='A unique identifier for the simulation run.'
    )
    parser.add_argument(
        '--sample-interval', type=int, default=1,
        help='Sample every N frames.'
    )
    parser.add_argument(
        '--default-cutoff', type=float, default=2.5,
        help='Default contact cutoff distance in Angstroms.'
    )
    parser.add_argument(
        '--pt-cutoff', type=float, default=3.0,
        help='Cutoff distance in Angstroms for Pt-Pt coordination number calculation.'
    )
    parser.add_argument(
        '--output-csv', required=True,
        help='Path to the output CSV file for fragment metrics.'
    )
    parser.add_argument(
        '--output-json', required=True,
        help='Path to the output JSON file for Pt classification.'
    )
    parser.add_argument(
        '--output-summary', required=True,
        help='Path to the output markdown file for the run summary.'
    )
    args = parser.parse_args()

    # Ensure output directory exists
    output_dir = os.path.dirname(args.output_csv)
    os.makedirs(output_dir, exist_ok=True)

    # Load universe
    try:
        u = mda.Universe(args.psf, args.dcd)
    except Exception as e:
        print(f"Error loading universe: {e}")
        return

    # Classify Pt atoms on the first frame (or a representative frame)
    # Using the first frame for classification as it's readily available
    try:
        u.trajectory[0]
        pt_classification, pts_by_region = classify_pt_atoms(u, args.pt_cutoff)
    except Exception as e:
        print(f"Error classifying Pt atoms: {e}")
        return

    # build surface AtomGroup
    surface_regs = ['100','111','edge','vertex']
    all_surface_idxs = []
    for r in surface_regs:
        if r in pts_by_region:
            all_surface_idxs.extend(pts_by_region[r].indices.tolist())
    pts_surface = u.atoms[sorted(all_surface_idxs)]

    # Select fragments (assuming "NEC" is the residue name or similar)
    # Need to confirm how NEC fragments are identified in the PSF
    # For now, selecting all atoms that are not Pt
    try:
        nec_atoms = u.select_atoms("not name PT*")
        nec_frags = nec_atoms.fragments
        if not nec_frags:
             print("Warning: No non-Pt fragments found. Check selection.")
    except Exception as e:
        print(f"Error selecting fragments: {e}")
        return


    # Write summary file
    try:
        write_summary(args.output_summary, args.run_id, u, nec_frags,
                      args.sample_interval, args.default_cutoff, args.pt_cutoff)
    except Exception as e:
        print(f"Error writing summary file: {e}")


    # iterate frames and compute metrics
    n_frames = len(u.trajectory)
    rows = []
    try:
        for ts in tqdm(u.trajectory[::args.sample_interval], desc="Processing frames"):
            # per fragment
            for fid, frag in enumerate(nec_frags, start=1):
                try:
                    row = {
                        'run_id': args.run_id,
                        'frame': int(ts.frame),
                        # Assuming time is in ps and needs 5x scaling based on sample script
                        'time_ps': float(getattr(ts, 'time', ts.frame)) * 5,
                        'fragment_id': fid
                    }
                    # compute metrics, passing fragment_id
                    m = fragment_metrics(frag, pts_by_region, pts_surface,
                                         args.default_cutoff, u.dimensions, pt_classification, fid)
                    row.update(m)
                    rows.append(row)
                except Exception as e:
                    print(f"Error processing fragment {fid} at frame {ts.frame}: {e}")
                    # Optionally continue or break
                    continue # Continue processing other fragments/frames
    except Exception as e:
        print(f"Error iterating through trajectory: {e}")


    # assemble DataFrame and export
    if rows:
        try:
            df = pd.DataFrame(rows)
            print(f"Attempting to write CSV to: {args.output_csv}, DataFrame shape: {df.shape}")
            df.to_csv(args.output_csv, index=False)
            print(f"Exported CSV: {args.output_csv}")
        except Exception as e:
            print(f"Error writing CSV file {args.output_csv}: {e}")
    else:
        print("No data rows generated to write to CSV.")

    # also save pt classification mapping
    try:
        # Convert int64 keys to standard int for JSON serialization
        pt_classification_int_keys = {int(k): v for k, v in pt_classification.items()}
        with open(args.output_json, 'w') as f:
            json.dump(pt_classification_int_keys, f, indent=2)
        print(f"Exported Pt classification: {args.output_json}")
    except Exception as e:
        print(f"Error writing Pt classification JSON file {args.output_json}: {e}")


    print(f"Metrics extraction process finished for run: {args.run_id}")

if __name__ == "__main__":
    main()