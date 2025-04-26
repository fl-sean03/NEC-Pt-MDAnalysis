#!/usr/bin/env python3
"""
Comprehensive script to extract per-fragment, per-frame NEC–Pt interaction metrics
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
def fragment_metrics(frag, pts_by_region, pts_surface, default_cutoff, box_dims):
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
    # per-atom min distances
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
        data[f'atom_mindists_{reg}'] = json.dumps(atom_mind.tolist())
    # global min_dist (to any surface Pt)
    surf_regs = ['100','111','edge','vertex']
    all_surface_mind = np.vstack([atom_dists[r] for r in surf_regs])
    global_atom_mind = all_surface_mind.min(axis=0)
    data['min_dist'] = float(global_atom_mind.min())
    # contact_frac at default cutoff
    data['contact_frac'] = float((global_atom_mind <= default_cutoff).sum() / n_atoms)
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
    # set up contacting atoms mask
    mask = global_atom_mind <= default_cutoff
    pts_contact = atom_positions[mask]
    # patch_area via convex hull
    if len(pts_contact) >= 3:
        xy = pts_contact[:,:2]
        hull = ConvexHull(xy)
        data['patch_area'] = float(hull.volume)
        cen = xy[hull.vertices].mean(axis=0)
        data['patch_centroid_x'], data['patch_centroid_y'] = cen
    else:
        data['patch_area'] = 0.0
        data['patch_centroid_x'] = np.nan
        data['patch_centroid_y'] = np.nan
    # pt_density: count surface Pt within 5Å of fragment COM
    tree = cKDTree(pts_surface.positions)
    data['pt_density'] = int(len(tree.query_ball_point(com, r=5.0)))
    # nearest_pt_idx and class
    # compute dists from each atom to surface Pt, find global min
    d_all = distances.distance_array(atom_positions, pts_surface.positions, box=box_dims)
    aidx, pidx = np.unravel_index(d_all.argmin(), d_all.shape)
    nearest_pt_idx = pts_surface.indices[pidx]
    data['nearest_pt_idx'] = int(nearest_pt_idx)
    data['nearest_pt_class'] = pt_classification[nearest_pt_idx]
    return data

def write_summary(output_dir, run_id, u, nec_frags, sample_interval, default_cutoff, pt_cutoff):
    """
    Writes a comprehensive summary of the run details and output files to a markdown file.
    """
    summary_file = os.path.join(output_dir, f"{run_id}_summary.md")
    with open(summary_file, 'w') as f:
        f.write(f"# Run Summary: {run_id}\n\n")

        f.write("## System Information\n\n")
        f.write(f"- System time step: 1 fs\n")
        f.write(f"- Total number of frames: {len(u.trajectory)}\n")
        f.write(f"- Number of NEC molecules (fragments): {len(nec_frags)}\n")
        f.write(f"- Number of Pt atoms: {len(u.select_atoms('name PT*'))}\n")
        f.write(f"- Box dimensions (Å): {u.dimensions.tolist()}\n\n")

        f.write("## Run Parameters\n\n")
        f.write(f"- Sample interval (frames): {sample_interval}\n")
        f.write(f"- Default contact cutoff (Å): {default_cutoff}\n")
        f.write(f"- Pt classification cutoff (Å): {pt_cutoff}\n")
        f.write(f"- ΔG_D calculation temperature: 453 K\n")
        f.write(f"- Mean residence time (τ) and standard deviation are calculated from event durations.\n")
        f.write(f"- Minimum event duration threshold for analysis: 0.1 ns\n")
        f.write(f"- Minimum event count threshold for facet inclusion in averages: N_events >= 10\n\n")
        # Note: Bins and max_bin are parameters for the analysis scripts, not this one.
        # We can mention the defaults or note they are set in analysis scripts.
        f.write("Parameters for analysis scripts (bins, max_bin) are configured within those scripts.\n\n")


        f.write("## Output Files\n\n")
        f.write(f"- Fragment metrics CSV: `{os.path.join(output_dir, f'{run_id}_fragment_metrics.csv')}`\n")
        f.write(f"- Pt classification JSON: `{os.path.join(output_dir, f'{run_id}_pt_classification.json')}`\n")
        f.write(f"- Run Summary (this file): `{summary_file}`\n")
        f.write(f"- Residence time histograms: `{os.path.join(output_dir, 'residence_time_histograms/')}`\n")
        f.write(f"- Max residence time plots: `{os.path.join(output_dir, 'max_residence_plots/')}`\n")
        f.write(f"- Convergence analysis plots (from `analyze_convergence.py`):\n")
        f.write(f"    - Unique-facet coverage: Tracks the number of unique facets visited over time.\n")
        f.write(f"    - Radius of Gyration (Rg): Shows the evolution of fragment size/shape.\n")
        f.write(f"    - Mean Squared Displacement (MSD): Indicates fragment mobility.\n")
        f.write(f"    - Event-count growth: Monitors the accumulation of binding events over time.\n")
        f.write(f"- Visualization plots:\n")
        f.write(f"    - Mean τ comparison: Compares mean residence times across different facets.\n")
        f.write(f"    - Spaghetti plot: Visualizes individual fragment residence time traces.\n\n")

        f.write("Plot files within the histogram directories are named based on the run ID and region.\n")


    print(f"Exported summary: {summary_file}")


# ----- Main workflow -----
def main(data_dir, psf_file, dcd_file,
         run_id, sample_interval, default_cutoff,
         pt_cutoff, output_csv):
    # Load
    u = mda.Universe(psf_file, dcd_file)
    # Classify Pt atoms on middle frame
    mid = len(u.trajectory) // 2
    u.trajectory[mid]
    global pt_classification, pts_by_region
    pt_classification, pts_by_region = classify_pt_atoms(u, pt_cutoff)
    # build surface AtomGroup
    surface_regs = ['100','111','edge','vertex']
    all_surface_idxs = []
    for r in surface_regs:
        all_surface_idxs.extend(pts_by_region[r].indices.tolist())
    pts_surface = u.atoms[sorted(all_surface_idxs)]

    # Fragments
    nec_atoms     = u.select_atoms("not name PT*")
    nec_frags     = nec_atoms.fragments

    # Write summary file
    # Write summary file
    write_summary(os.path.dirname(output_csv), run_id, u, nec_frags,
                  sample_interval, default_cutoff, pt_cutoff)

    # iterate frames
    n_frames = len(u.trajectory)
    rows = []
    for frame in tqdm(range(0, n_frames, sample_interval), desc="Sampling frames"):
        ts = u.trajectory[frame]
        # per fragment
        for fid, frag in enumerate(nec_frags, start=1):
            row = {
                'run_id': run_id,
                'frame': int(ts.frame),
                'time_ps': float(getattr(ts, 'time', ts.frame)) * 5, # Apply 5x scaling
                'fragment_id': fid
            }
            # compute metrics
            m = fragment_metrics(frag, pts_by_region, pts_surface,
                                 default_cutoff, u.dimensions)
            row.update(m)
            rows.append(row)

    # assemble DataFrame
    df = pd.DataFrame(rows)
    print(f"Attempting to write CSV to: {output_csv}, DataFrame shape: {df.shape}")
    try:
        df.to_csv(output_csv, index=False)
    except Exception as e:
        print(f"Error writing CSV file {output_csv}: {e}")
    # also save pt classification mapping
    cls_file = os.path.join(os.path.dirname(output_csv), f"{run_id}_pt_classification.json")
    # Convert int64 keys to standard int for JSON serialization
    pt_classification_int_keys = {int(k): v for k, v in pt_classification.items()}
    with open(cls_file, 'w') as f:
        json.dump(pt_classification_int_keys, f, indent=2)
    print(f"Exported CSV: {output_csv}")
    print(f"Exported Pt classification: {cls_file}")

# ----- User inputs -----
if __name__ == "__main__":
    # Adjust these paths and parameters
    data_dir        = "./data/0H"
    psf_file        = os.path.join(data_dir, "0HPt.psf")
    dcd_file        = os.path.join(data_dir, "out_eq.dcd")
    run_id          = "0H_full_run"
    sample_interval = 100     # sample every N frames (was 1000)
    default_cutoff  = 2.5     # Å, default contact
    pt_cutoff       = 3.0     # Å, for Pt-Pt coordination
    # write output CSV to current working directory
    output_dir      = os.path.join("2-1250Run/outputs", run_id)
    print(f"Attempting to create directory: {output_dir}")
    try:
        os.makedirs(output_dir, exist_ok=True)
    except Exception as e:
        print(f"Error creating output directory {output_dir}: {e}")
    output_csv      = os.path.join(output_dir, f"{run_id}_fragment_metrics.csv")

    main(data_dir, psf_file, dcd_file,
         run_id, sample_interval, default_cutoff,
         pt_cutoff, output_csv)
