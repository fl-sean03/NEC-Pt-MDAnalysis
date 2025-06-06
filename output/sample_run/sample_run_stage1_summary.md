# Run Summary: sample_run

## Project Overview

This run is part of the MDAnalysis Script Regeneration project, which aims to rebuild Python scripts for analyzing molecular dynamics simulation data. The objective is to create a robust and maintainable pipeline for calculating key metrics and performing convergence checks and visualizations.

## System Information

- System time step: 1 fs
- Total number of frames: 924
- Number of NEC molecules (fragments): 999
- Number of Pt atoms: 561
- Box dimensions (Å): [73.06975555419922, 73.06975555419922, 73.06975555419922, 90.0, 90.0, 90.0]

## Run Parameters

- Sample interval (frames): 50
- Default contact cutoff (Å): 2.5
- Pt classification cutoff (Å): 3.0
- ΔG_D calculation temperature: 453 K (used in analyze_binding_ratios.py)
- Minimum event duration threshold for analysis: 0.1 ns (used in analysis scripts)
- Minimum event count threshold for facet inclusion in averages: N_events >= 10 (used in analyze_binding_ratios.py)

Parameters for analysis scripts (bins, max_bin) are configured within those scripts.

## Analysis Pipeline Stages

The analysis pipeline consists of several modular stages:

### Stage 1: Data Preparation (`nec_pt_fragment_metrics.py`)
This script (this run) reads raw simulation data (PSF, DCD) and extracts per-fragment, per-frame interaction metrics, including center-of-mass, minimum distance to the Pt surface (atom-atom based), contact fraction, tilt angle, radius of gyration, asphericity, principal moments, nearest Pt atom index, and nearest Pt atom class. It also classifies Pt atoms into regions based on coordination number.

### Stage 2: Residence Time Analysis (`analyze_residence_times.py`, `analyze_max_residence_times.py`)
These scripts read the fragment metrics CSV generated in Stage 1 and calculate/plot residence time distributions and maximum residence times. Residence events are defined consistently using utility functions.

### Stage 3: Binding Metrics Analysis (`analyze_binding_ratios.py`)
This script reads the fragment metrics CSV generated in Stage 1 and calculates/plots binding metrics such as dissociation constant (K_D), free energy of dissociation (ΔG_D), mean residence time (τ), and standard deviation of τ.

### Stage 4: Convergence Analysis (`analyze_convergence.py`)
This script reads the fragment metrics CSV generated in Stage 1 and performs/plots various convergence checks to assess if the simulation has reached a stable state.

### Stage 5: Advanced Visualization (`plot_residence_time_traces.py`)
This script reads the fragment metrics CSV generated in Stage 1 and generates advanced visualizations, such as spaghetti plots showing individual fragment traces over time, with features like simulation box boundaries.

## Data Verification

The `verify_min_dist_and_region.py` script is used to verify the consistency of the calculated fragment metrics, particularly the minimum distance and residence event definition, against the raw trajectory data for specific fragments.

## Output Files

- Fragment metrics CSV: `sample_run_fragment_metrics.csv`
- Pt classification JSON: `sample_run_pt_classification.json`
- Run Summary (this file): `sample_run_stage1_summary.md`
- Residence time histograms: `residence_time_histograms/` (generated by analyze_residence_times.py)
- Max residence time plots: `max_residence_plots/` (generated by analyze_max_residence_times.py)
- Binding metrics plots: `binding_metrics_plots/` (generated by analyze_binding_ratios.py)
- Convergence analysis plots: `convergence_plots/` (generated by analyze_convergence.py)
- Visualization plots: `visualization_plots/` (generated by plot_residence_time_traces.py)

Plot files within the histogram and plot directories are named based on the run ID and region/metric.

## Known Issues

A persistent `KeyError: 'min_distance'` has been observed in `analyze_convergence.py` during testing with sample data, despite the column being present in the DataFrame. This issue requires further investigation in Phase 7.

