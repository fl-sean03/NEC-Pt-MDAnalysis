# TASK.md: MDAnalysis Script Regeneration Tasks

## Overview

This document outlines the detailed tasks for regenerating the MDAnalysis data analysis scripts, following the plan defined in PLANNING.md.

## Phases and Task Breakdown: Building the Analysis Pipeline

### Phase 1: Planning & Setup

- [x] Create PLANNING.md
- [x] Create TASK.md (this document)
- [x] Set up project directory structure (e.g., ensure `src/` and `outputs/` exist).
- [x] Review provided sample scripts in `src_sample_scripts/` for reference and ideas (remembering to generate new code originally).

### Phase 2: Pipeline Stage 1 - Data Preparation (`nec_pt_fragment_metrics.py`)

- [ ] Regenerate `nec_pt_fragment_metrics.py`:
    - [ ] Define clear inputs (PSF, DCD file paths, parameters like cutoffs).
    - [ ] Implement data loading using MDAnalysis.
    - [ ] Implement Pt atom classification.
    - [ ] Implement per-fragment, per-frame metric calculations (COM, min distances, Rg, etc.).
    - [ ] Implement robust error handling for file loading and calculations.
    - [ ] Define clear outputs (fragment metrics CSV path, Pt classification JSON path, summary markdown path).
    - [ ] Implement CSV, JSON, and basic summary markdown output.
    - [ ] Add command-line arguments for all inputs and outputs.

### Phase 3: Pipeline Stage 2 - Residence Time Analysis (`analyze_residence_times.py`, `analyze_max_residence_times.py`)

- [ ] Regenerate `analyze_residence_times.py`:
    - [ ] Define clear inputs (fragment metrics CSV path, parameters like cutoff, min duration).
    - [ ] Implement data loading (fragment metrics CSV).
    - [ ] Implement sampling interval calculation.
    - [ ] Implement robust residence event identification with minimum duration threshold (0.1 ns).
    - [ ] Implement residence time histogram plotting.
    - [ ] Define clear outputs (histogram plot file paths).
    - [ ] Add command-line arguments for all inputs and outputs.
- [ ] Regenerate `analyze_max_residence_times.py`:
    - [ ] Define clear inputs (fragment metrics CSV path, parameters like cutoff).
    - [ ] Implement data loading (fragment metrics CSV).
    - [ ] Implement sampling interval calculation.
    - [ ] Implement robust maximum contiguous residence time calculation per fragment per region.
    - [ ] Implement histogram plotting for maximum residence times.
    - [ ] Define clear outputs (histogram plot file paths).
    - [ ] Add command-line arguments for all inputs and outputs.

### Phase 4: Pipeline Stage 3 - Binding Metrics Analysis (`analyze_binding_ratios.py`)

- [ ] Regenerate `analyze_binding_ratios.py`:
    - [x] Define clear inputs (fragment metrics CSV path, parameters like cutoff, temperature, min event duration, min event threshold).
    - [x] Implement data loading (fragment metrics CSV).
    - [x] Implement sampling interval calculation.
    - [x] Implement robust overall and per-facet K_D calculation (Time Off / Time On).
    - [x] Implement robust ΔG_D calculation (using T=453 K).
    - [x] Implement robust event counting per facet.
    - [x] Implement filtering of facets with N_events < 10.
    - [x] Implement robust mean residence time (τ) and standard deviation calculation per facet.
    - [x] Implement plotting for average K_D, ΔG_D, and mean τ per facet.
    - [x] Define clear outputs (plot file paths).
    - [x] Add command-line arguments for all inputs and outputs.

### Phase 5: Pipeline Stage 4 - Convergence Analysis (`analyze_convergence.py`)

- [x] Regenerate `analyze_convergence.py`:
- [x] Define clear inputs (fragment metrics CSV path).
- [x] Implement data loading (fragment metrics CSV).
- [x] Implement robust unique-facet coverage calculation and plotting.
- [x] Implement robust average Rg over time plotting.
- [x] Implement robust MSD calculation and plotting.
- [x] Implement robust event-count growth calculation and plotting.
- [x] Define clear outputs (plot file paths).
- [x] Add command-line arguments for all inputs and outputs.

### Phase 6: Pipeline Stage 5 - Advanced Visualization (`plot_residence_time_traces.py`)

- [ ] Regenerate `plot_residence_time_traces.py`:
    - [ ] Define clear inputs (fragment metrics CSV path, parameters for molecule selection).
    - [ ] Implement data loading (fragment metrics CSV).
    - [ ] Implement robust selection of "high-mover" molecules (by facet coverage and spatial excursion - potentially simplified initially).
    - [ ] Implement extraction of COM positions.
    - [ ] Implement robust 3D trajectory plotting (spaghetti plot), colored by time or facet ID.
    - [ ] Define clear outputs (plot file paths).
    - [ ] Add command-line arguments for all inputs and outputs.

### Phase 7: Integration, Documentation & Testing

- [ ] Create a master script or workflow manager to orchestrate the execution of the pipeline stages in the correct order, passing outputs as inputs.
- [ ] Ensure all scripts are runnable independently and as part of the pipeline.
- [ ] Update summary markdown generation in `nec_pt_fragment_metrics.py` to include details on all calculated metrics, thresholds, and generated plots from all pipeline stages.
- [ ] Add comprehensive docstrings and comments to all functions and scripts.
- [ ] Write a brief README.md explaining the pipeline architecture and how to run the analysis workflow.
- [ ] Test the full analysis pipeline with the sample data in `data/0H/`.
- [ ] Review generated plots and output files for correctness and consistency.

## Backlog / Future Enhancements

- [ ] Implement more sophisticated MSD calculation (averaging over multiple time origins).
- [ ] Add bootstrapping for error bars on K_D and ΔG_D.
- [ ] Implement more advanced criteria for selecting "high-mover" molecules.
- [ ] Add support for different simulation file formats in the data loading stage.
- [ ] Explore using a dedicated workflow management tool (e.g., Snakemake, Nextflow) for the pipeline orchestration.