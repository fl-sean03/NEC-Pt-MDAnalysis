# PLANNING.md: MDAnalysis Script Regeneration

## 1. Project Overview

- **Project Name and Description:** Regeneration of MDAnalysis data analysis scripts for molecular dynamics simulations. This project aims to rebuild the Python scripts used for analyzing simulation data, focusing on clarity, modularity, and adherence to best practices, based on a provided detailed analysis plan.
- **Objective:** To create a robust and maintainable set of scripts for calculating key metrics (residence times, dissociation constants, free energy of dissociation) and performing convergence checks and visualizations for molecular dynamics simulation data.

## 2. Vision & Objectives

- **Vision Statement:** To develop a robust, modular, and easily extensible data analysis pipeline for molecular dynamics simulations that accurately implements the detailed analysis plan and facilitates reproducible research.
- **Key Objectives:**
    - Design and implement a clear pipeline architecture where each script represents a distinct, modular stage with well-defined inputs and outputs.
    - Accurately implement all calculations and metrics defined in the detailed analysis plan (residence times, K_D, ΔG_D, mean τ, std τ, event counts).
    - Implement all specified sampling and convergence checks (unique-facet coverage, Rg over time, MSD, event-count growth).
    - Implement all specified visualization recipes (K_D plot, ΔG_D plot, mean τ plot, spaghetti plot, convergence plots).
    - Ensure scripts are robust with appropriate error handling within and between pipeline stages.
    - Ensure scripts are modular, well-commented, and follow Python best practices, with clear interfaces.
    - Make scripts runnable from the command line with clear arguments for input data and parameters (cutoffs, temperatures, thresholds).
    - Generate comprehensive output files, including a detailed summary document.

## 3. Architecture Overview: Data Analysis Pipeline

- **Pipeline Stages (Core Components):**
    - **Stage 1: Data Preparation (`nec_pt_fragment_metrics.py`):** Reads raw simulation data (PSF, DCD) and generates the base fragment metrics CSV and Pt classification JSON.
    - **Stage 2: Residence Time Analysis (`analyze_residence_times.py`, `analyze_max_residence_times.py`):** Reads the fragment metrics CSV and calculates/plots residence time distributions.
    - **Stage 3: Binding Metrics Analysis (`analyze_binding_ratios.py`):** Reads the fragment metrics CSV and calculates/plots K_D, ΔG_D, mean τ, and std τ.
    - **Stage 4: Convergence Analysis (`analyze_convergence.py`):** Reads the fragment metrics CSV and performs/plots convergence checks.
    - **Stage 5: Advanced Visualization (`plot_residence_time_traces.py`):** Reads the fragment metrics CSV and generates advanced visualizations like spaghetti plots.
- **Modularity and Robustness:**
    - Each script acts as a distinct module with defined inputs (e.g., CSV files) and outputs (e.g., other CSVs, plots).
    - Clear interfaces (command-line arguments) for passing data and parameters between stages.
    - Implementation will include error handling within each script to ensure robustness.
    - The modular design allows for easy addition or modification of pipeline stages.
- **Deliverable Documents:**
    - PLANNING.md (this document)
    - TASK.md (detailed task breakdown)
    - Regenerated analysis scripts (as listed in Pipeline Stages).
    - Generated output files (CSV, JSON, plots, summary markdown).
- **Technology Stack:** Python, MDAnalysis, NumPy, Pandas, Matplotlib, SciPy, scikit-learn (for DBSCAN if needed).
- **Constraints & Considerations:**
    - New code must be generated originally, referencing sample scripts for ideas/reference only.
    - Scripts should be runnable from the command line.
    - Adherence to the detailed analysis plan's definitions, formulas, and thresholds (0.1 ns min event length, T=453 K, N_events >= 10).

## 4. Milestones & Roadmap

- **Phase 1: Planning & Setup:**
    - Create PLANNING.md (Completed).
    - Create TASK.md.
    - Set up project directory structure.
- **Phase 2: Core Metric Calculation Scripts:**
    - Regenerate `nec_pt_fragment_metrics.py`.
    - Regenerate `analyze_residence_times.py` (including event counting and duration calculation).
    - Regenerate `analyze_binding_ratios.py` (including K_D, ΔG_D, mean τ, std τ, and filtering).
- **Phase 3: Convergence & Visualization Scripts:**
    - Regenerate `analyze_convergence.py`.
    - Regenerate `analyze_max_residence_times.py`.
    - Regenerate `plot_residence_time_traces.py`.
- **Phase 4: Integration, Documentation & Testing:**
    - Ensure scripts can be run sequentially as a workflow.
    - Update summary file generation in `nec_pt_fragment_metrics.py`.
    - Add command-line arguments and help messages to all scripts.
    - Write basic usage documentation.
    - Perform testing with sample data.

## 5. Project Organization & Workflow

- **Documentation Structure:** PLANNING.md and TASK.md will be in the root directory. Analysis scripts will be in the `src/` directory. Output files will be in `outputs/<run_id>/` subdirectories.
- **Workflow Overview:** Follow the tasks outlined in TASK.md. Implement scripts incrementally, testing each component as it is developed. Use command-line arguments for flexibility.

```

**Instructions:**
- Create the file `PLANNING.md` with the content provided above.
- Signal completion by using the `attempt_completion` tool. The `result` parameter should confirm that `PLANNING.md` has been created.
- These specific instructions supersede any conflicting general instructions for the Code mode.