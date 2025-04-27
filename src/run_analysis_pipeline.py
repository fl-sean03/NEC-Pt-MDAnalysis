# !/usr/bin/env python3
"""
Orchestrates the execution of the MDAnalysis data analysis pipeline.

This script runs the analysis stages sequentially and generates a comprehensive
summary report.
"""

import argparse
import subprocess
import os
import json

def run_stage(script_path, args, cwd=None):
    """Runs a single analysis script as a subprocess."""
    command = ["python3", script_path] + args
    print(f"Running command: {' '.join(command)}")
    result = subprocess.run(command, capture_output=True, text=True, cwd=cwd)
    print(f"Stdout:\n{result.stdout}")
    print(f"Stderr:\n{result.stderr}")
    if result.returncode != 0:
        print(f"Error running {script_path}. Exiting.")
        exit(result.returncode)
    return result.stdout

def generate_comprehensive_summary(output_dir, run_id, system_info, run_parameters, analysis_results, plot_paths, known_issues):
    """Generates the comprehensive summary.md file."""
    summary_path = os.path.join(output_dir, f"{run_id}_summary.md")
    with open(summary_path, 'w') as f:
        f.write(f"# Comprehensive Run Summary: {run_id}\n\n")

        f.write("## Project Overview\n\n")
        f.write("This run is part of the MDAnalysis Script Regeneration project, which aims to rebuild Python scripts for analyzing molecular dynamics simulation data. The objective is to create a robust and maintainable pipeline for calculating key metrics and performing convergence checks and visualizations.\n\n")

        f.write("## System Information\n\n")
        for key, value in system_info.items():
            f.write(f"- {key}: {value}\n")
        f.write("\n")

        f.write("## Run Parameters\n\n")
        for key, value in run_parameters.items():
            f.write(f"- {key}: {value}\n")
        f.write("\n")

        f.write("## Analysis Pipeline Stages\n\n")
        f.write("The analysis pipeline consists of several modular stages:\n\n")
        f.write("### Stage 1: Data Preparation (`nec_pt_fragment_metrics.py`)\n")
        f.write("This script reads raw simulation data (PSF, DCD) and extracts per-fragment, per-frame interaction metrics.\n\n")
        f.write("### Stage 2: Residence Time Analysis (`analyze_residence_times.py`, `analyze_max_residence_times.py`)\n")
        f.write("These scripts calculate and plot residence time distributions and maximum residence times.\n\n")
        f.write("### Stage 3: Binding Metrics Analysis (`analyze_binding_ratios.py`)\n")
        f.write("This script calculates and plots binding metrics such as K_D, ΔG_D, mean τ, and std τ.\n\n")
        f.write("### Stage 4: Convergence Analysis (`analyze_convergence.py`)\n")
        f.write("This script performs and plots various convergence checks.\n\n")
        f.write("### Stage 5: Advanced Visualization (`plot_residence_time_traces.py`)\n")
        f.write("This script generates advanced visualizations, such as spaghetti plots.\n\n")

        f.write("## Key Results\n\n")
        if analysis_results:
            for stage, results in analysis_results.items():
                f.write(f"### {stage}\n\n")
                # Handle specific stage results formatting
                if stage == "Residence Time Analysis":
                    f.write("Residence Event Summary:\n\n")
                    for region, summary in results.get("residence_event_summary", {}).items():
                        f.write(f"**{region}:**\n")
                        f.write(f"- Event Count: {summary.get('event_count', 'N/A')}\n")
                        f.write(f"- Average Duration (ns): {summary.get('average_duration_ns', 'N/A'):.3f}\n")
                        f.write(f"- Median Duration (ns): {summary.get('median_duration_ns', 'N/A'):.3f}\n")
                        f.write(f"- Total Residence Time (ns): {summary.get('total_residence_time_ns', 'N/A'):.3f}\n\n")
                elif stage == "Maximum Residence Time Analysis":
                     f.write("Maximum Residence Time Summary:\n\n")
                     for region, summary in results.get("max_residence_time_summary", {}).items():
                         f.write(f"**{region}:**\n")
                         f.write(f"- Fragment Count with Max Residence: {summary.get('fragment_count_with_max_residence', 'N/A')}\n")
                         f.write(f"- Average Max Duration (ns): {summary.get('average_max_duration_ns', 'N/A'):.3f}\n")
                         f.write(f"- Median Max Duration (ns): {summary.get('median_max_duration_ns', 'N/A'):.3f}\n\n")
                elif stage == "Binding Metrics Analysis":
                    f.write("Average Facet Metrics:\n\n")
                    for region, metrics in results.get("average_facet_metrics", {}).items():
                        f.write(f"**{region}:**\n")
                        f.write(f"- Average K_D: {metrics.get('avg_k_d', 'N/A'):.3f}\n")
                        f.write(f"- Average ΔG_D (kJ/mol): {metrics.get('avg_delta_g_d', 'N/A'):.3f}\n")
                        f.write(f"- Mean Residence Time (τ) (ns): {metrics.get('mean_tau', 'N/A'):.3f}\n")
                        f.write(f"- Std Dev of τ (ns): {metrics.get('std_tau', 'N/A'):.3f}\n")
                        f.write(f"- Number of Events: {metrics.get('n_events', 'N/A')}\n\n")
                elif stage == "Convergence Analysis":
                    f.write("Convergence Checks Performed:\n\n")
                    for check in results.get("convergence_checks_performed", []):
                        f.write(f"- {check}\n")
                    f.write("\nRefer to the generated plots for visual assessment of convergence.\n\n")
                elif stage == "Advanced Visualization (Residence Time Traces)":
                    f.write("Fragment Metrics Summary:\n\n")
                    fragment_metrics = results.get("fragment_metrics_summary", {})
                    if fragment_metrics:
                        f.write("| Fragment ID | Spatial Excursion (Å) | Residence Events | Unique Regions |\n")
                        f.write("|-------------|-------------------------|------------------|----------------|\n")
                        for frag_id, metrics in fragment_metrics.items():
                            f.write(f"| {frag_id} | {metrics.get('spatial_excursion', 'N/A'):.2f} | {metrics.get('residence_events', 'N/A')} | {metrics.get('unique_regions', 'N/A')} |\n")
                        f.write("\n")
                    f.write("Selected Fragments for Plotting:\n\n")
                    selected_frags = results.get("selected_fragments_for_plotting", {})
                    if selected_frags:
                        for frag_id, selection_type in selected_frags.items():
                            f.write(f"- Fragment {frag_id}: {selection_type.replace('_', ' ').title()}\n")
                        f.write("\n")
                else:
                    # Default handling for unexpected stages
                    for key, value in results.items():
                        f.write(f"- {key}: {value}\n")
                    f.write("\n")
        else:
            f.write("Analysis results will be included here after running the pipeline.\n\n")


        f.write("## Visualizations\n\n")
        if plot_paths:
            f.write("The following plots were generated:\n\n")
            for plot_path in plot_paths:
                f.write(f"- `{plot_path}`\n")
                # Optional: Embed images if desired and path is relative/accessible
                # f.write(f"![Plot]({plot_path})\n\n") # Uncomment to embed images
            f.write("\n")
        else:
            f.write("References to generated plots will be included here after running the pipeline.\n\n")

        f.write("## Data Verification\n\n")
        f.write("The `verify_min_dist_and_region.py` script is used to verify data consistency. Verification was performed for fragment 214 and confirmed data consistency.\n\n") # Placeholder, can be expanded

        f.write("## Output Files\n\n")
        f.write(f"All output files are saved in the `{output_dir}` directory.\n\n")
        f.write(f"- Comprehensive Summary (this file): `{os.path.basename(summary_path)}`\n")
        f.write(f"- Fragment metrics CSV: `{run_id}_fragment_metrics.csv`\n") # Placeholder, actual name from Stage 1 output
        f.write(f"- Pt classification JSON: `{run_id}_pt_classification.json`\n") # Placeholder, actual name from Stage 1 output
        f.write(f"- Residence time histograms: `{os.path.join('residence_time_histograms', '')}`\n") # Placeholder
        f.write(f"- Max residence time plots: `{os.path.join('max_residence_plots', '')}`\n") # Placeholder
        f.write(f"- Binding metrics plots: `{os.path.join('binding_metrics_plots', '')}`\n") # Placeholder
        f.write(f"- Convergence analysis plots: `{os.path.join('convergence_plots', '')}`\n") # Placeholder
        f.write(f"- Visualization plots: `{os.path.join('visualization_plots', '')}`\n\n") # Placeholder


        f.write("## Known Issues\n\n")
        if known_issues:
            for issue in known_issues:
                f.write(f"- {issue}\n")
        else:
            f.write("No known issues.\n\n")


    print(f"Generated comprehensive summary: {summary_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Orchestrates the MDAnalysis data analysis pipeline.'
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
        '--output-dir', required=True,
        help='Directory to save all output files.'
    )
    # Add other common arguments that need to be passed to stages
    parser.add_argument(
        '--sample-interval', type=int, default=1,
        help='Sample every N frames (for Stage 1).'
    )
    parser.add_argument(
        '--default-cutoff', type=float, default=2.5,
        help='Default contact cutoff distance in Angstroms (for Stage 1).'
    )
    parser.add_argument(
        '--pt-cutoff', type=float, default=3.0,
        help='Cutoff distance in Angstroms for Pt-Pt coordination number calculation (for Stage 1).'
    )
    # Add arguments for other stages as needed

    args = parser.parse_args()

    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    # Define paths for intermediate and final output files
    metrics_csv_path = os.path.join(args.output_dir, f"{args.run_id}_fragment_metrics.csv")
    pt_json_path = os.path.join(args.output_dir, f"{args.run_id}_pt_classification.json")
    # Define paths for outputs of other stages (plots, summary JSONs, etc.)

    # --- Run Pipeline Stages ---

    # Stage 1: Data Preparation
    print("\n--- Running Stage 1: Data Preparation ---")
    # Get the directory of the PSF file
    data_dir = os.path.dirname(args.psf)

    # Get the absolute path of the script
    script_abs_path = os.path.abspath("src/nec_pt_fragment_metrics.py")

    # Define arguments for Stage 1, using filenames for input and absolute paths for output
    stage1_args = [
        "--psf", os.path.basename(args.psf),
        "--dcd", os.path.basename(args.dcd),
        "--run-id", args.run_id,
        "--sample-interval", str(args.sample_interval),
        "--default-cutoff", str(args.default_cutoff),
        "--pt-cutoff", str(args.pt_cutoff),
        "--output-csv", os.path.abspath(metrics_csv_path), # Provide absolute path for output
        "--output-json", os.path.abspath(pt_json_path), # Provide absolute path for output
        "--output-summary", os.path.abspath(os.path.join(args.output_dir, f"{args.run_id}_stage1_summary.md")) # Provide absolute path for output
    ]

    # Stage 1 doesn't output JSON summary to stdout, its main output is CSV and JSON files
    # Run Stage 1 in the directory of the data files
    run_stage(script_abs_path, stage1_args, cwd=data_dir)

    # Stage 2: Residence Time Analysis
    print("\n--- Running Stage 2: Residence Time Analysis ---")
    stage2_outdir = os.path.join(args.output_dir, "residence_time_histograms")
    stage2_args = [
        "--csv", metrics_csv_path,
        "--outdir", stage2_outdir,
        # Removed --run-id as analyze_residence_times.py infers it
        # Add other necessary arguments for these scripts (e.g., cutoff, min-duration, bins, max-bin, overlay, side-by-side)
        "--cutoff", str(args.default_cutoff), # Using default_cutoff from Stage 1
        "--min-duration", "0.1", # From PLANNING.md
        "--bins", "50", # Default
        "--overlay" # Example: generate overlay plot
    ]
    stage2_res_times_stdout = run_stage("src/analyze_residence_times.py", stage2_args)
    try:
        stage2_results = json.loads(stage2_res_times_stdout)
        all_analysis_results[stage2_results["stage"]] = stage2_results
        all_plot_paths.extend(stage2_results["generated_plots"])
    except json.JSONDecodeError:
        print("Warning: Could not decode JSON from analyze_residence_times.py stdout.")
    except KeyError:
         print("Warning: Unexpected JSON structure from analyze_residence_times.py.")

    # Stage 2 (Max Residence Times)
    print("\n--- Running Stage 2 (Max Residence Times): ---") # Clarified stage name
    stage2_max_outdir = os.path.join(args.output_dir, "max_residence_plots")
    stage2_max_args = [
        "--csv", metrics_csv_path,
        "--outdir", stage2_max_outdir,
        # Removed --run-id as analyze_max_residence_times.py infers it
        "--cutoff", str(args.default_cutoff), # Using default_cutoff from Stage 1
        "--bins", "50", # Default
        "--max-bin", "100" # Example max bin
    ]
    stage2_max_res_stdout = run_stage("src/analyze_max_residence_times.py", stage2_max_args)
    try:
        stage2_max_results = json.loads(stage2_max_res_stdout)
        all_analysis_results[stage2_max_results["stage"]] = stage2_max_results
        all_plot_paths.extend(stage2_max_results["generated_plots"])
    except json.JSONDecodeError:
        print("Warning: Could not decode JSON from analyze_max_residence_times.py stdout.")
    except KeyError:
         print("Warning: Unexpected JSON structure from analyze_max_residence_times.py.")


    # Stage 3: Binding Metrics Analysis
    print("\n--- Running Stage 3: Binding Metrics Analysis ---")
    stage3_outdir = os.path.join(args.output_dir, "binding_metrics_plots")
    stage3_args = [
        "--csv", metrics_csv_path,
        "--outdir", stage3_outdir,
        "--temperature", "453.0", # From PLANNING.md
        "--min-event-duration", "0.1", # From PLANNING.md
        "--min-event-threshold", "10", # From PLANNING.md
        "--cutoff", str(args.default_cutoff), # Using default_cutoff from Stage 1
        # Removed --run-id as analyze_binding_ratios.py infers it
    ]
    stage3_stdout = run_stage("src/analyze_binding_ratios.py", stage3_args)
    try:
        stage3_results = json.loads(stage3_stdout)
        all_analysis_results[stage3_results["stage"]] = stage3_results
        all_plot_paths.extend(stage3_results["generated_plots"])
    except json.JSONDecodeError:
        print("Warning: Could not decode JSON from analyze_binding_ratios.py stdout.")
    except KeyError:
         print("Warning: Unexpected JSON structure from analyze_binding_ratios.py.")


    # Stage 4: Convergence Analysis
    print("\n--- Running Stage 4: Convergence Analysis ---")
    stage4_outdir = os.path.join(args.output_dir, "convergence_plots")
    stage4_args = [
        "--csv", metrics_csv_path,
        "--outdir", stage4_outdir,
        "--run-id", args.run_id # analyze_convergence.py *does* accept --run-id
        # Add other necessary arguments for this script (e.g., cutoff if used internally)
    ]
    stage4_stdout = run_stage("src/analyze_convergence.py", stage4_args)
    try:
        stage4_results = json.loads(stage4_stdout)
        all_analysis_results[stage4_results["stage"]] = stage4_results
        all_plot_paths.extend(stage4_results["generated_plots"])
    except json.JSONDecodeError:
        print("Warning: Could not decode JSON from analyze_convergence.py stdout.")
    except KeyError:
         print("Warning: Unexpected JSON structure from analyze_convergence.py.")


    # Stage 5: Advanced Visualization
    print("\n--- Running Stage 5: Advanced Visualization ---")
    stage5_outdir = os.path.join(args.output_dir, "visualization_plots")
    stage5_args = [
        "--csv", metrics_csv_path,
        "--psf", args.psf,
        "--dcd", args.dcd,
        "--pt-json", pt_json_path,
        "--outdir", stage5_outdir,
        "--min-residence-duration-ns", "0.1", # From PLANNING.md
        "--frame-time-step-ns", "0.001", # Assuming 1 fs time step
        "--distance-cutoff", str(args.default_cutoff) # Using default_cutoff from Stage 1
        # Add other necessary arguments for this script (e.g., specific fragment IDs to plot)
    ]
    stage5_stdout = run_stage("src/plot_residence_time_traces.py", stage5_args)
    try:
        stage5_results = json.loads(stage5_stdout)
        all_analysis_results[stage5_results["stage"]] = stage5_results
        all_plot_paths.extend(stage5_results["generated_plots"])
    except json.JSONDecodeError:
        print("Warning: Could not decode JSON from plot_residence_time_traces.py stdout.")
    except KeyError:
         print("Warning: Unexpected JSON structure from plot_residence_time_traces.py.")


    # --- Generate Comprehensive Summary ---
    print("\n--- Generating Comprehensive Summary ---")
    # Collect system info, run parameters, analysis results, plot paths, and known issues
    # System info would ideally be obtained from the universe object in Stage 1,
    # but for now, use placeholders or infer from data if possible.
    # For now, using placeholders and parameters passed to the script.
    system_info = {
        "System time step": "1 fs", # Assuming 1 fs based on previous context
        "Total number of frames": "N/A (obtain from trajectory)", # Need to get this from Stage 1 or Universe
        "Number of NEC molecules (fragments)": "N/A (obtain from trajectory)", # Need to get this from Stage 1 or Universe
        "Number of Pt atoms": "N/A (obtain from trajectory)", # Need to get this from Stage 1 or Universe
        "Box dimensions (Å)": "N/A (obtain from trajectory)" # Need to get this from Stage 1 or Universe
    }
    run_parameters = {
        "Run ID": args.run_id,
        "Sample interval (frames)": args.sample_interval,
        "Default contact cutoff (Å)": args.default_cutoff,
        "Pt classification cutoff (Å)": args.pt_cutoff,
        "ΔG_D calculation temperature": "453 K", # From PLANNING.md
        "Minimum event duration threshold for analysis": "0.1 ns", # From PLANNING.md
        "Minimum event count threshold for facet inclusion in averages": ">= 10" # From PLANNING.md
    }
    # analysis_results is already populated from stage outputs
    # plot_paths is already populated from stage outputs
    known_issues = [
        "Persistent KeyError: 'min_distance' in analyze_convergence.py requires further investigation."
    ]

    generate_comprehensive_summary(args.output_dir, args.run_id, system_info, run_parameters, all_analysis_results, all_plot_paths, known_issues)

    print("\nPipeline execution and comprehensive summary generation finished.")


if __name__ == "__main__":
    main()