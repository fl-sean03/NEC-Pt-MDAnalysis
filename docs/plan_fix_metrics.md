# Plan to Fix and Debug Metrics Calculation

This document outlines the steps to address the discrepancy in the `min_dist` and `nearest_pt_class` metrics calculated by `nec_pt_fragment_metrics.py` and ensure they align with a center-of-mass (COM) based calculation. This is necessary to correct the trace coloring and contacted regions count in the `plot_residence_time_traces.py` script.

## Plan Steps:

1.  **Modify `nec_pt_fragment_metrics.py`:**
    *   Update the `fragment_metrics` function to calculate `min_dist` as the minimum distance from the fragment's COM to all Pt atoms (or specifically surface Pt atoms, depending on the desired definition).
    *   Update the logic for determining `nearest_pt_idx` and `nearest_pt_class` to be based on the Pt atom nearest to the fragment's COM.
    *   Ensure the updated calculations are correct and handle edge cases (e.g., no Pt atoms).

2.  **Regenerate Fragment Metrics CSV:**
    *   Run the modified `nec_pt_fragment_metrics.py` script with the sample data (`data/0H/0HPt.psf`, `data/0H/out_eq.dcd`) and output the new `fragment_metrics.csv` and `pt_classification.json` files to a new directory (e.g., `outputs/com_metrics_run/`).

3.  **Verify New Fragment Metrics:**
    *   Run the `verify_min_dist_and_region.py` script using the newly generated `fragment_metrics.csv` and `pt_classification.json` files as input.
    *   Analyze the output of `verify_min_dist_and_region.py` to confirm that the `min_dist` and `nearest_pt_class` values in the new CSV now match the COM-based calculations performed by the verification script. Address any remaining discrepancies.

4.  **Re-run `plot_residence_time_traces.py`:**
    *   Run the `plot_residence_time_traces.py` script using the newly generated and verified `fragment_metrics.csv` and `pt_classification.json` files as input.
    *   Use an appropriate `distance_cutoff` (e.g., 8.0 Å or a value determined from the new `min_dist` distribution).

5.  **Visually Inspect Plots:**
    *   Visually inspect the generated plots in the `outputs/com_metrics_run/visualization_plots/` directory to confirm that the trace coloring and contacted regions count are now correct and align with expectations.

6.  **Update Documentation:**
    *   Update `nec_pt_fragment_metrics.py` docstrings and comments to reflect the change in how `min_dist` and `nearest_pt_class` are calculated.
    *   Update `PLANNING.md` and `TASK.md` if necessary to clarify the definition of these metrics.

This plan will ensure that the input data for the plotting script is accurate according to a COM-based definition, which should resolve the issues observed in the plots.