import MDAnalysis as mda
import sys
import traceback

psf_path = 'data/0H/0HPt.psf'
dcd_path = 'data/0H/out_eq.dcd'
fragment_id_to_check = 20 # Fragment ID for "most unique regions" from previous run
sample_interval = 100

try:
    u = mda.Universe(psf_path, dcd_path)
    
    # Select the fragment by its fragment index (MDAnalysis fragment indices are 0-based)
    # The fragment_id in the CSV is 1-based, so subtract 1 for MDAnalysis selection
    try:
        fragment_to_check = u.select_atoms("not name PT*").fragments[fragment_id_to_check - 1]
    except IndexError:
        print(f"Error: Fragment with ID {fragment_id_to_check} not found.", file=sys.stderr)
        sys.exit(1)

    print(f"Checking location for Fragment ID {fragment_id_to_check} every {sample_interval} frames:")

    for ts in u.trajectory[::sample_interval]:
        com = fragment_to_check.center_of_mass()
        print(f"Frame {ts.frame}: COM = ({com[0]:.3f}, {com[1]:.3f}, {com[2]:.3f}) Å")

except FileNotFoundError as e:
    print(f"Error: File not found - {e}", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"An error occurred: {e}", file=sys.stderr)
    traceback.print_exc(file=sys.stderr)
    sys.exit(1)