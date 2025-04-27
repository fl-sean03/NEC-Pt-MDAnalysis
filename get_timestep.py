import MDAnalysis as mda
import sys

psf_path = 'data/0H/0HPt.psf'
dcd_path = 'data/0H/out_eq.dcd'

try:
    u = mda.Universe(psf_path, dcd_path)
    print(f"Frame time step (ns): {u.trajectory.dt}")
except FileNotFoundError as e:
    print(f"Error: File not found - {e}")
    sys.exit(1)
except Exception as e:
    print(f"An error occurred: {e}")
    sys.exit(1)