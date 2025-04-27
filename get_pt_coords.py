import MDAnalysis as mda
import sys

psf_path = 'data/0H/0HPt.psf'
pdb_path = 'data/0H/0HPt.pdb'

try:
    u = mda.Universe(psf_path, pdb_path)
    pt_atoms = u.select_atoms("name PT*")

    print("Pt atom indices and coordinates:")
    for atom in pt_atoms:
        print(f"Index {atom.index}: ({atom.position[0]:.3f}, {atom.position[1]:.3f}, {atom.position[2]:.3f}) Å")

except FileNotFoundError as e:
    print(f"Error: File not found - {e}", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"An error occurred: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc(file=sys.stderr)
    sys.exit(1)