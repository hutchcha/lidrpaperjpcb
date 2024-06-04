import argparse
import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import rms

def calculate_rmsd(psf_file, dcd_file, output_file, select, step=1, ref_frame=1):
    """
    Calculate the RMSD for selected atoms in a molecular dynamics simulation.

    Parameters:
    psf_file (str): Path to the PSF file.
    dcd_file (str): Path to the DCD file.
    output_file (str): Path to save the output RMSD file.
    select (str): Atom selection string.
    step (int): Step size for trajectory iteration.
    ref_frame (int): Reference frame for RMSD calculation.
    """
    universe = mda.Universe(psf_file, dcd_file)
    R = rms.RMSD(universe, universe, select=select, step=step, ref_frame=ref_frame, verbose=True)
    R.run()
    rmsd = R.results.rmsd.T
    np.savetxt(output_file, rmsd)

def main():
    parser = argparse.ArgumentParser(description="Calculate RMSD for selected atoms in a molecular dynamics simulation.")
    parser.add_argument('--psf_file', type=str, help='Path to the PSF file', required=True)
    parser.add_argument('--dcd_file', type=str, help='Path to the DCD file', required=True)
    parser.add_argument('--output_file', type=str, help='Path to save the output RMSD file', required=True)
    parser.add_argument('--select', type=str, help='Atom selection string', required=True)
    parser.add_argument('--step', type=int, help='Step size for trajectory iteration', default=1)
    parser.add_argument('--ref_frame', type=int, help='Reference frame for RMSD calculation', default=1)

    args = parser.parse_args()
    calculate_rmsd(args.psf_file, args.dcd_file, args.output_file, args.select, args.step, args.ref_frame)

if __name__ == "__main__":
    main()
