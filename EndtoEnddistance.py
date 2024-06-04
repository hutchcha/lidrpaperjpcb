import argparse
import MDAnalysis as mda
import numpy as np
from tqdm import tqdm
from MDAnalysis.analysis import distances

def calculate_end_to_end_distance(psf_file, dcd_file, output_file, resid1, resid2, step=1):
    """
    Calculate the end-to-end distance between two residues in a molecular dynamics simulation.

    Parameters:
    psf_file (str): Path to the PSF file.
    dcd_file (str): Path to the DCD file.
    output_file (str): Path to save the output distance file.
    resid1 (int): Residue ID of the first residue.
    resid2 (int): Residue ID of the last residue.
    step (int): Step size for trajectory iteration.
    """
    universe = mda.Universe(psf_file, dcd_file)
    distance_array = np.array([])

    ag1 = universe.select_atoms(f"(protein and resid {resid1}) and name CA")
    print(f"Number of atoms in ag1: {ag1.n_atoms}")
    ag2 = universe.select_atoms(f"resname {resid2} and name CA")
    print(f"Number of atoms in ag2: {ag2.n_atoms}")

    for ts in tqdm(universe.trajectory[::step], desc="Distance Measurement", unit="frame"):
        dis = distances.distance_array(ag1.positions, ag2.positions)
        avg_distance = np.mean(dis)
        distance_array = np.concatenate((distance_array, [avg_distance]), axis=0)
    
    np.savetxt(output_file, distance_array)

def main():
    parser = argparse.ArgumentParser(description="Calculate end-to-end distance in a molecular dynamics simulation.")
    parser.add_argument('--psf_file', type=str, help='Path to the PSF file', required=True)
    parser.add_argument('--dcd_file', type=str, help='Path to the DCD file', required=True)
    parser.add_argument('--output_file', type=str, help='Path to save the output distance file', required=True)
    parser.add_argument('--resid1', type=int, help='Residue ID of the first residue', required=True)
    parser.add_argument('--resid2', type=int, help='Residue ID of the second residue', required=True)
    parser.add_argument('--step', type=int, help='Step size for trajectory iteration', default=1)

    args = parser.parse_args()
    calculate_end_to_end_distance(args.psf_file, args.dcd_file, args.output_file, args.resid1, args.resid2, args.step)

if __name__ == "__main__":
    main()