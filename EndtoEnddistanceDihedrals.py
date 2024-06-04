import argparse
import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import distances
from MDAnalysis.analysis.dihedrals import Dihedral
import pandas as pd
from tqdm import tqdm

def measure_end_to_end_distance(universe, resid1, resid2, step=1):
    """
    Measure end-to-end distance between two residues in a molecular dynamics simulation.

    Parameters:
    universe (MDAnalysis.Universe): The MDAnalysis universe.
    resid1 (str): Selection string for the first residue.
    resid2 (str): Selection string for the second residue.
    step (int): Step size for trajectory iteration.

    Returns:
    np.ndarray: Array of end-to-end distances.
    """
    distance_array = np.array([])
    ag1 = universe.select_atoms(resid1)
    ag2 = universe.select_atoms(resid2)
    for ts in tqdm(universe.trajectory[::step], desc="End-to-End Distance Measurement", unit="frame"):
        dis = distances.distance_array(ag1.positions, ag2.positions)
        avg_distance = np.mean(dis)
        distance_array = np.concatenate((distance_array, [avg_distance]), axis=0)
    return distance_array

def measure_dihedrals(universe, resid_range, step=1):
    """
    Measure dihedrals for a chosen residue sequence in a molecular dynamics simulation.

    Parameters:
    universe (MDAnalysis.Universe): The MDAnalysis universe.
    resid_range (str): Selection string for the residue range.
    step (int): Step size for trajectory iteration.

    Returns:
    np.ndarray: Array of dihedral angles.
    """
    ags = [universe.select_atoms(resid_range)]
    R = Dihedral(ags).run(step=step, verbose=True)
    angles = np.reshape(R.results.angles, len(R.results.angles))
    return angles

def main():
    parser = argparse.ArgumentParser(description="Measure end-to-end distances and dihedrals in a molecular dynamics simulation.")
    parser.add_argument('--psf_file', type=str, help='Path to the PSF file', required=True)
    parser.add_argument('--dcd_file', type=str, help='Path to the DCD file', required=True)
    parser.add_argument('--output_file', type=str, help='Path to save the output file', required=True)
    parser.add_argument('--measure', type=str, choices=['distance', 'dihedral'], help='Measurement to perform: "distance" or "dihedral"', required=True)
    parser.add_argument('--resid1', type=str, help='Selection string for the first residue (for distance measurement)', required=False)
    parser.add_argument('--resid2', type=str, help='Selection string for the second residue (for distance measurement)', required=False)
    parser.add_argument('--resid_range', type=str, help='Selection string for the residue range (for dihedral measurement)', required=False)
    parser.add_argument('--step', type=int, help='Step size for trajectory iteration', default=1)

    args = parser.parse_args()

    universe = mda.Universe(args.psf_file, args.dcd_file)

    if args.measure == 'distance':
        if not args.resid1 or not args.resid2:
            print("Error: --resid1 and --resid2 must be provided for distance measurement")
            return
        distances = measure_end_to_end_distance(universe, args.resid1, args.resid2, args.step)
        np.savetxt(args.output_file, distances)
    elif args.measure == 'dihedral':
        if not args.resid_range:
            print("Error: --resid_range must be provided for dihedral measurement")
            return
        dihedrals = measure_dihedrals(universe, args.resid_range, args.step)
        np.savetxt(args.output_file, dihedrals)

if __name__ == "__main__":
    main()
