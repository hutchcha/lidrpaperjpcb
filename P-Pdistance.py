import argparse
import MDAnalysis as mda
import numpy as np
from tqdm import tqdm

def phosphorous_distance(psf_file, dcd_file, output_file, step=1):
    """
    Calculate the distance between phosphorous groups in a lipid membrane.

    Parameters:
    psf_file (str): Path to the PSF file.
    dcd_file (str): Path to the DCD file.
    output_file (str): Path to save the output distance file.
    step (int): Step size for trajectory iteration.
    """
    universe = mda.Universe(psf_file, dcd_file)
    d = np.array([])

    for ts in tqdm(universe.trajectory[::step], desc="P-P Distance Progress"):
        memcen = universe.select_atoms("resname POPC POPE POPS SAPI CHL1 PSM PLA18").center_of_mass()[2]
        top = universe.select_atoms(f"name P and prop z > {memcen}")
        bot = universe.select_atoms(f"name P and prop z < {memcen}")
        avgz1 = top.center_of_mass()
        avgz2 = bot.center_of_mass()
        dis = avgz1 - avgz2
        z = dis[2]
        a = np.array(z).reshape((1,))
        d = np.concatenate((d, a), axis=0)

    np.savetxt(output_file, d)

def main():
    parser = argparse.ArgumentParser(description="Calculate P-P distance in a lipid membrane.")
    parser.add_argument('--psf_file', type=str, help='Path to the PSF file', required=True)
    parser.add_argument('--dcd_file', type=str, help='Path to the DCD file', required=True)
    parser.add_argument('--output_file', type=str, help='Path to save the output distance file', required=True)
    parser.add_argument('--step', type=int, help='Step size for trajectory iteration', default=1)

    args = parser.parse_args()
    phosphorous_distance(args.psf_file, args.dcd_file, args.output_file, args.step)

if __name__ == "__main__":
    main()
