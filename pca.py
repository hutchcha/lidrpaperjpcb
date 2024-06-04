import argparse
import matplotlib.pyplot as plt
import numpy as np
import pyemma

def calculate_pca(top_file, traj_file, output_file, frames=None, dim=2):
    """
    Calculate PCA using PyEMMA for a given molecular dynamics trajectory.

    Parameters:
    top_file (str): Path to the topology file.
    traj_file (str): Path to the trajectory file.
    output_file (str): Path to save the PCA results.
    frames (list of int, optional): Specific frames to load.
    dim (int): Number of PCA dimensions to calculate.
    """
    positions_feat = pyemma.coordinates.featurizer(top_file)
    positions_feat.add_distances_ca()

    if frames is not None:
        positions_data = pyemma.coordinates.load(traj_file, features=positions_feat, frames=frames)
    else:
        positions_data = pyemma.coordinates.load(traj_file, features=positions_feat)

    pca = pyemma.coordinates.pca(positions_data, dim=dim)
    pca_output = pca.get_output()
    pca_con = np.concatenate(pca_output)

    np.save(output_file, pca_con)

def main():
    parser = argparse.ArgumentParser(description="Calculate PCA for a molecular dynamics trajectory using PyEMMA.")
    parser.add_argument('--top_file', type=str, help='Path to the topology file', required=True)
    parser.add_argument('--traj_file', type=str, help='Path to the trajectory file', required=True)
    parser.add_argument('--output_file', type=str, help='Path to save the PCA results', required=True)
    parser.add_argument('--frames', type=int, nargs='*', help='Specific frames to load', default=None)
    parser.add_argument('--dim', type=int, help='Number of PCA dimensions to calculate', default=2)

    args = parser.parse_args()
    calculate_pca(args.top_file, args.traj_file, args.output_file, args.frames, args.dim)

if __name__ == "__main__":
    main()
