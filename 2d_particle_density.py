import argparse
import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from MDAnalysis.analysis.leaflet import LeafletFinder
import pickle
from scipy.stats import gaussian_kde

def load_universes(psf_files, dcd_files):
    """Load MDAnalysis universes for different systems."""
    universes = []
    for psf, dcd in zip(psf_files, dcd_files):
        universes.append(mda.Universe(psf, dcd))
    return universes

def get_lipid_Pdistances(lipid, universe, frame):
    """
    Calculate the 2D positions of lipid phosphate (P) atoms relative to the protein center of mass.
    
    Parameters:
    lipid (str): The lipid type.
    universe (MDAnalysis.Universe): The MDAnalysis universe.
    frame (int): The frame number in the trajectory.

    Returns:
    np.ndarray: 2D positions of lipid P atoms.
    """
    # Calculate membrane center of mass, replace with desired lipids
    memcom = universe.select_atoms("resname POPC POPE POPS PSM PLA18 CHL1 and not name H*").center_of_mass()[2]
    # Select protein atoms, replace with proper protein selection
    protein = universe.select_atoms("resname CYSF CYSG GLYM and name CA")
    pro_com = protein.center_of_mass()
    # Select lipid P atoms in the upper leaflet
    ag1 = universe.select_atoms(f"resname {lipid} and name P O3 and prop z > {memcom} and prop z > {pro_com[2] + 10}")
    # Calculate positions relative to the protein center of mass
    p_xypositions = ag1.positions[:, 0:2] - pro_com[0:2]
    
    return p_xypositions

def compute_2Dgaussian_KDE(dataset, grid_length, x_min, x_max, y_min, y_max):
    """
    Compute the 2D Gaussian KDE for a given dataset.
    
    Parameters:
    dataset (np.ndarray): The dataset for which to compute the KDE.
    grid_length (int): The number of grid points along each axis.
    x_min (float): Minimum value for x-axis grid.
    x_max (float): Maximum value for x-axis grid.
    y_min (float): Minimum value for y-axis grid.
    y_max (float): Maximum value for y-axis grid.

    Returns:
    np.ndarray: X, Y grid coordinates and the KDE values Z.
    """
    xy_values = dataset.T
    kde = gaussian_kde(xy_values)
    x_grid = np.linspace(x_min, x_max, grid_length)
    y_grid = np.linspace(y_min, y_max, grid_length)
    X, Y = np.meshgrid(x_grid, y_grid)
    Z = kde(np.vstack([X.ravel(), Y.ravel()]))
    Z = np.reshape(Z, X.shape)
    return X, Y, Z

def analyze_universe(universe, lipids, name, grid_length):
    """
    Analyze lipid distributions in a given universe and compute KDE.
    
    Parameters:
    universe (MDAnalysis.Universe): The MDAnalysis universe.
    lipids (list): List of lipid types.
    name (str): The name of the universe/system for output.
    grid_length (int): The number of grid points for KDE computation.
    """
    clustering_dic = {}
    kde_dic = {}

    # Calculate average box vectors
    x_lengths = []
    y_lengths = []
    for ts in universe.trajectory[::1]:
        x_lengths.append(universe.dimensions[0])
        y_lengths.append(universe.dimensions[1])
    avg_x = np.mean(x_lengths)
    avg_y = np.mean(y_lengths)
    x_min, x_max = -avg_x / 2, avg_x / 2
    y_min, y_max = -avg_y / 2, avg_y / 2

    print(name)
    for lipid in lipids:
        print(f"Analyzing lipid: {lipid}")
        lipid_atoms = universe.select_atoms(f"resname {lipid} and name P O3")
        if len(lipid_atoms.resnames) == 0:
            print(f"No {lipid} detected, skipping...")
            continue
        xy_array = []
        for ts in tqdm(universe.trajectory[::1], desc=f"{lipid}"):
            xy_pos = get_lipid_Pdistances(lipid, universe, ts)
            xy_array.append(xy_pos)
        combined_xy = np.vstack(xy_array)
        clustering_dic[lipid] = combined_xy
        X, Y, Z = compute_2Dgaussian_KDE(combined_xy, grid_length, x_min, x_max, y_min, y_max)
        kde_dic[lipid] = (X, Y, Z)
    save_results(clustering_dic, kde_dic, name)

def save_results(clustering_dic, kde_dic, name):
    """
    Save the clustering results and KDE data to pickle files.
    
    Parameters:
    clustering_dic (dict): Dictionary of clustering data.
    kde_dic (dict): Dictionary of KDE data.
    name (str): The name of the universe/system for output.
    """
    with open(f"{name}_upperclusterdata.pkl", 'wb') as file:
        pickle.dump(clustering_dic, file)
    
    with open(f"{name}_kde_data.pkl", 'wb') as file:
        pickle.dump(kde_dic, file)

def main():
    parser = argparse.ArgumentParser(description="Analyze MD simulation data.")
    parser.add_argument('--psf_files', nargs='+', help='List of PSF files', required=True)
    parser.add_argument('--dcd_files', nargs='+', help='List of DCD files', required=True)
    parser.add_argument('--names', nargs='+', help='Names for the systems', required=True)
    parser.add_argument('--lipids', nargs='+', help='List of lipid types', required=True)
    parser.add_argument('--grid_length', type=int, help='Grid length for KDE computation', default=100)

    args = parser.parse_args()
    
    if len(args.psf_files) != len(args.dcd_files) or len(args.psf_files) != len(args.names):
        print("Error: The number of PSF files, DCD files, and names must be the same.")
        return

    universes = load_universes(args.psf_files, args.dcd_files)
    names = dict(zip(universes, args.names))
    
    for universe in universes:
        analyze_universe(universe, args.lipids, names[universe], args.grid_length)

if __name__ == "__main__":
    main()
