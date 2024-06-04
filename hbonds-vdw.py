import argparse
import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from MDAnalysis.analysis import contacts
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import pickle

def contact_time_series(u, ag1, ag2, radius, lipid, step=1):
    """
    Calculate contacts over time between two atom groups within a given radius.

    Parameters:
    u (MDAnalysis.Universe): The MDAnalysis universe.
    ag1 (MDAnalysis.AtomGroup): The first atom group.
    ag2 (MDAnalysis.AtomGroup): The second atom group.
    radius (float): The cutoff radius for contacts.
    lipid (str): The lipid type.
    step (int): Step size for trajectory iteration.

    Returns:
    np.ndarray: Array of contact counts per frame.
    """
    contact_array = []
    for ts in tqdm(u.trajectory[::step], desc=f"{lipid} Contacts Progress"):
        dis = contacts.distance_array(ag1.positions, ag2.positions)
        ncon = contacts.contact_matrix(dis, radius).sum()
        contact_array.append(ncon)
    return np.array(contact_array)

def analyze_model(universe, lipids, model_name, step=1):
    """
    Analyze van der Waals contacts and hydrogen bonds between protein and lipids in a given universe.

    Parameters:
    universe (MDAnalysis.Universe): The MDAnalysis universe.
    lipids (list): List of lipid types.
    model_name (str): The name of the model for output.
    step (int): Step size for trajectory iteration.
    """
    vdw_mean = {}
    vdw_standard_deviation = {}
    vdw_contacts_dictionary = {}
    hbonds_bytime_dictionary = {}

    for lip in lipids:
        #selects hydrophobic residues for VDW contacts
        vdw_prot = universe.select_atoms("(protein and resname ALA LEU ILE PHE VAL PRO MET and not backbone and name C*) or resname GLYM and name C* and not name C N O CA H* S")
        lipid_carbons = universe.select_atoms(f"(resname {lip} and name C1 C2* C3*)")
        
        # Calculate van der Waals contacts
        contact_array = contact_time_series(universe, vdw_prot, lipid_carbons, 4, lip, step)
        vdw_contacts_dictionary[lip] = contact_array
        contact_mean = np.mean(contact_array)
        contact_std = np.std(contact_array)
        vdw_mean[lip] = contact_mean
        vdw_standard_deviation[lip] = contact_std

        # Calculate hydrogen bonds
        hbonds = HBA(universe=universe)
        #replace with desired protein selection
        hbonds.hydrogens_sel = hbonds.guess_hydrogens("protein or resname GLYM CYSF CYSG")
        #replace with desired lipid selection around protein
        hbonds.acceptors_sel = hbonds.guess_acceptors(f"resname {lip} and around 20 (protein or resname GLYM CYSF CYSG)")
        hbonds.run(step=step, verbose=True)
        by_time = hbonds.count_by_time()
        hbonds_bytime_dictionary[lip] = by_time

    save_results(vdw_contacts_dictionary, hbonds_bytime_dictionary, model_name)

def save_results(vdw_contacts_dictionary, hbonds_bytime_dictionary, model_name):
    """
    Save the VDW contacts and hydrogen bonds results to text files.

    Parameters:
    vdw_contacts_dictionary (dict): Dictionary of VDW contacts.
    hbonds_bytime_dictionary (dict): Dictionary of hydrogen bonds by time.
    model_name (str): The name of the model for output.
    """
    for lip, value in vdw_contacts_dictionary.items():
        np.savetxt(f"{model_name}_{lip}_vdwcontacts.txt", value)
    
    for lip, value in hbonds_bytime_dictionary.items():
        np.savetxt(f"{model_name}_{lip}_hbonds_bytime.txt", value)

def main():
    parser = argparse.ArgumentParser(description="Analyze MD simulation data for VDW contacts and hydrogen bonds.")
    parser.add_argument('--psf_files', nargs='+', help='List of PSF files', required=True)
    parser.add_argument('--dcd_files', nargs='+', help='List of DCD files', required=True)
    parser.add_argument('--model_names', nargs='+', help='Names for the models', required=True)
    parser.add_argument('--lipids', nargs='+', help='List of lipid types', required=True)
    parser.add_argument('--step', type=int, help='Step size for trajectory iteration', default=1)

    args = parser.parse_args()

    if len(args.psf_files) != len(args.dcd_files) or len(args.psf_files) != len(args.model_names):
        print("Error: The number of PSF files, DCD files, and model names must be the same.")
        return

    for psf_file, dcd_file, model_name in zip(args.psf_files, args.dcd_files, args.model_names):
        universe = mda.Universe(psf_file, dcd_file)
        analyze_model(universe, args.lipids, model_name, args.step)

if __name__ == "__main__":
    main()
