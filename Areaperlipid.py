import argparse
import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.leaflet import AssignLeaflets
from MDAnalysis.analysis.hole import AreaPerLipid
from tqdm import tqdm

def calculate_area_per_lipid(psf_file, dcd_file, lipid_sel, output_file, step=1):
    """
    Calculate the average area per lipid for selected lipids in a molecular dynamics simulation.

    Parameters:
    psf_file (str): Path to the PSF file.
    dcd_file (str): Path to the DCD file.
    lipid_sel (str): Lipid selection string.
    output_file (str): Path to save the output area file.
    step (int): Step size for trajectory iteration.
    """
    universe = mda.Universe(psf_file, dcd_file)
    
    # Assign leaflets
    leaf = AssignLeaflets(
        universe=universe,
        lipid_sel=lipid_sel
    )
    leaf.run(
        start=None,
        stop=None,
        step=step,
        verbose=True
    )

    # Calculate area per lipid
    areas = AreaPerLipid(
        universe=universe,
        lipid_sel=lipid_sel,
        leaflets=leaf.leaflets
    )
    areas.run(
        start=None,
        stop=None,
        step=step,
        verbose=True
    )

    # Calculate average area per lipid
    average_area = areas.areas.mean(axis=0)
    np.savetxt(output_file, average_area)

def main():
    parser = argparse.ArgumentParser(description="Calculate average area per lipid in a molecular dynamics simulation.")
    parser.add_argument('--psf_file', type=str, help='Path to the PSF file', required=True)
    parser.add_argument('--dcd_file', type=str, help='Path to the DCD file', required=True)
    parser.add_argument('--lipid_sel', type=str, help='Lipid selection string', required=True)
    parser.add_argument('--output_file', type=str, help='Path to save the output area file', required=True)
    parser.add_argument('--step', type=int, help='Step size for trajectory iteration', default=1)

    args = parser.parse_args()
    calculate_area_per_lipid(args.psf_file, args.dcd_file, args.lipid_sel, args.output_file, args.step)

if __name__ == "__main__":
    main()
