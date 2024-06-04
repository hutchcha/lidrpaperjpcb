import argparse
import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.analysis import rdf
from tqdm import tqdm

def calculate_rdf(psf_file, dcd_file, output_file, step=20):
    """
    Calculate the radial distribution function (RDF) for specified atom groups in a molecular dynamics simulation.

    Parameters:
    psf_file (str): Path to the PSF file.
    dcd_file (str): Path to the DCD file.
    output_file (str): Path to save the RDF results.
    step (int): Step size for trajectory iteration.
    """
    u = mda.Universe(psf_file, dcd_file)

    ng = ["LYS", "ARG", "ASP"]
    ng2 = ["POPC", "POPE", "SAPI", "CHL1"]

    rdf_results = []

    for n in ng:
        ag = u.select_atoms(f"resname {n} and name CA")

        ids = ag.resids
        num = len(ids)
        
        for i in ids:
            if num != 0:
                if n == "LYS":
                    ag1 = u.select_atoms(f"resname {n} and resid {i} and name NZ")
                    name = "Lysine"
                elif n == "ASP":
                    ag1 = u.select_atoms(f"resname {n} and resid {i} and name OD1 OD2")
                    name = "Aspartic Acid"
                else:
                    continue

                for l in ng2:
                    if l == "CHL1":
                        ag2 = u.select_atoms(f"resname {l} and name O3")
                    else:
                        ag2 = u.select_atoms(f"resname {l} and name O11 O12 O13 O14")

                    rd = rdf.InterRDF(ag1, ag2, step=step)
                    rd.run(verbose=True)

                    rdf_results.append({
                        "resname1": n,
                        "resid1": i,
                        "resname2": l,
                        "rdf": rd.results.rdf,
                        "bins": rd.results.bins
                    })

                    if n == "ASP" and l == "POPE":
                        ag2 = u.select_atoms(f"resname {l} and name N")
                        rd = rdf.InterRDF(ag1, ag2, step=step)
                        rd.run(verbose=True)

                        rdf_results.append({
                            "resname1": n,
                            "resid1": i,
                            "resname2": f"{l}-N",
                            "rdf": rd.results.rdf,
                            "bins": rd.results.bins
                        })

    # Save RDF results to a file
    with open(output_file, 'wb') as f:
        np.save(f, rdf_results)

def main():
    parser = argparse.ArgumentParser(description="Calculate RDF for specified atom groups in a molecular dynamics simulation.")
    parser.add_argument('--psf_file', type=str, help='Path to the PSF file', required=True)
    parser.add_argument('--dcd_file', type=str, help='Path to the DCD file', required=True)
    parser.add_argument('--output_file', type=str, help='Path to save the RDF results', required=True)
    parser.add_argument('--step', type=int, help='Step size for trajectory iteration', default=20)

    args = parser.parse_args()
    calculate_rdf(args.psf_file, args.dcd_file, args.output_file, args.step)

if __name__ == "__main__":
    main()
