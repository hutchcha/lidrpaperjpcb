import argparse
import mdtraj as md
import numpy as np

# Mapping from secondary structure characters to numeric values
ss_to_numeric = {
    'H': 0,  # Alpha helix
    'B': 1,  # Beta bridge
    'E': 2,  # Beta strand
    'G': 3,  # 3-10 helix
    'I': 4,  # Pi helix
    'T': 5,  # Turn
    'S': 6,  # Bend
    'C': 7,  # Coil
    ' ': 8   # Unstructured
}

def ss_numeric(char):
    return ss_to_numeric.get(char)

vectorized = np.vectorize(ss_numeric)

def calculate_secondary_structure(top_file, traj_file, output_file):
    """
    Calculate the secondary structure of a protein and save the results as numeric values.

    Parameters:
    top_file (str): Path to the topology file.
    traj_file (str): Path to the trajectory file.
    output_file (str): Path to save the output secondary structure file.
    """
    trajectory = md.load(traj_file, top=top_file)
    dssp = md.compute_dssp(trajectory, simplified=False)
    dssp_numeric = vectorized(dssp)
    np.savetxt(output_file, dssp_numeric)

def main():
    parser = argparse.ArgumentParser(description="Calculate secondary structure from a molecular dynamics trajectory.")
    parser.add_argument('--top_file', type=str, help='Path to the topology file', required=True)
    parser.add_argument('--traj_file', type=str, help='Path to the trajectory file', required=True)
    parser.add_argument('--output_file', type=str, help='Path to save the output secondary structure file', required=True)

    args = parser.parse_args()
    calculate_secondary_structure(args.top_file, args.traj_file, args.output_file)

if __name__ == "__main__":
    main()
