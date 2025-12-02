"""

"""
### MODULES ###########################################################################################

import os
import csv
import warnings
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import argparse
from pathlib import Path

### FUNCTIONS #########################################################################################
def get_sequence_from_pdb(file_path, parser=None, ppb=None):
    """
    Extracts the amino acid sequence from a single PDB file.
    Returns the sequence string (chains separated by /).
    """
    # Initialize tools if not provided (allows standalone use)
    if parser is None:
        parser = PDBParser(QUIET=True)
    if ppb is None:
        ppb = PPBuilder()

    # Use filename as structure ID
    structure_id = os.path.basename(file_path)
    
    # Parse the structure
    structure = parser.get_structure(structure_id, file_path)
    
    # Extract Sequence
    # PPBuilder gets the sequence from the ATOM coordinates.
    peptides = ppb.build_peptides(structure)
    
    # If there are multiple chains (e.g., Chain A and Chain B), 
    # we join them with a slash '/' to keep them distinct but in one cell.
    sequence_list = [str(pp.get_sequence()) for pp in peptides]
    return "/".join(sequence_list)

def pdb_folder_to_csv(input_folder, output_csv_name):
    """
    Reads all .pdb files in a folder, extracts sequences, and writes to CSV.
    """
    
    # Initialize PDB parser once for efficiency
    parser = PDBParser(QUIET=True)
    ppb = PPBuilder()
    
    # Check if folder exists
    if not os.path.exists(input_folder):
        print(f"Error: The folder '{input_folder}' does not exist.")
        return

    # Prepare list to store data
    data_rows = []
    
    print(f"Scanning folder: {input_folder}...")

    # Iterate through files
    for filename in os.listdir(input_folder):
        if filename.endswith(".pdb"):
            
            # 1. Extract Protein Name (Split by '.' and take first element)
            protein_name = filename.split('.')[0]
            file_path = os.path.join(input_folder, filename)
            
            try:
                # 2. Get Sequence using helper function
                full_sequence = get_sequence_from_pdb(file_path, parser, ppb)
                
                if full_sequence:
                    data_rows.append([protein_name, full_sequence])
                else:
                    data_rows.append([protein_name, "No Protein Sequence Found"])
                    
            except Exception as e:
                print(f"Error processing {filename}: {e}")
                data_rows.append([protein_name, "Error Reading File"])

    # 4. Write to CSV
    with open(output_csv_name, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Header
        writer.writerow(["Protein Name", "Sequence"])
        # Data
        writer.writerows(data_rows)

    print(f"Done! Successfully processed {len(data_rows)} PDB files.")
    print(f"Results saved to: {output_csv_name}")

### INPUT CHECK ####################################################################################

parser = argparse.ArgumentParser(description="Generate codon optimized DNA sequences from protein sequences using IDT's API.")
parser.add_argument("--input_folder", type=str, help="Path to the input CSV file containing protein sequences.")
parser.add_argument("--output_file", type=str, help="Directory to save the output multifasta files.")
args = parser.parse_args()

input_folder = args.input_folder
output_csv_name = args.output_file
#Path(output_csv_name).mkdir(parents=True, exist_ok=True)

### MAIN EXECUTION #################################################################################
pdb_folder_to_csv(input_folder, output_csv_name)