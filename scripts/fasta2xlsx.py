import pandas as pd
import random
import sys
import os

def parse_fasta(fasta_file):
    """
    Parses a FASTA file manually to avoid BioPython dependency
    if you only need simple header/sequence extraction.
    """
    sequences = []
    current_header = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue

            if line.startswith(">"):
                if current_header:
                    sequences.append((current_header, "".join(current_seq)))
                current_header = line[1:] # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)
        
        # Append the last record
        if current_header:
            sequences.append((current_header, "".join(current_seq)))
    
    return sequences

def generate_96_well_map():
    """Generates a list ['A1', 'A2', ... 'H12']"""
    rows = list("ABCDEFGH")
    wells = []
    for r in rows:
        for c in range(1, 13):
            wells.append(f"{r}{c}")
    return wells

def process_sequence(sequence):
    nterm = "GGCTACGGTCTCGAGGA"
    cterm = "TTCCTGAGACCCATCAT"
    # Remove start codon if needed
    if sequence.startswith("ATG"):
        sequence = sequence[3:]
    # Add N termini
    sequence = nterm + sequence
    # Perform C termini check and addition
    if sequence[-1] == "T":
        sequence = sequence[:-1]
        sequence += cterm
    else:
        sequence += "GG" + cterm
    
    # Pad the beginning with random nucleotides if under 300 nts until 300 nts
    while len(sequence) < 300:
        random_nt = random.choice(['A', 'T', 'C', 'G'])
        sequence = random_nt + sequence  # You can randomize this if needed
    
    return sequence

def fasta_to_xlsx_plate(fasta_file, excel_file):
    # 1. Get Data
    seq_data = parse_fasta(fasta_file)
    well_map = generate_96_well_map()
    
    # 2. Map sequences to wells
    final_data = []
    
    for i, (header, sequence) in enumerate(seq_data):
        # Calculate well index (0-95). 
        # The % 96 ensures that if you have >96 sequences, it restarts at A1.
        well_index = i % 96
        well_position = well_map[well_index]

        sequence = process_sequence(sequence)

        final_data.append({
            "Well Position": well_position,
            "Header": header,
            "Sequence": sequence
        })

    # 3. Create DataFrame and Save
    df = pd.DataFrame(final_data)
    
    # Reorder columns explicitly to match your request
    df = df[["Well Position", "Header", "Sequence"]]
    
    try:
        df.to_excel(excel_file, index=False)
        print(f"Success! Processed {len(seq_data)} sequences.")
        print(f"Saved to: {excel_file}")
    except PermissionError:
        print(f"Error: Could not write to {excel_file}. Make sure the file is closed.")

# --- Usage ---
if __name__ == "__main__":
    input_fasta = "/home/eduardo/final_DNA_sequences.fasta"   # Change this to your file name
    output_xlsx = "plate_layout.xlsx"
    
    if os.path.exists(input_fasta):
        fasta_to_xlsx_plate(input_fasta, output_xlsx)
    else:
        print(f"File not found: {input_fasta}")