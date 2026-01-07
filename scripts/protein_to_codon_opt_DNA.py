"""
This script, given a csv file that contains protein sequences under the column 'seq', returns two multifasta files one with the protein sequences and one with the codon optimised 
DNA sequence, using IDT's API

"""
### MODULES #####################################################################################
import requests
import pandas as pd
from Bio.Data import CodonTable
from pathlib import Path
import textwrap
import argparse
import math
import random
### PARAMS ######################################################################################

id_col_name = "Protein Name" # Indicate the name of the column in your csv that contains the protein IDs
seq_col_name = "Sequence" # indicate the name of the column in your csv that contains the protein sequence 

organism = "Escherichia coli K12" # Select one of the organisms that appear in IDT's users interface
avoid_restriction_sites = "BsaI" # Include the name of the enzyme to avoid. Must be in IDT's interface aswell

remove_aa = [('R','K')] # None or empty list if no changes. If you wanna replace ARG with LYS for example: (R,K)
limit_cysteine = 1 # float or None. If float, in sequences with more than x cysteines, cysteines will be replaced with serine until the limit is reached.

limit_net_charge = -3.0 # Protein sequence net charge.None or float value. If float is given, sequences with net charge above this value will edited to reduce net charge 
pH = 7.4 # pH at which net charge is calculated

terminal_T= False # True for making the last nucleotide of the sequence to end with T

include_stop_codon= False # Ensures that the sequences in the fasta will start with ATG. IMPORTANT! If false, stop codons will be removed if they are there
include_start_codon= True # Ensures that the sequences in the fasta will finish with a stop codon. Default TAA
 
DNA_prefix = "GGCTACGGTCTCGAGGA" # The string of nucleotides indicated here will be appended to the beginning of every coding sequence. i.e a restriction site
DNA_suffix = "TTCCTGAGACCCATCAT" # The string of nucleotides indicated here will be appended to the end of every coding sequence. i.e a restriction site
min_DNA_length = 300 # If the final DNA sequence is shorter than this value, random nucleotides will be added to the beginning of the sequence until the length is reached.

fasta_single_line = True # FASTA output will be single line, overrides any line length
fasta_line_length = 60 # FASTA output will be in lines of this length
### FUNCTIONS ###################################################################################

# Function to authenticate and get the bearer token
def get_bearer_token(client_id, client_secret, username, password):
    url = "https://eu.idtdna.com/Identityserver/connect/token"
    headers = {
        "Content-Type": "application/x-www-form-urlencoded"
    }
    data = {
        "grant_type": "password",
        "username": username,
        "password": password,
        "scope": "test",
        "client_id": client_id,
        "client_secret": client_secret
    }

    response = requests.post(url, headers=headers, data=data)

    if response.status_code == 200:
        return response.json()['access_token']
    else:
        print(f"Failed to authenticate: {response.status_code}, {response.text}")
        return None

# Function to optimize codons
def codoonopt_call(bearer_token, sequence, name="MyAminoAcidSequence"):
    url = "https://eu.idtdna.com/restapi/v1/CodonOpt/Optimize"
    headers = {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {bearer_token}"
    }
    payload = {
        "Organism": organism,
        "optimizationItems": [
            {
                "Name": name,
                "Sequence": sequence
            }
        ],
        "sequenceType": "Amino Acids",
        "productType": "gblock",
    }

    response = requests.post(url, headers=headers, json=payload)

    if response.status_code == 200:
        return response.json()
    else:
        print(f"Failed to optimize codons: {response.status_code}, {response.text}")
        return None

def handle_restriction_sites_output(restriction_site_str):
    """The str typically looks like this:
    'The following restriction enzyme sites have been found in the selected reading frame:\nEcoRI (GAATTC)\nEcoRV
    (GATATC)\nFspI (TGCGCA)\nMfeI (CAATTG)\nMlyI (GAGTC)'.
    """
    restrict_enz = [i.split()[0] for i in restriction_site_str.split("\n")[1:]]
    restrict_recog_site = [i.split()[1] for i in restriction_site_str.split("\n")[1:]]
    return restrict_enz, restrict_recog_site

# Main function to execute the script
def optimize_codon_sequence(aa_sequence:str, avoid_restriction_enzyme=None, max_attempts=100):
    client_id = "madjepidt"
    client_secret = "43ea15e0-2686-45fb-a54c-5ff6c8b57f6a"
    username = "madsjeppesen"
    password = "sypwureg14"

    # Get bearer token
    token = get_bearer_token(client_id, client_secret, username, password)

    if token:
        # Perform codon optimization
        optimized_result = codoonopt_call(token, aa_sequence)

        if optimized_result:
            optresults = optimized_result[0]["OptResult"]
            restrict_enz, restrict_recog_site = handle_restriction_sites_output(optresults["RestrictionSites"])
            full_sequence = optresults["FullSequence"]
            complexity_score = optresults["ComplexityScore"]
            # if need to optimize
            attempt = 1
            if avoid_restriction_enzyme is not None:
                while avoid_restriction_enzyme in restrict_enz and attempt <= max_attempts:
                    print(f"Found {avoid_restriction_enzyme} in the codon optimized sequence. "
                          f"Will attempt to optimize again (Attempts = {attempt}/{max_attempts})")

                    optimized_result = codoonopt_call(token, aa_sequence)
                    optresults = optimized_result[0]["OptResult"]
                    restrict_enz, restrict_recog_site = handle_restriction_sites_output(optresults["RestrictionSites"])
                    full_sequence = optresults["FullSequence"]
                    complexity_score = optresults["ComplexityScore"]
                    attempt += 1
                if avoid_restriction_enzyme in restrict_enz:
                    print(f"Could not create a codon sequence void of {avoid_restriction_enzyme} after {max_attempts}.")
                    return None
            # 'The following restriction enzyme sites have been found in the selected reading frame:\nAfeI (AGCGCT)'
            return {"sequence": full_sequence,
                    "restriction_enzymes": restrict_enz,
                    "restriction_recognition_sites": restrict_recog_site,
                    "complexity_score": complexity_score}
        else:
            print("Optimization failed.")
    else:
        print("Authentication failed.")

def get_codon_triplets_from_DNA(seq):
    codons = [seq[i:i + 3] for i in range(0, len(seq), 3)]
    assert all(len(c) == 3 for c in codons)
    return codons

def translate_DNA_to_protein(seq, end_with_stop_codon = True):
    ct = CodonTable.unambiguous_dna_by_name["Standard"]
    protein_seq = ""
    for codon in get_codon_triplets_from_DNA(seq):
        if end_with_stop_codon and codon in ct.stop_codons:
            print(f"Found Stop Codon: {codon}")
            return protein_seq
        else:
            protein_seq += ct.forward_table[codon]
    return protein_seq

def calculate_net_charge(sequence, pH=7.4):
    """
    Calculates the net charge of a protein sequence at a specific pH 
    using the Henderson-Hasselbalch equation and standard pKa values (EMBOSS scale).
    """
    # Standard pKa values (EMBOSS scale)
    pKa_values = {
        # Positively charged (Basic)
        'K': 10.8, # Lysine
        'R': 12.5, # Arginine
        'H': 6.5,  # Histidine
        # Negatively charged (Acidic)
        'D': 3.9,  # Aspartic Acid
        'E': 4.1,  # Glutamic Acid
        'C': 8.5,  # Cysteine (Weak)
        'Y': 10.1, # Tyrosine (Weak)
        # Termini
        'N_TERM': 8.6, # Amino terminus
        'C_TERM': 3.6  # Carboxyl terminus
    }

    seq = sequence.upper()
    net_charge = 0.0

    # 1. Calculate charge of the N-terminus (Positive)
    net_charge += 1.0 / (1.0 + math.pow(10, (pH - pKa_values['N_TERM'])))

    # 2. Calculate charge of the C-terminus (Negative)
    net_charge -= 1.0 / (1.0 + math.pow(10, (pKa_values['C_TERM'] - pH)))

    # 3. Calculate internal residues
    for aa in seq:
        if aa in pKa_values:
            pka = pKa_values[aa]
            
            # Basic amino acids (Positive contribution: K, R, H)
            if aa in ['K', 'R', 'H']:
                net_charge += 1.0 / (1.0 + math.pow(10, (pH - pka)))
            
            # Acidic amino acids (Negative contribution: D, E, C, Y)
            else:
                net_charge -= 1.0 / (1.0 + math.pow(10, (pka - pH)))

    return net_charge

def adjust_net_charge(seq, target_charge):
    """
    Iteratively replaces positively charged amino acids (K, R, H) with 
    Glutamic Acid (E) until the net charge is <= target_charge.
    """
    # 1. Check if adjustment is even needed
    current_charge = calculate_net_charge(seq)
    if current_charge <= target_charge:
        return seq, 0 # Return original sequence and 0 replacements

    # 2. Identify positions of positive residues (K, R, H)
    # We store indices to replace them specifically
    positive_aa_positions = [i for i, aa in enumerate(seq) if aa in ['K', 'R', 'H']]
    
    # 3. Shuffle to randomize which residues get replaced first
    # This avoids clustering mutations at the N-terminus
    random.shuffle(positive_aa_positions)

    # 4. Convert to list for mutable editing
    seq_list = list(seq)
    replacements_count = 0

    # 5. Iterate and mutate
    for pos in positive_aa_positions:
        # Break if we have reached the target
        if current_charge <= target_charge:
            break
        
        # Perform replacement: Positive AA -> Glutamic Acid (E)
        # This causes a large charge drop (approx -2 per change at pH 7.4)
        seq_list[pos] = 'E'
        replacements_count += 1
        
        # Re-calculate charge
        current_seq = "".join(seq_list)
        current_charge = calculate_net_charge(current_seq)

    final_seq = "".join(seq_list)

    return final_seq, current_charge

def generate_multifasta(df, id_column, seq_column, output, line_width=60, single_line=False):
    """
    Generates a multi-FASTA file from a pandas DataFrame.

    Args:
        df (pd.DataFrame): Input dataframe containing sequences.
        id_column (str): Name of the column to use as the sequence header/ID.
        seq_column (str): Name of the column containing the DNA or Protein sequence.
        output (str): File path for the output FASTA file.
        line_width (int, optional): Number of characters per line for sequence wrapping. 
                                    Defaults to 60. Ignored if single_line is True.
        single_line (bool, optional): If True, writes the entire sequence on one line.
                                      Defaults to False.
    """
    # Basic validation
    if id_column not in df.columns:
        raise ValueError(f"Column '{id_column}' not found in DataFrame.")
    if seq_column not in df.columns:
        raise ValueError(f"Column '{seq_column}' not found in DataFrame.")

    try:
        with open(output, 'w') as f:
            for index, row in df.iterrows():
                # Extract ID and Sequence
                seq_id = str(row[id_column])
                sequence = str(row[seq_column])

                # Write Header
                f.write(f">{seq_id}\n")

                # Write Sequence
                if single_line:
                    f.write(sequence + "\n")
                elif line_width:
                    # Wrap sequence into fixed-width lines
                    wrapped_seq = textwrap.fill(sequence, width=line_width)
                    f.write(wrapped_seq + "\n")
                else:
                    # Fallback if line_width is 0 or None
                    f.write(sequence + "\n")
        
        print(f"Successfully generated {output} with {len(df)} sequences.")

    except IOError as e:
        print(f"Error writing to file {output}: {e}")

def process_aa_sequence(seq, remove_aa, name):
    # 1. Safety: Ensure seq is a string and remove hidden whitespace
    seq = str(seq).strip()
    
    total_replacements = 0  # Initialize accumulator
    
    # 2. Iterate through replacements
    if remove_aa:
        for aa_from, aa_to in remove_aa:
            # Count current occurrence
            count = seq.count(aa_from)
            
            if count > 0:
                # Accumulate the total count (Fixes the bug)
                total_replacements += count 
                print(f"⚠️  REPLACEMENT WARNING: Replaced {count} instances of '{aa_from}' with '{aa_to}' in {name}.")
                
                # Perform replacement
                seq = seq.replace(aa_from, aa_to)
    
    # 3. Compute net charge 
    charge = calculate_net_charge(seq)
    
    if limit_net_charge is not None and charge > limit_net_charge:
        print(f"⚠️  CHARGE WARNING: Net charge of {name} is {charge}, exceeding limit of {limit_net_charge}. Adjusting sequence.")
        print(f"Original Sequence: {seq}")
        seq, charge = adjust_net_charge(seq, limit_net_charge)
    
    # 4. Cysteine check
    if limit_cysteine is not None:
        cysteine_count = seq.count('C')
        if cysteine_count > limit_cysteine:
            print(f"⚠️  CYSTEINE WARNING: {cysteine_count} cys in {name}, exceeding limit of {limit_cysteine}. Replacing C with S.")
            while cysteine_count > limit_cysteine:
                # Replace first occurrence of C with S
                seq = seq.replace('C', 'S', 1)
                cysteine_count -= 1
                total_replacements += 1
            
    cysteine_count = seq.count('C')
    # 4. Return the list
    return [seq, total_replacements, charge]

def mutate_to_T(seq,name):
    
    # 1. If the last letter is T, return unchanged
    if seq.endswith('T'):
        return seq

    # 2. Extract the last codon and the prefix
    if len(seq) < 3:
        return seq # Safety for short sequences
        
    old_codon = seq[-3:]
    prefix_seq = seq[:-3]

    # Data: Genetic Code (Codon -> AA)
    genetic_code = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    
    # Reverse Map (AA -> List of Codons ending in T)
    # We pre-calculate this to make the function faster
    aa_to_T_codons = {
        'I': ['ATT'], 'T': ['ACT'], 'N': ['AAT'], 'S': ['AGT', 'TCT'],
        'R': ['CGT'], 'L': ['CTT'], 'P': ['CCT'], 'H': ['CAT'],
        'V': ['GTT'], 'A': ['GCT'], 'D': ['GAT'], 'G': ['GGT'],
        'F': ['TTT'], 'Y': ['TAT'], 'C': ['TGT'], '_': ['TAA'] # Stop
    }

    # Similarity Map (If AA has no T-codon, swap to this AA)
    # Logic: Charge > Polarity > Structure
    similarity_substitutions = {
        'M': 'I', # Methionine -> Isoleucine (Both Hydrophobic, maintain structure)
        'K': 'R', # Lysine -> Arginine (Both Positively Charged)
        'Q': 'H', # Glutamine -> Histidine (Both Polar/Positively charged context)
        'E': 'D', # Glutamic Acid -> Aspartic Acid (Both Negatively Charged)
        'W': 'F', # Tryptophan -> Phenylalanine (Both Aromatic/Hydrophobic)
    }

    # 3. Find the Amino Acid
    current_aa = genetic_code.get(old_codon)
    
    if not current_aa:
        return seq + " [Error: Invalid Codon]"

    # 4. Find a codon for this AA that ends in T
    if current_aa in aa_to_T_codons:
        # We found a synonymous codon ending in T!
        new_codon = aa_to_T_codons[current_aa][0] # Pick the first available
        return prefix_seq + new_codon

    # 5. If no T-ending codon exists, swap Amino Acid
    else:
        # Find the best substitute
        alt_aa = similarity_substitutions.get(current_aa)
        
        if alt_aa:
            new_codon = aa_to_T_codons[alt_aa][0]
            print(f"⚠️  WARNING: {name} Mutated Amino Acid from {current_aa} to {alt_aa} to satisfy 'T' ending.")
            return prefix_seq + new_codon
        else:
            print(f"⚠️  WARNING: {name} Could not find a suitable substitution for {current_aa}.")
            return seq

def process_dna_sequence(seq,name):
    seq = seq.upper()  # Ensure consistency

    # Convert last nucleotide to T if needed
    if terminal_T:
        seq = mutate_to_T(seq,name)
    
    # Ensure there is an starting codon
    if include_start_codon:
        if not seq.startswith("ATG"):
            seq = "ATG" + seq
    # Ensure there is a finish codon
    if include_stop_codon:
        if not seq.endswith(("TAA","TAG","TGA")):
            seq = seq + "TAA"
    else:
        if seq.endswith(("TAA","TAG","TGA")):
            seq = seq[:-3]
    # Add suffix and prefix
    seq = DNA_prefix + seq

    # Perform C termini check and addition
    # QUICK:
    # To match the readingframe the T of the c_term flanking sequence has to be either merged into the amino acid
    # if the amino acid has a T at the 3rd position or else glycine is added
    # DETAILS:
    # The c termninus flanking sequence is 17 bp and therefor is not divisible by 3 (16 or 19 is). Therefore we either have to
    # remove 1 bp to get 16 or add 2 to get 19 so that the reading frame fits. if removing 1 bp, we can merge it into
    # the flanking sequence. The flanking sequence starts with T so the c terminus codon needs to end with T as well.
    # If it does not we just add a GLY in there. GLY codon  is GG[A/T/G/C] where in [] anything can be used. Since the
    # flanking sequence starts with T we just add GG, so it become GGT. Also GGT is used 48% of the time in E coli:
    # https://bmcbiotechnol.biomedcentral.com/articles/10.1186/1472-6750-9-36. SEE TABLE 3.

    if DNA_suffix:
        if sequence[-1] == "T":
            sequence = sequence[:-1]
            sequence += DNA_suffix
        else:
            sequence += "GG" + DNA_suffix
    # Add random padding if needed
    if min_DNA_length is not None:
        while len(seq) < min_DNA_length:
            random_nt = random.choice(['A', 'T', 'C', 'G'])
            seq = random_nt + seq 
    return seq

def generate_96_well_map():
    """Generates a list ['A1', 'A2', ... 'H12']"""
    rows = list("ABCDEFGH")
    wells = []
    for r in rows:
        for c in range(1, 13):
            wells.append(f"{r}{c}")
    return wells

    
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

def generate_DNA_cds_multifasta(sequence_df, general_output_path):
    # Initialise DNA df
    df_opt_codon_global = pd.DataFrame()

    # Iterate the sequence df filling the DNA df
    for _, row in sequence_df.iterrows():
        name, AA_seq = row[id_col_name], row[seq_col_name]
        print(f"> {name}")
        print(f"Original:{AA_seq}")
        aa_result = process_aa_sequence(AA_seq, remove_aa, name)
        AA_seq, replacements, charge = aa_result[0], aa_result[1], aa_result[2]
        
        optresults = optimize_codon_sequence(AA_seq, avoid_restriction_enzyme=avoid_restriction_sites)
        
        # Check if optimization was successful before proceeding
        if optresults is None:
            print(f"⚠️ Skipping {name}: No optimized codon sequence found.")
            continue 
            
        DNA_seq = optresults["sequence"]

        # Translation check
        translated_prot = str(translate_DNA_to_protein(DNA_seq))
        if translated_prot != str(AA_seq):
            print(f"⚠️ MISMATCH WARNING: Translation mismatch for {name}. Processing continues...")
            print(f"Translated: {translated_prot}")
            print(f"Input: {AA_seq}")
        
        print(f"Final:{AA_seq}")
        DNA_seq = process_dna_sequence(DNA_seq, name)

        # Prepare restriction enzyme strings
        enzyme_string = "|".join(optresults["restriction_enzymes"])
        recognition_sites_string = "|".join(optresults["restriction_recognition_sites"])
        
        # Generate the individual dictionary
        df_opt_codon_indiv = {
            "name": name, 
            "aa_sequence": AA_seq, 
            "codon_sequence": DNA_seq,
            "restriction_enzymes": enzyme_string,
            "restriction_recognition_sites": recognition_sites_string,
            "complexity_score": optresults["complexity_score"],
            "NetCharge": charge,
            "fasta_header": f"{name}|CDS|{optresults['complexity_score']}|{enzyme_string}|{recognition_sites_string}",
            "protein_header": f"{name}|PROT|{len(AA_seq)} aa |{replacements} AA_replacements | NetCharge:{charge}"
        }
        
        # Convert dict to DataFrame and append
        df_opt_codon_indiv = pd.DataFrame([df_opt_codon_indiv])
        df_opt_codon_global = pd.concat([df_opt_codon_global, df_opt_codon_indiv], ignore_index=True)
    
    # Generate the protein multifasta
    generate_multifasta(df_opt_codon_global, "protein_header", "aa_sequence", f"{general_output_path}/protein_sequences.fasta", single_line=fasta_single_line, line_width=fasta_line_length)
    
    # Generate the DNA multifasta
    generate_multifasta(df_opt_codon_global, "fasta_header", "codon_sequence", f"{general_output_path}/codon_optimized_DNA_sequences.fasta", single_line=fasta_single_line, line_width=fasta_line_length)
    
    # Save global file as csv
    df_opt_codon_global.to_csv(f"{general_output_path}/codon_optimized_DNA_sequences.csv", index=False)
    
    return

def generate_DNA_insert_multifasta():
    return
### INPUT CHECK #################################################################################
parser = argparse.ArgumentParser(description="Generate codon optimized DNA sequences from protein sequences using IDT's API.")
parser.add_argument("--input_csv", type=str, help="Path to the input CSV file containing protein sequences.")
parser.add_argument("--output_dir", type=str, help="Directory to save the output multifasta files.")
args = parser.parse_args()

input_csv = args.input_csv
output_dir = args.output_dir
Path(output_dir).mkdir(parents=True, exist_ok=True)
### MAIN EXECUTION ##############################################################################
generate_DNA_cds_multifasta(pd.read_csv(input_csv),output_dir)