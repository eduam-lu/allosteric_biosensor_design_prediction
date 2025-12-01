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

### PARAMS ######################################################################################

remove_aa = remove_aa = [('R','K')] # None or empty list if no changes. If you wanna replace ARG with LYS for example: [R,K]

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
        "Organism": "Escherichia coli K12",
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

def process_aa_sequence(seq,remove_aa):
    # Replace aa if needed
    if remove_aa:
        for aa_from, aa_to in remove_aa:
            seq = seq.replace(aa_from, aa_to)
    return seq

def process_dna_sequence(seq):
    # Ensure there is an starting codon
    if not seq.startswith("ATG"):
        seq = "ATG" + seq
    # Ensure there is a finish codon
    if not seq.endswith(("TAA","TAG","TGA")):
        seq = seq + "TAA"
    # If insert mode, apply insert treatment
    # under development
    return seq

def generate_DNA_cds_multifasta(sequence_df,general_output_path):

    # Initialise DNA df
    df_opt_codon_global = pd.DataFrame()

    # Iterate the sequence df filling the DNA df
    for _, row in sequence_df.iterrows():
        name, AA_seq = row["fileID"], row["seq"]
        AA_seq = process_aa_sequence(AA_seq,remove_aa)
        optresults = optimize_codon_sequence(AA_seq, avoid_restriction_enzyme="BsaI")
        DNA_seq = optresults["sequence"]
        #assert str(translate_DNA_to_protein(DNA_seq)) == str(AA_seq)
        DNA_seq = process_dna_sequence(DNA_seq)
        if optresults is None:
            raise ValueError(f"No optimized codon sequence found for {name}")
        else:
            # Prepare restriction enzyme string
            enzyme_string = ",".join(optresults["restriction_enzymes"])
            recognition_sites_string = ",".join(optresults["restriction_recognition_sites"])
            # Generate the novel dataframe
            df_opt_codon_indiv = {"name": name, "aa_sequence": AA_seq, "codon_sequence": DNA_seq,
                            "restriction_enzymes": enzyme_string,
                            "restriction_recognition_sites": recognition_sites_string,
                            "complexity_score": optresults["complexity_score"],
                            "fasta_header": f"{name}|CDS|{optresults['complexity_score']}|{enzyme_string}|{recognition_sites_string}"}
            df_opt_codon_indiv = pd.DataFrame(df_opt_codon_indiv, index=[0])
        df_opt_codon_global = pd.concat([df_opt_codon_global,df_opt_codon_indiv])
    
    # Generate the protein multifasta
    generate_multifasta(df_opt_codon_global,"name", "aa_sequence",f"{general_output_path}/protein_sequences.fasta")
    # Generate the DNA multifasta
    generate_multifasta(df_opt_codon_global,"fasta_header", "codon_sequence",f"{general_output_path}/codon_optimized_DNA_sequences.fasta")
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