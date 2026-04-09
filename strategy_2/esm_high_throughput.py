import torch
import esm
import pandas as pd
import argparse
import os
import sys
from pathlib import Path
from tqdm import tqdm  
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(f"[gpu-pin] esm_high_throughput visible={os.environ.get('CUDA_VISIBLE_DEVICES', 'ALL')} device={device}")

# Load model
model = esm.pretrained.esmfold_v1()
model = model.eval().to(device)

# Input check
parser = argparse.ArgumentParser(
    description="Given an input csv, it predicts all the sequences in it with ESMfold",
)
parser.add_argument('--input_csv', help="Path to the input structure file (PDB format)", type=str, required=True)
parser.add_argument('--output_folder', help="Output folder for PDBs", type=str, required=True)

args = parser.parse_args()

# Check if the input csv exists
if not Path(args.input_csv).exists():
    print(f'File: {args.input_csv} could not be found')
    sys.exit(1)

# Create the output folder if it doesn't exist
# Note: Fixed bug here (changed args.output to args.output_folder)
Path(args.output_folder).mkdir(parents=True, exist_ok=True)

# Load MPNN_csv as a df
MPNN_df = pd.read_csv(args.input_csv)

# Parse data frame
# <--- 2. Wrap the iterator with tqdm and provide the total count
for i, row in tqdm(MPNN_df.iterrows(), total=len(MPNN_df), desc="Predicting Structures"):
    sequence = row['seq']
    file_name = row['file_ID']
    
    with torch.no_grad():
        output = model.infer_pdb(sequence)

    with open(f'{args.output_folder}/{file_name}.pdb', "w") as f:
        f.write(output)

