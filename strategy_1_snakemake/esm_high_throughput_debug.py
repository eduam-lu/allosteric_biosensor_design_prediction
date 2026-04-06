import torch
import esm
import pandas as pd
import argparse
import os
import sys
from pathlib import Path
from tqdm import tqdm  
import ssl
import time
ssl._create_default_https_context = ssl._create_unverified_context

# DIAGNOSTIC: Capture environment and startup state
print("="*80)
print("ESM DIAGNOSTICS - JOB INITIALIZATION")
print("="*80)
print(f"PID: {os.getpid()}")
print(f"CUDA_VISIBLE_DEVICES: {os.environ.get('CUDA_VISIBLE_DEVICES', 'NOT SET')}")
print(f"TORCH_HOME: {os.environ.get('TORCH_HOME', 'NOT SET')}")
print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
print(f"PyTorch version: {torch.__version__}")
print(f"CUDA available: {torch.cuda.is_available()}")
print(f"CUDA device count: {torch.cuda.device_count()}")
if torch.cuda.is_available():
    print(f"CUDA current device: {torch.cuda.current_device()}")
    print(f"CUDA device name: {torch.cuda.get_device_name(0)}")
    
    # Check for GPU memory issues
    for i in range(torch.cuda.device_count()):
        props = torch.cuda.get_device_properties(i)
        print(f"  GPU {i}: {props.name}, Total Memory: {props.total_memory / 1e9:.1f} GB")

print("\nAttempting model load...")
sys.stdout.flush()

# Load model with detailed error handling
try:
    print(f"[{time.strftime('%H:%M:%S')}] Loading ESM model...")
    sys.stdout.flush()
    
    model = esm.pretrained.esmfold_v1()
    print(f"[{time.strftime('%H:%M:%S')}] Model loaded, dtype: {next(model.parameters()).dtype}")
    print(f"[{time.strftime('%H:%M:%S')}] Model device before to(): {next(model.parameters()).device}")
    
    # Check model default precision
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    print(f"[{time.strftime('%H:%M:%S')}] Target device: {device}")
    
    # Move to device - THIS IS WHERE THE ERROR OCCURS
    print(f"[{time.strftime('%H:%M:%S')}] Moving model to device...")
    sys.stdout.flush()
    
    model = model.eval().to(device)
    print(f"[{time.strftime('%H:%M:%S')}] Model on device, dtype: {next(model.parameters()).dtype}")
    
    # Check if model is in half precision
    is_half = next(model.parameters()).dtype == torch.float16
    print(f"[{time.strftime('%H:%M:%S')}] Model is using Half (float16): {is_half}")
    
    if is_half:
        print("WARNING: Model is in Half precision (float16)")
        print("LayerNorm operations may fail on certain GPU configurations")
        
except Exception as e:
    print(f"\n{'='*80}")
    print(f"ERROR during model loading: {e}")
    print(f"{'='*80}")
    print("\nDEBUG INFO:")
    print(f"Current device: {torch.cuda.current_device() if torch.cuda.is_available() else 'CPU'}")
    print(f"CUDA_VISIBLE_DEVICES: {os.environ.get('CUDA_VISIBLE_DEVICES', 'NOT SET')}")
    sys.stdout.flush()
    raise

print("\n" + "="*80)
print("Model loaded successfully!")
print("="*80 + "\n")

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
Path(args.output_folder).mkdir(parents=True, exist_ok=True)

# Load MPNN_csv as a df
MPNN_df = pd.read_csv(args.input_csv)

# Parse data frame
for i, row in tqdm(MPNN_df.iterrows(), total=len(MPNN_df), desc="Predicting Structures"):
    sequence = row['seq']
    file_name = row['file_ID']
    
    with torch.no_grad():
        output = model.infer_pdb(sequence)

    with open(f'{args.output_folder}/{file_name}.pdb', "w") as f:
        f.write(output)

print(f"\n[{time.strftime('%H:%M:%S')}] Job completed successfully")
