"""
binding_site_designer.py --structure --ligand --output --help

DESCRIPTION

This script, given an input structure with a ligand bound to it, performs a binding site redesign around a given new ligand, in order to accept it.
Binding site redesign is performed through RF diffusion all atom and LigandMPNN. Candidates are filtered based on prediction, 3D and binding metrics

FLAGS

--structure:

--ligand:

--output:

--help:

WORKFLOW


Eduardo Amo González
2025-2026
"""
### IMPORT MODULES ###########################################################################################################################

import functions_bsd as func
import argparse
import pandas
from pathlib import Path
import sys
import logging

### PARAMS ####################################################################################################################################

### INPUT CHECK ###############################################################################################################################

parser = argparse.ArgumentParser(
    description="This script runs the pipeline shapedesign with an additional RF_diffusion step"
)
parser.add_argument('--folder', help="Folder that contains all the input structures", type=str)
parser.add_argument('--detailed-help', action='store_true', help="Show detailed help message and exit")

args = parser.parse_args()


### SET UP LOG FILE #############################################################################################################################
logging.basicConfig(
    filename=f"{str(global_output_name)}/improved_shapedesign.log",       # Log file path
    filemode='a',                   # Append mode
    format='%(asctime)s - %(message)s',
    level=logging.INFO
)

def save_to_log(message):
    logging.info(message)

### MAIN EXECUTION ##############################################################################################################################
