"""

"""
### MODULES ############################################################################
import pyrosetta
from pyrosetta.rosetta.core.pack.guidance_scoreterms.sap import calculate_sap
from pyrosetta.rosetta.core.select.residue_selector import TrueResidueSelector

### INPUT CHECK #######################################################################

### MAIN ##############################################################################
# 1. Initialize Rosetta (standard init)
pyrosetta.init("-mute all")

def get_sap_score(pdb_path, radius=5.0):
    """
    Calculates the Spatial Aggregation Propensity (SAP) score.
    
    Args:
        pdb_path (str): Path to PDB file
        radius (float): The radius (Angstroms) to consider for hydrophobic neighbors.
                        5.0A is standard for antibody developability.
    
    Returns:
        float: The global SAP score.
    """
    try:
        # Load the structure
        pose = pyrosetta.pose_from_pdb(pdb_path)
        
        # 2. Setup the Selector
        # SAP requires a residue selector to know which residues to score.
        # TrueResidueSelector selects the whole protein.
        selector = TrueResidueSelector()
        
        # 3. Calculate Score
        # calculate_sap(pose, residue_selector, radius)
        # Note: This function allows you to mask specific regions if needed.
        sap_score = calculate_sap(pose, selector, radius)
        
        return sap_score
        
    except Exception as e:
        print(f"Error on {pdb_path}: {e}")
        return None

# Example Usage
if __name__ == "__main__":
    score = get_sap_score("my_design.pdb")
    print(f"SAP Score: {score}")