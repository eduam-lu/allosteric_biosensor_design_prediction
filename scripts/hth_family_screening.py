"""
Script to screen the LacI family (PF00356) for candidates with significant conformational changes.
We mine the RCSB PDB for all structures with this Pfam family, group them by UniProt,
and analyze pairs of structures for hinge shifts at the DNA-binding domain.
"""
### MODULES ###################################################################################################################
import requests
import numpy as np
import warnings
import os
import itertools
from Bio.PDB import PDBList, PDBParser, Superimposer
import json
import sys
# Suppress warnings
warnings.filterwarnings('ignore')

### PARAMETERS ################################################################################################################

### FUNCTIONS #################################################################################################################

def mine_laci_family_data():
    """
    MINING VERSION:
    1. Searches RCSB for all structures with Pfam PF00356 (LacI Family).
    2. Maps PDB IDs to UniProt Accessions to group them by biological protein.
    3. Filters out singletons (proteins with only 1 structure), as they can't form pairs.
    
    Returns:
        dict: { 'UniProtID (Organism)': [List of PDB IDs] }
    """
    print("Step 1: Searching RCSB for all LacI family members (PF00356)...")
    
    # --- 1. SEARCH QUERY ---
    # We ask for "polymer_entity" to get chain-specific hits (e.g., 1EFA_1)
    search_query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                # Requirement 1: Must be LacI Head
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_annotation.annotation_id",
                        "operator": "exact_match",
                        "value": "PF00356"
                    }
                },
                # Requirement 2: Must be one of the Valid Bodies
                {
                    "type": "group",
                    "logical_operator": "or",
                    "nodes": [
                        {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_polymer_entity_annotation.annotation_id", "operator": "exact_match", "value": "PF13377"}}, # Peripla_BP_3 (LacI body)
                        {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_polymer_entity_annotation.annotation_id", "operator": "exact_match", "value": "PF00532"}}, # Peripla_BP_1
                        {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_polymer_entity_annotation.annotation_id", "operator": "exact_match", "value": "PF13407"}}  # Peripla_BP_4
                    ]
                },
                # Requirement 3: Good Resolution
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {"attribute": "rcsb_entry_info.resolution_combined", "operator": "less", "value": 3.5}
                }
            ]
        },
        "return_type": "polymer_entity",
        "request_options": {"return_all_hits": True}
    }
    
    # Execute Search
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    try:
        resp = requests.post(search_url, json=search_query)
        resp.raise_for_status()
        # Extract Entity IDs (e.g., "1EFA_1")
        entity_ids = [x['identifier'] for x in resp.json().get('result_set', [])]
    except Exception as e:
        print(f"Search failed: {e}")
        return {}

    print(f"Found {len(entity_ids)} matching entities. Grouping by UniProt...")

    # --- 2. BATCH FETCH METADATA (GraphQL) ---
    # We need to map 1EFA_1 -> UniProt ID. Doing this one-by-one is slow.
    # We use RCSB GraphQL to do it in batches of 50.
    
    grouped_data = {} # { 'UniProt_ID': {'organism': 'Name', 'pdbs': set()} }
    
    def chunks(lst, n):
        """Yield successive n-sized chunks from lst."""
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    graphql_url = "https://data.rcsb.org/graphql"
    
    # GraphQL Query to fetch UniProt ID and Organism Name
    gql_query = """
    query($ids: [String!]!) {
      polymer_entities(entity_ids: $ids) {
        rcsb_id
        rcsb_polymer_entity_container_identifiers {
          entry_id
          reference_sequence_identifiers {
            database_accession
            database_name
          }
        }
        rcsb_entity_source_organism {
          scientific_name
        }
      }
    }
    """

    for chunk in chunks(entity_ids, 50):
        try:
            r = requests.post(graphql_url, json={'query': gql_query, 'variables': {'ids': chunk}})
            data = r.json().get('data', {}).get('polymer_entities', [])
            
            for entity in data:
                if not entity: continue
                
                pdb_id = entity['rcsb_polymer_entity_container_identifiers']['entry_id']
                
                # Get UniProt ID
                refs = entity['rcsb_polymer_entity_container_identifiers'].get('reference_sequence_identifiers', [])
                uniprot_id = next((x['database_accession'] for x in refs if x['database_name'] == 'UniProt'), None)
                
                # Get Organism Name
                try:
                    org_name = entity['rcsb_entity_source_organism'][0]['scientific_name']
                except:
                    org_name = "Unknown"

                if uniprot_id:
                    if uniprot_id not in grouped_data:
                        grouped_data[uniprot_id] = {'organism': org_name, 'pdbs': set()}
                    grouped_data[uniprot_id]['pdbs'].add(pdb_id)
                    
        except Exception as e:
            continue

    # --- 3. FORMAT & FILTER ---
    final_candidates = {}
    
    for uid, data in grouped_data.items():
        pdb_list = list(data['pdbs'])
        
        # FILTER: Only keep candidates with at least 2 structures (Apo/Holo potential)
        if len(pdb_list) >= 2:
            key_name = f"{uid} ({data['organism']})"
            final_candidates[key_name] = pdb_list

    print(f"Screening complete. Identified {len(final_candidates)} unique candidates with >=2 structures.")
    return final_candidates

def get_laci_family_data():
    """
    Fetches lists of PDBs for key LacI family members.
    In a full production run, you would mine UniProt for all PF00356 members.
    Here we seed the most common candidates to ensure the script is immediately runnable.
    """
    # Dictionary of { 'CommonName': [List of PDBs] }
    candidates = {
        "LacI (E.coli)": ["1EFA", "1LBI", "1LBG", "2PFA", "1JWL"],
        "PurR (E.coli)": ["1WET", "1J5N", "1QPZ", "2PUA", "1PNR"],
        "CcpA (B.subtilis)": ["1RQP", "1QVZ", "1Z04", "3OQM"],
        "FruR (E.coli)": ["1LKP", "1LKD", "1L8E"],
        "TreR (E.coli)": ["1BYI", "2JF9"], 
        "GalR (E.coli)": ["1B7A"] # Limited structures available
    }
    return candidates

def download_and_parse(pdb_id, pdir='pdbs'):
    """Robust PDB downloader/parser."""
    pdbl = PDBList(verbose=False)
    parser = PDBParser(QUIET=True)
    try:
        if not os.path.exists(pdir): os.makedirs(pdir)
        fname = pdbl.retrieve_pdb_file(pdb_id, pdir=pdir, file_format='pdb')
        return parser.get_structure(pdb_id, fname)
    except Exception as e:
        return None

def get_hinge_residue_id(name):
    """Returns the approximate Hinge residue index for known families."""
    if "PurR" in name: return 54
    if "CcpA" in name: return 58
    return 62  # Default for LacI and generic homologs

def analyze_candidate(name, pdb_list):
    """
    Analyzes a list of PDBs for a single protein.
    Returns the pair with the MAXIMAL hinge shift (The Switch).
    """
    if len(pdb_list) < 2: return None

    print(f"Analyzing {name} ({len(pdb_list)} structures)...")
    structures = []
    
    # 1. Load structures
    for pid in pdb_list:
        s = download_and_parse(pid)
        if s: structures.append((pid, s))
    
    if len(structures) < 2: return None

    hinge_res_id = get_hinge_residue_id(name)
    widths = {}
    valid_structs = []
    
    # 2. Measure Hinge Width for every structure
    # Defined as distance between CA of Chain A and Chain B at the hinge point
    for pid, s in structures:
        try:
            chains = [c for c in s[0]] # Get chains in Model 0
            if len(chains) < 2: continue # Skip monomers
            
            c1, c2 = chains[0], chains[1] 
            
            if hinge_res_id in c1 and hinge_res_id in c2:
                v1 = c1[hinge_res_id]['CA'].get_vector()
                v2 = c2[hinge_res_id]['CA'].get_vector()
                dist = (v1 - v2).norm()
                widths[pid] = dist
                valid_structs.append((pid, s))
        except:
            continue
            
    # 3. Find the pair with Max Delta
    max_delta = 0
    best_pair = None
    best_bfactor = 100
    
    for (pid1, s1), (pid2, s2) in itertools.combinations(valid_structs, 2):
        w1 = widths[pid1]
        w2 = widths[pid2]
        delta = abs(w1 - w2)
        
        if delta > max_delta:
            max_delta = delta
            best_pair = (pid1, pid2)
            
            # 4. Calculate Rigidity (B-Factor) of the N-subdomain
            # We check the "Wider" (Induced) structure, as this is often the template for RFdiffusion
            target_struc = s1 if w1 > w2 else s2
            atoms = []
            # N-subdomain range approx 60-150
            for r in target_struc[0]['A']:
                if 60 <= r.get_id()[1] <= 150: 
                     for a in r: atoms.append(a.get_bfactor())
            
            best_bfactor = np.mean(atoms) if atoms else 99

    if not best_pair: return None

    return {
        "Name": name,
        "Max_Shift_Angstrom": max_delta,
        "Pair": best_pair,
        "Anchor_B_Factor": best_bfactor
    }

### EXECUTION ##############################################################################################################
candidates = mine_laci_family_data()
json.dump(candidates, open('candidates.json', 'w'))
sys.exit()
results = []

for name, pdbs in candidates.items():
    res = analyze_candidate(name, pdbs)
    if res: results.append(res)

# Sort by Shift Magnitude (The most important metric for redesign)
results.sort(key=lambda x: x['Max_Shift_Angstrom'], reverse=True)

# Output Table
print("\n" + "="*75)
print(f"{'CANDIDATE':<20} | {'SHIFT (A)':<10} | {'B-FACTOR':<10} | {'BEST PAIR'}")
print("="*75)

for r in results:
    print(f"{r['Name']:<20} | {r['Max_Shift_Angstrom']:.2f}       | {r['Anchor_B_Factor']:.2f}       | {r['Pair'][0]} / {r['Pair'][1]}")
print("="*75)