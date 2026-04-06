#!/usr/bin/env python3
"""
Validation script for RF Diffusion Binding Site Designer Snakemake Pipeline
Checks configuration, dependencies, and paths before running the pipeline
"""

import os
import sys
import json
import argparse
from pathlib import Path

try:
    import yaml
    YAML_OK = True
except ImportError:
    YAML_OK = False

def print_section(title):
    """Print a formatted section header."""
    print(f"\n{'='*60}")
    print(f"  {title}")
    print(f"{'='*60}")

def check_mark(condition, message):
    """Print check result with mark."""
    symbol = "✓" if condition else "✗"
    color = "\033[92m" if condition else "\033[91m"
    reset = "\033[0m"
    print(f"{color}{symbol}{reset} {message}")
    return condition

def check_dependencies():
    """Check if required Python packages are installed."""
    print_section("CHECKING PYTHON DEPENDENCIES")
    
    dependencies = {
        'Bio': 'Biopython',
        'numpy': 'NumPy',
        'pandas': 'Pandas',
        'rdkit': 'RDKit',
        'tmtools': 'TMTools',
        'yaml': 'PyYAML'
    }
    
    all_ok = True
    for module, name in dependencies.items():
        try:
            __import__(module)
            check_mark(True, f"{name} installed")
        except ImportError:
            check_mark(False, f"{name} NOT FOUND - Install with: pip install {name.lower()}")
            all_ok = False
    
    # Check optional dependencies
    print("\nOptional dependencies:")
    optional = {
        'pymol': 'PyMOL (for ligand detection)',
        'snakemake': 'Snakemake (required for running)',
        'biotite': 'Biotite (for structure I/O)'
    }
    
    for module, name in optional.items():
        try:
            __import__(module)
            check_mark(True, f"{name} available")
        except ImportError:
            check_mark(False, f"{name} not available")
    
    return all_ok

def check_config(config_file):
    """Validate configuration file."""
    print_section("CHECKING CONFIGURATION FILE")
    
    # Check file exists
    if not check_mark(os.path.exists(config_file), f"Config file exists: {config_file}"):
        return False
    
    # Try to load YAML
    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        check_mark(True, "Config YAML is valid")
    except Exception as e:
        check_mark(False, f"Config YAML error: {e}")
        return False
    
    # Check required parameters
    required_params = ['structure_path', 'ligand_smiles', 'ligand_name']
    missing = []
    
    for param in required_params:
        if param not in config:
            missing.append(param)
            check_mark(False, f"Missing required parameter: {param}")
        else:
            value = config[param]
            if value is None or str(value).startswith('path/to/'):
                check_mark(False, f"Parameter not set or is placeholder: {param} = {value}")
            else:
                check_mark(True, f"Parameter set: {param}")
    
    return len(missing) == 0

def check_paths(config_file):
    """Check if all specified paths exist."""
    print_section("CHECKING PATHS IN CONFIGURATION")
    
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    paths_to_check = [
        'structure_path',
        'path_to_RFAA_apptainer',
        'path_to_RFAA_script',
        'path_to_RFAA_weights',
        'path_to_RF1_script',
        'path_to_RF1_env',
        'path_to_ligand_MPNN_script',
        'path_to_ligand_MPNN_env',
        'path_to_ESM_env',
        'path_to_ESM_script',
        'path_to_boltz_env'
    ]
    
    warnings = []
    
    for path_key in paths_to_check:
        if path_key in config and config[path_key]:
            path = config[path_key]
            if str(path).startswith('path/to/'):
                check_mark(False, f"Placeholder path: {path_key}")
                warnings.append(f"  → Update {path_key} to actual path")
            elif os.path.exists(path):
                check_mark(True, f"Path exists: {path_key} = {path}")
            else:
                check_mark(False, f"Path NOT FOUND: {path_key} = {path}")
                warnings.append(f"  → {path_key} points to non-existent path")
    
    if warnings:
        print("\n⚠️  Warnings to address:")
        for w in warnings:
            print(w)
    
    return len(warnings) == 0

def check_files():
    """Check if required files exist in pipeline directory."""
    print_section("CHECKING PIPELINE FILES")
    
    required_files = [
        'rf_diffusion_bsd_snakemake.smk',
        'config_rfdiffusion_bsd.yml',
        'functions_bsd.py',
        'ligand_sampling.py'
    ]
    
    all_ok = True
    for filename in required_files:
        exists = os.path.exists(filename)
        check_mark(exists, f"File present: {filename}")
        all_ok = all_ok and exists
    
    return all_ok

def check_external_tools():
    """Check if external tools are available."""
    print_section("CHECKING EXTERNAL TOOLS")
    
    tools = [
        ('dssp', 'mkdssp --version', 'DSSP (secondary structure)'),
        ('gnina', 'gnina --help', 'GNINA (molecular docking)'),
    ]
    
    all_ok = True
    for tool_name, command, description in tools:
        result = os.system(f"{command} > /dev/null 2>&1")
        is_available = result == 0
        check_mark(is_available, f"{description}: {tool_name}")
        if not is_available:
            print(f"  → Install with: apt install {tool_name} (or conda install {tool_name})")
        all_ok = all_ok and is_available
    
    return all_ok

def generate_summary(results):
    """Generate summary of checks."""
    print_section("VALIDATION SUMMARY")
    
    checks = [
        ("Dependencies", results.get('dependencies', False)),
        ("Configuration", results.get('config', False)),
        ("Paths", results.get('paths', False)),
        ("Files", results.get('files', False)),
    ]
    
    passed = sum(1 for _, result in checks if result)
    total = len(checks)
    
    for name, result in checks:
        symbol = "✓" if result else "✗"
        color = "\033[92m" if result else "\033[91m"
        reset = "\033[0m"
        print(f"{color}{symbol}{reset} {name}")
    
    print(f"\n{passed}/{total} checks passed")
    
    if passed == total:
        print("\n" + "="*60)
        print("🎉 All checks passed! Ready to run the pipeline!")
        print("="*60)
        print("\nRun the pipeline with:")
        print("  snakemake --snakefile rf_diffusion_bsd_snakemake.smk \\")
        print("    --configfile config_rfdiffusion_bsd.yml \\")
        print("    -c 1 -p")
        return True
    else:
        print("\n" + "="*60)
        print("⚠️  Some checks failed. Please fix the issues above.")
        print("="*60)
        return False

def main():
    """Main validation routine."""
    parser = argparse.ArgumentParser(
        description="Validate RF Diffusion Binding Site Designer Snakemake Pipeline"
    )
    parser.add_argument(
        '--config',
        default='config_rfdiffusion_bsd.yml',
        help='Path to configuration file (default: config_rfdiffusion_bsd.yml)'
    )
    parser.add_argument(
        '--skip-external',
        action='store_true',
        help='Skip external tool checks'
    )
    
    args = parser.parse_args()
    
    print("\n" + "="*60)
    print("  RF DIFFUSION BINDING SITE DESIGNER")
    print("  Pipeline Validation Script")
    print("="*60)
    
    results = {}
    
    # Run checks
    results['dependencies'] = check_dependencies()
    
    if os.path.exists(args.config):
        results['config'] = check_config(args.config)
        results['paths'] = check_paths(args.config) if results['config'] else False
    else:
        check_mark(False, f"Config file not found: {args.config}")
        results['config'] = False
        results['paths'] = False
    
    results['files'] = check_files()
    
    if not args.skip_external:
        check_external_tools()
    
    # Generate summary and exit
    success = generate_summary(results)
    sys.exit(0 if success else 1)

if __name__ == '__main__':
    main()
