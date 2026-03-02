#!/bin/bash

# Usage: ./export_conda.sh /path/to/destination
TARGET_DIR=${1:-"./conda_exports"}

# Create the directory if it doesn't exist
mkdir -p "$TARGET_DIR"

echo "Exporting conda environments to: $TARGET_DIR"

# Get a list of all conda environments (names only)
# We filter out lines starting with # and empty lines
envs=$(conda env list | grep -v '^#' | awk '{print $1}')

for env in $envs; do
    # Skip 'base' or any empty strings if necessary
    if [ -z "$env" ]; then
        continue
    fi

    echo "Exporting $env..."
    
    # Export the environment to a yaml file
    # --no-builds is optional but makes files more portable across OS versions
    conda env export -n "$env" --no-builds > "$TARGET_DIR/${env}.yml"
    
    if [ $? -eq 0 ]; then
        echo "Successfully exported $env"
    else
        echo "Failed to export $env"
    fi
done

echo "Done! All environments are located in $TARGET_DIR"