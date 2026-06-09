#!/bin/bash
# Date
# 25/03/2026

# Description
# Launch 12 jobs with Samuel pdbs in 3 batches of 4 parallel jobs
# Each batch uses all 16 GPUs (4 GPUs per job, 4 jobs in parallel)


#SBATCH -t 0-18:00:00
#SBATCH -J trial_sbatch_run
#SBATCH -A Berzelius-2025-416  # REQUIRED on Berzelius
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:A100-SXM4-40GB:8
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --exclude=node014,node020
#SBATCH --array=1-3

#SBATCH -o trial_sbatch_run_%j.out
#SBATCH -e trial_sbatch_run_%j.err
#SBATCH --mail-user=edu.amgo@gmail.com
#SBATCH --mail-type=ALL

# Print script to log for debugging
cat $0

set -eo pipefail

# Load environment
module load Miniforge3/24.7.1-2-hpc1-bdist
conda activate allostery 
set -u

# Define paths as Shell Variables (No spaces around '=')
pipeline_path="/proj/berzelius-2023-361/users/x_eduam/allosteric_biosensor_design_prediction/strategy_1_snakemake/bsd_def.smk"
config_base="/proj/berzelius-2023-361/users/x_eduam/strategy_1_campaign/configs/config"
report_base="/proj/berzelius-2023-361/users/x_eduam/strategy_1_campaign/reports/report"

# Calculate which jobs to run in this batch
# Batch 1 (array_id=1): jobs 1,2,3,4
# Batch 2 (array_id=2): jobs 5,6,7,8  
# Batch 3 (array_id=3): jobs 9,10,11,12
job_start=$(( (SLURM_ARRAY_TASK_ID - 1) * 4 + 1 ))
job_end=$(( SLURM_ARRAY_TASK_ID * 4 ))

echo "========================================="
echo "Running batch ${SLURM_ARRAY_TASK_ID}: jobs ${job_start}-${job_end}"
echo "========================================="

# NOTE: Do NOT copy cache here on master node - each compute node has its own /tmp!
# We'll copy on each compute node via srun instead
echo "Will cache ESM models to local /tmp on each compute node during job execution"

# CRITICAL: Clean up old GPU leases that might block even-numbered jobs
# GPU lock files from previous runs can cause resource contention
for job_id in $(seq $job_start $job_end); do
    config_path="${config_base}_${job_id}.yml"
    # Extract output_dir from config file dynamically
    output_dir=$(grep "^output_dir:" "$config_path" | sed 's/output_dir: "\(.*\)"/\1/')
    gpu_lease_dir="${output_dir}/.gpu_leases"
    if [ -d "$gpu_lease_dir" ]; then
        rm -rf "$gpu_lease_dir"
        echo "Cleaned old GPU leases for job ${job_id}: $gpu_lease_dir"
    fi
done

# Run 4 jobs in parallel using background processes with GPU pinning
# CRITICAL: Each job gets its own work directory AND explicit GPU assignment
# GPU pinning prevents resource contention between parallel jobs
declare -a pids
declare -a job_ids
for job_id in $(seq $job_start $job_end); do
    config_path="${config_base}_${job_id}.yml"
    work_dir="/tmp/snakemake_job_${job_id}_$$"
    
    # Calculate which GPUs this job should use (0-3, 4-7, 8-11, or 12-15)
    # Based on job position within the batch
    job_offset=$(( job_id - job_start ))
    gpu_start=$(( job_offset * 4 ))
    gpu_end=$(( gpu_start + 3 ))
    cuda_visible_devices="${gpu_start},$(( gpu_start + 1 )),$(( gpu_start + 2 )),$(( gpu_start + 3 ))"
    
    echo "Starting job ${job_id} with config: $config_path"
    echo "  Work directory: $work_dir"
    echo "  GPU Assignment: devices ${cuda_visible_devices}"
    
    mkdir -p "$work_dir"
    
    # Launch each pipeline with GPU pinning using CUDA_VISIBLE_DEVICES
    # CRITICAL: Pass CUDA_VISIBLE_DEVICES directly to srun (not via export)
    srun --exclusive --ntasks=1 --nodes=1 --cpus-per-task="${SLURM_CPUS_PER_TASK}" --gres=gpu:A100-SXM4-40GB:4 \
        --export=ALL,CUDA_VISIBLE_DEVICES="${cuda_visible_devices}" \
        bash -c "set -euo pipefail; \
        LOCAL_NODE_CACHE=/tmp/torch_cache_${job_id}; \
        mkdir -p \$LOCAL_NODE_CACHE; \
        cp -r '$TORCH_HOME/hub' \$LOCAL_NODE_CACHE/ 2>/dev/null || true; \
        export OMP_NUM_THREADS=1 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=2; \
        export TORCH_HOME=\$LOCAL_NODE_CACHE; \
        echo 'CUDA_VISIBLE_DEVICES in subprocess: '\$CUDA_VISIBLE_DEVICES && \
        snakemake -s '$pipeline_path' --configfile '$config_path' --cores '${SLURM_CPUS_PER_TASK}' --resources gpu=4 --directory '$work_dir'; \
        rm -rf \$LOCAL_NODE_CACHE; \
        touch '$work_dir/.pipeline_ok'" &

    pids+=("$!")
    job_ids+=("$job_id")
done

# Wait for all background jobs to complete
echo "Waiting for all ${#pids[@]} background jobs to complete..."
pipeline_failed=0
for idx in "${!pids[@]}"; do
    pid="${pids[$idx]}"
    job_id="${job_ids[$idx]}"
    if wait "$pid"; then
        echo "✓ Pipeline job ${job_id} completed successfully"
    else
        exit_code=$?
        echo "✗ Pipeline job ${job_id} failed with exit code ${exit_code}"
        pipeline_failed=1
    fi
done

echo "Batch ${SLURM_ARRAY_TASK_ID} pipelines completed"

# CRITICAL: Force file system sync and clear caches
# This ensures all GPU locks and job outputs are fully written to disk
echo "Syncing file system and clearing caches..."
sync
sync
sleep 2

# Verify that output directories exist before report generation
echo "Verifying output directories for batch ${SLURM_ARRAY_TASK_ID}..."
for job_id in $(seq $job_start $job_end); do
    config_path="${config_base}_${job_id}.yml"
    # Extract output_dir from config file dynamically
    output_dir=$(grep "^output_dir:" "$config_path" | sed 's/output_dir: "\(.*\)"/\1/')
    if [ -d "$output_dir" ]; then
        echo "  ✓ Job ${job_id}: $output_dir exists"
    else
        echo "  ✗ Job ${job_id}: $output_dir MISSING"
    fi
done

# Generate reports for all jobs in this batch
# CRITICAL: Use the same work directory where the pipeline ran to access checkpoints!
echo "Generating reports for batch ${SLURM_ARRAY_TASK_ID}..."
report_failed=0
for job_id in $(seq $job_start $job_end); do
    config_path="${config_base}_${job_id}.yml"
    report_output="${report_base}_${job_id}_2403.html"
    work_dir="/tmp/snakemake_job_${job_id}_$$"

    echo "  Job ${job_id}: Generating report..."
    
    # Always attempt to generate reports - Snakemake will handle jobs that didn't complete
    # In parallel execution, file visibility can be delayed, so don't skip based on marker files
    if snakemake -s "$pipeline_path" --configfile "$config_path" --report "$report_output" \
        --directory "$work_dir" 2>&1 | tail -5; then
        
        if [ -f "$report_output" ]; then
            report_size=$(du -h "$report_output" | cut -f1)
            echo "  ✓ Job ${job_id}: Report generated (${report_size})"
        else
            echo "  ✗ Job ${job_id}: Report command passed but file not created"
            report_failed=1
        fi
    else
        echo "  ✗ Job ${job_id}: Report generation failed"
        report_failed=1
    fi
done

echo "All reports for batch ${SLURM_ARRAY_TASK_ID} generation completed"

if [ "$report_failed" -ne 0 ]; then
    echo "WARNING: Some reports were not generated successfully in batch ${SLURM_ARRAY_TASK_ID}"
fi

if [ "$pipeline_failed" -ne 0 ]; then
    echo "One or more pipeline jobs failed in batch ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi