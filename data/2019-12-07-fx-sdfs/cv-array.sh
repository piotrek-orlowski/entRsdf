#!/bin/bash
#SBATCH --account=def-pior
#SBATCH --time=24:00:00
#SBATCH --array=1-4
#SBATCH --checkpoint=2:0:0
#SBATCH --ntasks=12
#SBATCH --mem-per-cpu=512M

export NODELIST=$(echo $(srun hostname | cut -f 1 -d '.'))

echo "Starting task $SLURM_ARRAY_TASK_ID"
ARGS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" array-cv-args)

Rscript --no-save run-full-sample-sdf.R $ARGS
