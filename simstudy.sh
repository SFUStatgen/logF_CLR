#!/bin/bash
#SBATCH --account=def-mcneney
#SBATCH --array=1-480
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=20:00:00
module load r
echo "This is job $SLURM_ARRAY_TASK_ID out of $SLURM_ARRAY_TASK_COUNT jobs."
Rscript simstudy.R
