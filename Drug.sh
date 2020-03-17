#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=16000
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=72:00:00
#SBATCH --ntasks-per-node=1
#SBATCH -A dutta_tumor
#SBATCH --partition=standard

module load gcc R
echo ${SLURM_ARRAY_TASK_ID}
Rscript RivannaDrugAnalysis.R ${SLURM_ARRAY_TASK_ID}
