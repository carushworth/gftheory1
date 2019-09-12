#!/bin/bash -l
#SBATCH -D /home/crusher/projects/GFs/scripts
#SBATCH -J gf_initialrun
#SBATCH -o /home/crusher/projects/GFs/slurm-log/GFs%A_%a.out
#SBATCH -e /home/crusher/projects/GFs/slurm-log/GFs%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --array=1-15
#SBATCH --time=7-00:00
#SBATCH --mail-user=crushworth@ucdavis.edu
set -e
set -u

INPUT=$(awk "NR==$SLURM_ARRAY_TASK_ID" scripts/param_combos_initial.txt)

Rscript paramspace_test2.R $INPUT