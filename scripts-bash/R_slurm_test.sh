#!/bin/bash
#SBATCH --job-name R_job
#SBATCH --output R_job_output.log
#SBATCH --error R_job_error.log
#SBATCH --cpus-per-task 1
#SBATCH --mem 4G
 
module load R/4.4.1
 
Rscript /lustre/home/tpyork/projects/preg-maternal/scripts/00_test-bash.R