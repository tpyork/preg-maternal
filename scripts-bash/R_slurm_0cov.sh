#!/bin/bash
#SBATCH --job-name R_job_0cov
#SBATCH --output R_job_output_0cov.log
#SBATCH --error R_job_error_0cov.log
#SBATCH --cpus-per-task 1
#SBATCH --mem 256G
 
module load R/4.4.1
 
Rscript /lustre/home/tpyork/projects/preg-maternal/scripts/05_mwas-timepoint-day_0cov.R