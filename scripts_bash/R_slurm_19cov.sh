#!/bin/bash
#SBATCH --job-name R_job_19cov
#SBATCH --output R_job_output_19cov.log
#SBATCH --error R_job_error_19cov.log
#SBATCH --cpus-per-task 1
#SBATCH --mem 512G
 
module load R/4.4.1
 
Rscript /lustre/home/tpyork/projects/preg-maternal/scripts/05_mwas-timepoint-day_19cov.R