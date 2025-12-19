#!/bin/bash
#SBATCH --job-name R_job_24cov
#SBATCH --output R_job_output_24cov.log
#SBATCH --error R_job_error_24cov.log
#SBATCH --cpus-per-task 1
#SBATCH --mem 1024G
 
module load R/4.4.1
 
Rscript /lustre/home/tpyork/projects/preg-maternal/scripts/05_mwas-timepoint-day_24cov_625n.R