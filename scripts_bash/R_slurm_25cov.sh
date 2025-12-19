#!/bin/bash
#SBATCH --job-name R_job_25cov
#SBATCH --output R_job_output_25cov.log
#SBATCH --error R_job_error_25cov.log
#SBATCH --cpus-per-task 1
#SBATCH --mem 1024G
 
module load R/4.4.1
 
Rscript /lustre/home/tpyork/projects/preg-maternal/scripts/05_mwas-timepoint-day_25cov_578n_p.R