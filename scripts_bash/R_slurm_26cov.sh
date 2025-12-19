#!/bin/bash
#SBATCH --job-name R_job_26cov
#SBATCH --output R_job_output_26cov.log
#SBATCH --error R_job_error_26cov.log
#SBATCH --cpus-per-task 1
#SBATCH --mem 1024G
 
module load R/4.4.1
 
Rscript /lustre/home/tpyork/projects/preg-maternal/scripts/05_mwas-timepoint-day_26cov_578n_p.R