#!/bin/bash
#SBATCH --job-name R_job_resid
#SBATCH --output R_job_output_resid.log
#SBATCH --error R_job_error_resid.log
#SBATCH --cpus-per-task 1
#SBATCH --mem 1024G
 
module load R/4.4.1
 
Rscript /lustre/home/tpyork/projects/preg-maternal/scripts/03_pca-data.R