#!/bin/bash
#SBATCH --job-name R_job_cv_9
#SBATCH --output R_job_output_cv_9.log
#SBATCH --error R_job_error_cv_9.log
#SBATCH --cpus-per-task 1
#SBATCH --mem 1024G
 
module load R/4.4.1
 
Rscript /lustre/home/tpyork/projects/preg-maternal/scripts_pred/05_mwas_cv_9.R