#!/bin/bash
#SBATCH --job-name R_job_assess
#SBATCH --output R_job_output_assess.log
#SBATCH --error R_job_error_assess.log
#SBATCH --cpus-per-task 1
#SBATCH --mem 1024G
 
module load R/4.4.1
 
Rscript /lustre/home/tpyork/projects/preg-maternal/scripts_cv/10_cv-assessment-fxn-qval.R