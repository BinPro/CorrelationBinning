#!/bin/bash
#SBATCH -A b2010008
#SBATCH -p devel
#SBATCH -t 1:00:00
#SBATCH -J score_combined_factor_10_devel

make score_combined_cv_46_columns_6_10.tsv