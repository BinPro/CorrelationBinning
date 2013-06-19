#!/bin/bash
#SBATCH -A b2010008
#SBATCH -p node
#SBATCH -t 1:00:00
#SBATCH -J score_combined_factor_100_core_real

make score_combined_real_46_columns_6_100.tsv