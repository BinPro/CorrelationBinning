#!/bin/bash
#SBATCH -A b2010008
#SBATCH -p node
#SBATCH -t 12:00:00
#SBATCH -J score_time_series_real

make score_time_series_real_only_iso.tsv