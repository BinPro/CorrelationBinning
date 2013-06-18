#!/bin/bash
#SBATCH -A b2010008
#SBATCH -p node
#SBATCH -t 12:00:00
#SBATCH -J score_time_series_test

make score_time_series_even_only_iso_2_columns.tsv