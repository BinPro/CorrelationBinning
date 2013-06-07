#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p node
#SBATCH -t 1-00:00:00
#SBATCH -J score_cv_mul

CMD="make score_mul_0_0_cv_"$1".tsv"
$CMD
