#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p node
#SBATCH -t 1-00:00:00
#SBATCH -J score_mock_devel

CMD="make score_mul_0_0_mock_"$1".tsv"
$CMD
