#!/usr/bin/env bash
#SBATCH -A b2010008
#SBATCH -p node
#SBATCH -C mem72GB
#SBATCH -n 8
#SBATCH -t 1-00:00:00
#SBATCH -J ProBin
#SBATCH --mail-type=ALL
#SBATCH --mail-user=brynjar.bjarnason@scilifelab.se

make mock-even-old 
