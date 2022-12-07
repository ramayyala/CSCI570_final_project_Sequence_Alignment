#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20GB
#SBATCH --time=00:01:00
python3 ../code/efficient.py ../data/datapoints/in9.txt ../outputs/in9.txt_output.txt
