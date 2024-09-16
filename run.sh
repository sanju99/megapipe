#!/bin/bash 
#SBATCH -c 8
#SBATCH -t 1-00:00
#SBATCH -p medium
#SBATCH --mem=100G
#SBATCH -o /home/sak0914/Errors/zerrors_%j.out 
#SBATCH -e /home/sak0914/Errors/zerrors_%j.err 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skulkarni@g.harvard.edu

# snakemake --configfile config.yaml --cores 8 --rerun-incomplete --keep-going
snakemake --cores 8 --use-conda --conda-frontend conda --rerun-incomplete --keep-going