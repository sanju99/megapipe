#!/bin/bash 
#SBATCH -c 1
#SBATCH -t 0-06:00
#SBATCH -p short
#SBATCH --mem=30G
#SBATCH -o /home/sak0914/Errors/zerrors_%j.out 
#SBATCH -e /home/sak0914/Errors/zerrors_%j.err 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skulkarni@g.harvard.edu

# snakemake --configfile config.yaml --cores 8 --rerun-incomplete --keep-going
snakemake --configfile config.yaml --cores 1 --use-conda  --conda-frontend conda --keep-going 