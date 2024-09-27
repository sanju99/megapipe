#!/bin/bash 
#SBATCH -c 8
#SBATCH -t 0-11:59
#SBATCH -p short
#SBATCH --mem=50G
#SBATCH -o /home/sak0914/Errors/zerrors_%j.out 
#SBATCH -e /home/sak0914/Errors/zerrors_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skulkarni@g.harvard.edu

snakemake --snakefile snakefile_2 --use-conda --conda-frontend conda --rerun-incomplete --keep-going --config isolates_to_run="TRUST_samples_quick.tsv" --cores 8 --directory snakemake_small