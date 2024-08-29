import numpy as np
import pandas as pd

configfile: "config.yaml"

# Define the list of samples, either directly or through config
output_dir = config["output_dir"]

# dataframe of isolates to run. First column is the sample ID, second column is a comma-separated string of the sequencing run IDs
df_samples_runs = pd.read_csv(config['isolates_to_run'], sep='\t', header=None)
sample_run_dict = dict(zip(df_samples_runs[0], df_samples_runs[1].str.split(',')))

# make a dictionary of binary values indicating whether (1) or not (0) to download publicly available FASTQ files
assert df_samples_runs[2].isin([0, 1]).all()
download_public_FASTQ_dict = dict(zip(df_samples_runs[0], df_samples_runs[2]))

include: "rules.smk"

# Define rules in the Snakefile

# rule all:
#     input:
#         [f"{config['output_dir']}/{config['sample_ID']}/{run_ID}/lineage/fastlin_output.txt" for run_ID in config["run_IDs"]],
rule all:
    input:
        # [f"{output_dir}/{sample_ID}/{run_ID}/{run_ID}_R{rep}.fastq.gz" for sample_ID in sample_run_dict.keys() for run_ID in sample_run_dict[sample_ID] for rep in [1, 2]]
        # [f"{output_dir}/{sample_ID}/bam/pass_run_IDs.txt" for sample_ID in sample_run_dict.keys()],
        # [f"{output_dir}/{sample_ID}/{run_ID}/bam/{run_ID}.dedup.bam" for sample_ID in sample_run_dict.keys() for run_ID in sample_run_dict[sample_ID]],
        # [f"{output_dir}/{sample_ID}/bam/{sample_ID}.dedup.bam" for sample_ID in sample_run_dict.keys()],
        # [f"{output_dir}/{sample_ID}/pilon/{sample_ID}_variants_combinedCodons.vcf" for sample_ID in sample_run_dict.keys()],
        [f"{output_dir}/{sample_ID}/pilon/{sample_ID}_variants_combinedCodons.eff.vcf" for sample_ID in sample_run_dict.keys()],
        [f"{output_dir}/{sample_ID}/lineage/{sample_ID}/F2_Coll2014.txt" for sample_ID in sample_run_dict.keys()]