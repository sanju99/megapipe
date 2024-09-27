# Mtb Megapipe: Illumina Short-Read WGS Variant Calling Pipeline for <i>Mycobacterium tuberculosis</i>

To run the pipeline:

```
snakemake --configfile config.yaml --cores 1 --keep-going --use-conda --conda-frontend conda
```

Add the flag `--rerun-incomplete` if a snakemake workflow was interrupted and there are incomplete files.

To create the DAG of jobs, run either of the following. The first will list every sample, the second will list just the steps. 

```bash
snakemake --dag | dot -Tsvg > dag.svg

snakemake --rulegraph | dot -Tsvg > dag.svg
```

<!-- rule repair_reads_bbmap:
    input:
        fastq1=f"{run_out_dir}/fastq/{{run_ID}}_R1.fastq.gz",
        fastq2=f"{run_out_dir}/fastq/{{run_ID}}_R2.fastq.gz"
    output:
        fastq1_fixed=f"{run_out_dir}/fastq/{{run_ID}}.R1.fixed.fastq",
        fastq2_fixed=f"{run_out_dir}/fastq/{{run_ID}}.R2.fixed.fastq",
    conda:
        "./envs/bioinformatics.yaml"
    shell:
        """
        bash $CONDA_PREFIX/bin/repair.sh in={input.fastq1} in2={input.fastq2} out={output.fastq1_fixed} out2={output.fastq2_fixed}
        """ -->

```bash
# get all human-adapted MTBC taxids
esearch -db taxonomy -query "txid1773[Subtree]" | efetch -format uid > ./references/phylogeny/human_MTBC_taxids.txt
```