# Mtb Megapipe: Illumina Short-Read WGS Variant Calling Pipeline for <i>Mycobacterium tuberculosis</i>

To run the pipeline:

```
snakemake --configfile config.yaml --cores 1  --keep-going --use-conda --conda-frontend conda
```

To create the DAG of jobs, run either of the following. The first will list every sample, the second will list just the steps. 
```
snakemake --dag | dot -Tpng > dag.png

snakemake --rulegraph | dot -Tpng > dag.png
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


        # find {params.sample_out_dir} -path "*/bam/*.dedup.bam" -type f > {output.run_IDs_file}
        # samtools depth -a -Q 1 -f {output.run_IDs_file} > {output.depth_file}
