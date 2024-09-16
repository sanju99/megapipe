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

To make the custom kraken database (assuming yours is named 20240915_humanAdaptedMTBC_DB), follow these instructions:

```bash
DB_NAME="20240915_humanAdaptedMTBC_DB"

kraken2-build --download-taxonomy --db $DB_NAME

# use FTP to avoid rsync connection problems to NCBI
kraken2-build --download-library bacteria --db $DB_NAME --use-ftp
```

To add your own genomes to the database, they must have a taxid in the FASTA file header, so add the human-adapted MTBC taxid (1773)

```python
for fName in glob.glob("./kraken_db/genomes/*.fasta"):

    # added these manually by searching their exact taxids
    if not os.path.basename(fName).startswith('GCF'):
   
        sequences = [(seq.id, seq.seq) for seq in SeqIO.parse(fName, "fasta")]
        assert len(sequences) == 1
            
        # make the header the sample name
        header = os.path.basename(fName).split('.')[0]

        # add |kraken:taxid|1773 to the header for each one
        header += '|kraken:taxid|1773'
    
        with open(fName, "w") as file:
            file.write(f">{header}\n")
            file.write(str(sequences[0][1]) + "\n")
```

```bash
for file in "./kraken_db/genomes/*.fasta"
do
    kraken2-build --add-to-library $file --db $DBNAME
done
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