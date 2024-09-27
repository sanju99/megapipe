To make the custom kraken database (assuming yours is named 20240915_humanAdaptedMTBC_DB), follow these instructions, taken from the <a href="https://github.com/DerrickWood/kraken2/wiki/Manual" target="_blank">kraken2 documentaion</a>:

```bash
DBNAME="custom_human_MTBC"

kraken2-build --download-taxonomy --db $DBNAME
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

For the four reference genomes, I manully appended the taxids to each header. 

```bash
# Add all 155 genomes in the `genomes/` directory to the database
find ./genomes/ -name '*.fasta' -print0 | xargs -0 -I{} -n1 kraken2-build --add-to-library {} --db $DBNAME

# build the database. Took approximately 5 minutes for 155 genomes
kraken2-build --build --db $DBNAME

# verify that the database was built correctly
kraken2-inspect --db $DBNAME

# clean the database to reduce disk usage
kraken2-build --clean --db $DBNAME
```