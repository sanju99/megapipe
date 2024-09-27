import os, glob
import numpy as np
import pandas as pd

# define some paths to make the path names more readable
sample_out_dir = f"{output_dir}/{{sample_ID}}"
run_out_dir = f"{output_dir}/{{sample_ID}}/{{run_ID}}"

scripts_dir = config["scripts_dir"]
references_dir = config["references_dir"]


rule get_input_FASTQ_files:
    output:
        fastq1 = f"{run_out_dir}/{{run_ID}}_R1.fastq.gz",
        fastq2 = f"{run_out_dir}/{{run_ID}}_R2.fastq.gz",

        fastq1_unzipped = temp(f"{run_out_dir}/{{run_ID}}_1.fastq"),
        fastq2_unzipped = temp(f"{run_out_dir}/{{run_ID}}_2.fastq"),
    params:
        sample_out_dir = sample_out_dir,
        fastq_dir = config["fastq_dir"],
    run:        
        if download_public_FASTQ_dict[wildcards.sample_ID] == 1:
            shell("""
                module load sratoolkit/2.10.7

                # the script deletes the unzipped FASTQ files, so don't need to do it in the rule
                bash scripts/download_FASTQ.sh {params.sample_out_dir} {wildcards.run_ID}
            """)
        elif download_public_FASTQ_dict[wildcards.sample_ID] == 0:
            shell("""
                # copy the FASTQ files from the directory specified in the config file to the sample directory
                # they will be deleted in the next rule after performing adapter trimming, so they won't be doubly stored
                cp {params.fastq_dir}/{wildcards.run_ID}/{wildcards.run_ID}_R1.fastq.gz {output.fastq1}
                cp {params.fastq_dir}/{wildcards.run_ID}/{wildcards.run_ID}_R2.fastq.gz {output.fastq2}

                gunzip -c {output.fastq1} > {output.fastq1_unzipped}
                gunzip -c {output.fastq2} > {output.fastq2_unzipped}

                # first check that the original FASTQ files have the same numbers of lines
                FQ1_line_count=$(wc -l {output.fastq1_unzipped} | awk '{{print $1}}')
                FQ2_line_count=$(wc -l {output.fastq2_unzipped} | awk '{{print $1}}')
                
                # check that neither FASTQ file has no reads
                if [ $FQ1_line_count -eq 0 ] || [ $FQ2_line_count -eq 0 ]; then
                    echo "Error: At least one of the FASTQ files for $sample_ID/$run_ID has no reads"
                    exit 1
                # Compare the counts and raise an error if they are not equal 
                elif [ "$FQ1_line_count" -ne "$FQ2_line_count" ]; then
                    echo "Error: FASTQ files for $sample_ID/$run_ID have different line counts: $FQ1_line_count and $FQ2_line_count"
                    exit 1
                fi
                
                # compare paired end read files. If they are the same, then add to error list. Suppress output with -s tag, so it doesn't print out the differences
                # If the files are identical, the exit status is 0, and the condition is considered true, so an error will be returned.
                if cmp -s {output.fastq1_unzipped} {output.fastq2_unzipped}; then
                   echo "Error: {output.fastq1_unzipped} and {output.fastq2_unzipped} are duplicates"
                   exit 1
                fi
            """)


rule trim_adapters:
    input:
        fastq1 = f"{run_out_dir}/{{run_ID}}_R1.fastq.gz",
        fastq2 = f"{run_out_dir}/{{run_ID}}_R2.fastq.gz",
    output:
        fastq1_trimmed = f"{run_out_dir}/fastp/{{run_ID}}.R1.trimmed.fastq",
        fastq2_trimmed = f"{run_out_dir}/fastp/{{run_ID}}.R2.trimmed.fastq",
        fastp_html = f"{run_out_dir}/fastp/fastp.html",
        fastp_json = f"{run_out_dir}/fastp/fastp.json"
    conda:
        "./envs/read_processing_aln.yaml"
    params:
        min_read_length = config["min_read_length"]
    shell:
        """
        fastp -i {input.fastq1} -I {input.fastq2} -o {output.fastq1_trimmed} -O {output.fastq2_trimmed} -h {output.fastp_html} -j {output.fastp_json} --length_required {params.min_read_length} --dedup --thread 8

        rm {input.fastq1} {input.fastq2}
        """

rule kraken_classification:
    input:
        fastq1_trimmed = f"{run_out_dir}/fastp/{{run_ID}}.R1.trimmed.fastq",
        fastq2_trimmed = f"{run_out_dir}/fastp/{{run_ID}}.R2.trimmed.fastq",
    output:
        fastq1_trimmed_classified = f"{run_out_dir}/kraken/{{run_ID}}_1.kraken.filtered.fastq",
        fastq2_trimmed_classified = f"{run_out_dir}/kraken/{{run_ID}}_2.kraken.filtered.fastq",
        kraken_report = f"{run_out_dir}/kraken/kraken_report",
        kraken_classifications = temp(f"{run_out_dir}/kraken/kraken_classifications"),
    conda:
        "./envs/read_processing_aln.yaml"
    params:
        kraken_db = config["kraken_db"],
        output_dir = output_dir,
        classified_out_string = f"{run_out_dir}/kraken/{{run_ID}}#.kraken.filtered.fastq"
    shell:
        """
        kraken2 --db {params.kraken_db} --threads 8 --paired {input.fastq1_trimmed} {input.fastq2_trimmed} --report {output.kraken_report} --classified-out {params.classified_out_string} > {output.kraken_classifications}
        
        rm {input.fastq1_trimmed} {input.fastq2_trimmed}
        """

# Checkpoint to run Kraken classification
checkpoint kraken_classify:
    input:
        "path/to/input_file"  # Your input data (e.g., FASTQ file)
    output:
        kraken_report="path/to/kraken_report.txt"
    shell:
        """
        # Run Kraken classification and generate report
        kraken2 --db {params.kraken_db} --output {output.kraken_report} {input}
        """

# Function to check if Kraken classification passes the threshold
def check_kraken_classification(wildcards):
    checkpoint_output = checkpoints.kraken_classify.get(**wildcards).output.kraken_report
    
    # Read the unclassified percentage from the Kraken report
    with open(checkpoint_output) as f:
        for line in f:
            if 'unclassified' in line:
                unclassified_percent = float(line.split()[0])
                break
    
    if unclassified_percent > kraken_unclassified_max:
        # Skip the next step if the unclassified percentage is too high
        raise ValueError(f"Unclassified percentage {unclassified_percent}% exceeds the threshold")
    
    # Return the path to the next output if the threshold passes
    return "path/to/next_output_file"

# Rule for the next step, e.g., downstream analysis (only runs if classification passes)
rule downstream_analysis:
    input:
        # Use the function to check the threshold and proceed conditionally
        check_kraken_classification
    output:
        "path/to/next_output_file"
    shell:
        """
        # Perform downstream analysis
        echo "Running downstream analysis on {input}"
        """


rule fastlin_typing:
    input:
        fastq1_trimmed_classified=f"{run_out_dir}/kraken/{{run_ID}}_1.kraken.filtered.fastq",
        fastq2_trimmed_classified=f"{run_out_dir}/kraken/{{run_ID}}_2.kraken.filtered.fastq",
    output:
        fastq1_trimmed_classified_gzipped = temp(f"{run_out_dir}/fastlin/{{run_ID}}_1.fastq.gz"),
        fastq2_trimmed_classified_gzipped = temp(f"{run_out_dir}/fastlin/{{run_ID}}_2.fastq.gz"),
        fastlin_dir = directory(f"{run_out_dir}/fastlin"),
        fastlin_output = f"{run_out_dir}/fastlin/output.txt"
    conda:
        "./envs/read_processing_aln.yaml"
    params:
        fastlin_barcodes = os.path.join(references_dir, "phylogeny", "MTBC_barcodes.tsv"),
    shell:
        """
        gzip -c {input.fastq1_trimmed_classified} > {output.fastq1_trimmed_classified_gzipped}
        gzip -c {input.fastq2_trimmed_classified} > {output.fastq2_trimmed_classified_gzipped}
        
        fastlin -d {output.fastlin_dir} -b {params.fastlin_barcodes} -o {output.fastlin_output} -x 150
        """


def does_sample_pass_fastlin(output_dir, sample_ID):
    """
    Use this function only for cases when there are multiple distinct sequencing runs for a single sample. Before merging the individual BAMs, check that they are assigned the same lineage according to fastlin. This is to mitigate the risk of potentially mislabeled WGS runs being merged together. 
    """    
    fastlin_outputs = pd.concat([pd.read_csv(fName, sep='\t') for fName in glob.glob(f"{output_dir}/{sample_ID}/*/fastlin/output.txt")]).reset_index(drop=True)

    if len(fastlin_outputs) == 0:
        raise ValueError(f"There are no fastlin outputs for {output_dir}/{sample_ID}")

    # if there is only 1 WGS run, continue with that one
    elif len(fastlin_outputs) == 1:
        return True

    # split lineage from median k-mer occurrence
    for i, row in fastlin_outputs.iterrows():
    
        if ',' not in row['lineages']:

            # most common lineage is the same as the only lineage present
            fastlin_outputs.loc[i, ['lineage', 'lineage_koccur', 'most_common_lineage', 'most_common_lineage_koccur']] = [row['lineages'].split(' ')[0], row['lineages'].split(' ')[1].replace('(', '').replace(')', ''), row['lineages'].split(' ')[0], row['lineages'].split(' ')[1].replace('(', '').replace(')', '')]
        
        else:
    
            fastlin_lineage_lst = []
            fastlin_median_occur_lst = []
            
            for single_lineage in row['lineages'].split(', '):
                fastlin_lineage_lst.append(single_lineage.split(' ')[0])
                fastlin_median_occur_lst.append(single_lineage.split(' ')[1].replace('(', '').replace(')', ''))

            # sort by alpha-numeric order
            # need indices to sort the k-mer occurrences too
            alpha_sorted_idx = np.argsort(fastlin_lineage_lst)
            alpha_sorted_koccur_lst = [fastlin_median_occur_lst[idx] for idx in alpha_sorted_idx]

            # convert values in fastlin_median_occur_lst from strings to ints
            fastlin_median_occur_lst = np.array(fastlin_median_occur_lst).astype(int)

            # also add a column for the more common lineage in cases where there are multiple lineages assigned
            most_common_lineage = fastlin_lineage_lst[np.argmax(fastlin_median_occur_lst)]
            
            fastlin_outputs.loc[i, ['lineage', 'lineage_koccur', 'most_common_lineage', 'most_common_lineage_koccur']] = [','.join(np.sort(fastlin_lineage_lst)), ','.join(alpha_sorted_koccur_lst), most_common_lineage, np.max(fastlin_median_occur_lst)]

    fastlin_outputs['most_common_lineage_koccur'] = fastlin_outputs['most_common_lineage_koccur'].astype(int)
    
    # if the most common lineages across the runs match, then it probably indicates low-level contamination by another lineage, not a lineage mixture. We want to keep these
    if fastlin_outputs['most_common_lineage'].nunique() > 1:
        return False

    return True


rule check_fastlin:
    input:
        fastlin_output = f"{run_out_dir}/fastlin/output.txt"
    run:
        if not does_sample_pass_fastlin(params.output_dir, wildcards.sample_ID):
            print(f"Halting pipeline for {wildcards.sample_ID} because the different WGS runs for it have different lineages assigned by fastlin")
            exit()


rule align_reads_mark_duplicates:
    input:
        # require the fastlin output file as an input so that the fastlin rule gets run
        fastlin_output = f"{run_out_dir}/fastlin/output.txt",
        kraken_report = f"{run_out_dir}/kraken/kraken_report",
        fastq1_trimmed_classified=f"{run_out_dir}/kraken/{{run_ID}}_1.kraken.filtered.fastq",
        fastq2_trimmed_classified=f"{run_out_dir}/kraken/{{run_ID}}_2.kraken.filtered.fastq",
    output:
        sam_file = temp(f"{run_out_dir}/bam/{{run_ID}}.sam"),
        bam_file = temp(f"{run_out_dir}/bam/{{run_ID}}.bam"),
        bam_index_file = temp(f"{run_out_dir}/bam/{{run_ID}}.bam.bai"),
        bam_file_dedup = f"{run_out_dir}/bam/{{run_ID}}.dedup.bam",
        bam_file_dedup_metrics = f"{run_out_dir}/bam/{{run_ID}}.dedup.bam.metrics",
        bam_index_file_dedup = f"{run_out_dir}/bam/{{run_ID}}.dedup.bam.bai",
    params:
        output_dir = output_dir,
        kraken_unclassified_max = config["kraken_unclassified_max"],
        ref_genome = os.path.join(references_dir, "ref_genome", "H37Rv_NC_000962.3.fna"),
    conda:
        "./envs/read_processing_aln.yaml"
    shell:
        """
        # only include runs with kraken classified proportion of at least 75%. Otherwise, they are highly contaminated
        unclassified_percent=$(cat {input.kraken_report} | grep unclassified  | awk '{{print $1}}')

        # unclassified_percent can be None if there are 0 unclassified reads (so all the reads map to MTBC)
        # in that case, the if statement below, will fail, so have to check if it's None
        if [ -z "$unclassified_percent" ]; then
            unclassified_percent=0
        fi

        # stop if kraken-classified percentage is too high
        if [ "$(awk 'BEGIN{print ('$unclassified_percent' > {params.kraken_unclassified_max})}')" -eq 1 ]; then
            echo "$unclassified_percent of reads do not map to MTBC, which is below the threshold {params.kraken_unclassified_max} for inclusion. Halting this sample"
            exit
        else
            echo "$unclassified_percent of reads do not map to MTBC, which passes the threshold {params.kraken_unclassified_max} for inclusion"
        fi

        # index reference genome (which is required before aligning reads)
        bwa-mem2 index {params.ref_genome}

        # align reads to the reference genome sequence. The RG name specifies the read group name, which is necessary if you are merging multiple WGS runs into a single BAM file
        bwa-mem2 mem -M -R "@RG\\tID:{wildcards.run_ID}\\tSM:{wildcards.run_ID}" -t 8 {params.ref_genome} {input.fastq1_trimmed_classified} {input.fastq2_trimmed_classified} > {output.sam_file}

        # sort alignment and convert to bam file
        samtools view -b {output.sam_file} | samtools sort > {output.bam_file}

        # index alignment, which creates a .bai index file
        samtools index {output.bam_file}

        # -Xmx6g specifies to allocate 6 GB
        picard -Xmx30g MarkDuplicates I={output.bam_file} O={output.bam_file_dedup} REMOVE_DUPLICATES=true M={output.bam_file_dedup_metrics} ASSUME_SORT_ORDER=coordinate READ_NAME_REGEX='(?:.*.)?([0-9]+)[^.]*.([0-9]+)[^.]*.([0-9]+)[^.]*$'

        # index the deduplicated alignment with samtools, which will create a dedup_bam_file.bai file
        samtools index {output.bam_file_dedup}

        # delete the FASTQ files because they are no longer needed
        rm {input.fastq1_trimmed_classified} {input.fastq2_trimmed_classified}
        """


rule get_BAM_file_depths:
    input:
        bam_file_dedup = lambda wildcards: [f"{output_dir}/{wildcards.sample_ID}/{run_ID}/bam/{run_ID}.dedup.bam" for run_ID in sample_run_dict[wildcards.sample_ID]],
    params:
        ref_genome = os.path.join(references_dir, "ref_genome", "H37Rv_NC_000962.3.fna"),
        sample_out_dir = sample_out_dir,
    output:
        depth_file = f"{sample_out_dir}/bam/{{sample_ID}}.depth.tsv",
    conda:
        "./envs/read_processing_aln.yaml"
    shell:
        """
        # Write all .bam files to a text file

        # get all runs associated with this sample_ID and compute depth
        # -a computes depth at all positions, not just those with non-zero depth
        # -Q is for minimum mapping quality: use 1, so that multiply mapped reads aren't counted. These have mapping quality of 0
        samtools depth -a -Q 1 {input.bam_file_dedup} > {output.depth_file}

        # get the length of the reference genome
        genome_length=$(tail -n +2 {params.ref_genome} | tr -d '\n' | wc -c) # remove first line (FASTA header) and newline characters, then count characters to get ref genome length

        # when there are multiple bam files, each one is its own column in the depth file.
        num_sites_H37Rv=$(wc -l {output.depth_file} | awk '{{print $1}}')
    
        if [ ! "$num_sites_H37Rv" -eq "$genome_length" ]; then
            echo "Check that all $genome_length sites in the H37Rv reference genome are in {output.depth_file}, which currently has $num_sites_H37Rv sites"
            exit 1
        fi
        """


rule get_BAMs_passing_QC_thresholds:
    input:
        depth_file = f"{sample_out_dir}/bam/{{sample_ID}}.depth.tsv", # contains depths for all BAM files for all WGS runs
        bam_file_dedup = lambda wildcards: [f"{output_dir}/{wildcards.sample_ID}/{run_ID}/bam/{run_ID}.dedup.bam" for run_ID in sample_run_dict[wildcards.sample_ID]],
    output:
        pass_BAMs_file = f"{sample_out_dir}/bam/pass_BAMs.txt",
        depth_file_gzip = f"{sample_out_dir}/bam/{{sample_ID}}.depth.tsv.gz",
    params:
        sample_out_dir = sample_out_dir,
        BAM_depth_QC_script = os.path.join(scripts_dir, "BAM_depth_QC.py"),
        median_depth = config["median_depth"],
        min_cov = config["min_cov"],
        genome_cov_prop = config["genome_cov_prop"],
    shell:
        """
        # run the script to determine which runs pass the BAM depth criteria
        python3 -u {params.BAM_depth_QC_script} -i {input.depth_file} -b {input.bam_file_dedup} -o {output.pass_BAMs_file} --median-depth {params.median_depth} --min-cov {params.min_cov} --genome-cov-prop {params.genome_cov_prop}

        # finally, gzip the depth file because it is very large
        gzip {input.depth_file}
        """


rule merge_BAMs:
    input:
        pass_BAMs_file = f"{sample_out_dir}/bam/pass_BAMs.txt",
    output:
        merged_bam_file = f"{sample_out_dir}/bam/{{sample_ID}}.dedup.bam",
        merged_bam_index_file = f"{sample_out_dir}/bam/{{sample_ID}}.dedup.bam.bai",
    conda:
        "./envs/read_processing_aln.yaml"
    params:
        ref_genome = os.path.join(references_dir, "ref_genome", "H37Rv_NC_000962.3.fna"),
        sample_out_dir = sample_out_dir,
        median_depth = config["median_depth"],
        min_cov = config["min_cov"],
        genome_cov_prop = config["genome_cov_prop"],
    shell:
        """
        num_runs_passed=$(wc -l {input.pass_BAMs_file} | awk '{{print $1}}')

        # stop processing samples that don't pass the BAM coverage requirements
        if [ $num_runs_passed -eq 0 ]; then
            echo "No BAM files for {wildcards.sample_ID} passed the minimum coverage requirements. Halting pipeline for this sample"
            exit
            
        else 
            # if only one BAM file passed, or there is only one sequencing run for this isolate, just use that BAM file for variant calling
            echo "$num_runs_passed WGS runs for {wildcards.sample_ID} have median depth ≥ {params.median_depth} and at least {params.genome_cov_prop} of sites with {params.min_cov}x coverage"

            # merge them using samtools. works because the original bam files were sorted prior to running picard and dropping duplicates (after which they remain sorted)
            samtools merge -b {input.pass_BAMs_file} {output.merged_bam_file}

            if [ $num_runs_passed -eq 1 ]; then

                # delete the original BAM file to reduce disk space usage because it's a duplicate of the merged BAM file
                for file_path in $(cat {input.pass_BAMs_file}); do
                    rm "$file_path" "$file_path.bai"
                done

            fi

            # index the merged BAM file for variant calling
            samtools index {output.merged_bam_file}

        fi
        """


rule pilon_variant_calling:
    input:
        merged_bam_file = f"{sample_out_dir}/bam/{{sample_ID}}.dedup.bam",
    output:
        vcf_file = temp(f"{sample_out_dir}/pilon/{{sample_ID}}.vcf"),
        vcf_file_gzip = f"{sample_out_dir}/pilon/{{sample_ID}}_full.vcf.gz",
        vcf_file_variants_only = f"{sample_out_dir}/pilon/{{sample_ID}}_variants.vcf",
        fasta_file = temp(f"{sample_out_dir}/pilon/{{sample_ID}}.fasta"),        
    params:
        ref_genome = os.path.join(references_dir, "ref_genome", "H37Rv_NC_000962.3.fna"),
        sample_pilon_dir = f"{sample_out_dir}/pilon",
    conda:
        "./envs/variant_calling.yaml"
    shell:
        """
        pilon -Xmx30g --minmq 1 --genome {params.ref_genome} --bam {input.merged_bam_file} --output {wildcards.sample_ID} --outdir {params.sample_pilon_dir} --variant
            
        # then gzip the full VCF file and delete the unzipped version. Also delete the FASTA file because it's not needed
        gzip -c {output.vcf_file} > {output.vcf_file_gzip}

        # save the variants only (non-REF calls) to another VCF file
        bcftools view --types snps,indels,mnps,other {output.vcf_file_gzip} > {output.vcf_file_variants_only}
        """



# rule freebayes_variant_calling:
#     input:
#         merged_bam_file = f"{sample_out_dir}/bam/{{sample_ID}}.dedup.bam",
#     output:
#         vcf_file = temp(f"{sample_out_dir}/freebayes/{{sample_ID}}.vcf"),
#         vcf_file_gzip = f"{sample_out_dir}/freebayes/{{sample_ID}}_full.vcf.gz",
#         vcf_file_SNPs_only = f"{sample_out_dir}/freebayes/{{sample_ID}}_SNPs.vcf",
#     params:
#         ref_genome = os.path.join(references_dir, "ref_genome", "H37Rv_NC_000962.3.fna"),
#     conda:
#         "./envs/variant_calling.yaml"
#     shell:
#         """
#         # use ploidy of 1
#         freebayes -f {params.ref_genome} -p 1 {input.merged_bam_file} > {output.vcf_file}

#         # then gzip the full VCF file and delete the unzipped version. Also delete the FASTA file because it's not needed
#         gzip -c {output.vcf_file} > {output.vcf_file_gzip}

#         # save SNPs only to another VCF file
#         bcftools view --types snps,mnps {output.vcf_file_gzip} > {output.vcf_file_SNPs_only}
#         """



rule create_lineage_helper_files:
    input:
        vcf_file_gzip = f"{sample_out_dir}/pilon/{{sample_ID}}_full.vcf.gz",
    params:
        lineage_pos_for_F2 = os.path.join(references_dir, "phylogeny", "Coll2014_positions_all.txt"),
        output_dir = output_dir,
    output:
        bcf_file = f"{sample_out_dir}/lineage/{{sample_ID}}.bcf",
        bcf_index_file = f"{sample_out_dir}/lineage/{{sample_ID}}.bcf.csi",
        vcf_lineage_positions = f"{sample_out_dir}/lineage/{{sample_ID}}_lineage_positions.vcf",
    conda:
        "./envs/variant_calling.yaml"
    shell:
        """
        # convert the full VCF file to a BCF fileto get only the lineage-defining positions according to the Coll 2014 scheme
        bcftools view {input.vcf_file_gzip} -O b -o {output.bcf_file}

        # index bcf file
        bcftools index {output.bcf_file}

        # create VCF file of just the lineage positions, which will be used by the F2 metric script. Per the documentation, if --regions-file is a tab-delimited file, then it needs two columns (CHROM and POS), and POS is 1-indexed and inclusive
        # THIS IS DIFFERENT BEHAVIOR FROM IF IT WAS A BED FILE OR IF YOU USE BEDTOOLS. IN BOTH OF THOSE CASES, YOU NEED THREE COLUMNS (CHROM, BEG, AND END), AND THEY ARE 0-INDEXED WITH END BEING EXCLUSIVE (I.E. HALF-OPEN)
        bcftools view {output.bcf_file} --regions-file {params.lineage_pos_for_F2} -O v -o {output.vcf_lineage_positions}   
        """


rule lineage_typing:
    input:
        bcf_file = f"{sample_out_dir}/lineage/{{sample_ID}}.bcf",
        bcf_index_file = f"{sample_out_dir}/lineage/{{sample_ID}}.bcf.csi",
        vcf_lineage_positions = f"{sample_out_dir}/lineage/{{sample_ID}}_lineage_positions.vcf",
        vcf_file_variants_only = f"{sample_out_dir}/pilon/{{sample_ID}}_variants.vcf",
    params:
        lineage_SNP_info = os.path.join(references_dir, "phylogeny", "Coll2014_SNPs_all.csv"),
        F2_metric_script = os.path.join(scripts_dir, "calculate_F2_metric.py"),
        output_dir = output_dir,        
    output:
        F2_metric_output = f"{sample_out_dir}/lineage/F2_Coll2014.txt",
        fast_lineage_caller_output = f"{sample_out_dir}/lineage/fast_lineage_caller_output.txt",
    shell:
        """
        python3 -u {params.F2_metric_script} -i {params.output_dir}/{wildcards.sample_ID} -o {output.F2_metric_output} --lineage-file {params.lineage_SNP_info}

        rm {input.bcf_file} {input.bcf_index_file} {input.vcf_lineage_positions}

        fast-lineage-caller {input.vcf_file_variants_only} --pass --out {output.fast_lineage_caller_output}
        """


rule combine_codon_variants:
    input:
        vcf_file_variants_only = f"{sample_out_dir}/pilon/{{sample_ID}}_variants.vcf",
    output:
        vcf_file_variants_combinedCodons = f"{sample_out_dir}/pilon/{{sample_ID}}_variants_combinedCodons.vcf",
    params:
        combine_codon_variants_script = os.path.join(scripts_dir, "combine_codon_variants.py"),
    shell:
        """
        # this creates _variants_combined_codons.vcf in the same directory as the input VCF. You can also specify the -o flag if you want to write the output file somewhere else
        python3 -u {params.combine_codon_variants_script} -i {input.vcf_file_variants_only}        
        """



rule annotate_variants_snpEff:
    input:
        vcf_file_variants_combinedCodons = f"{sample_out_dir}/pilon/{{sample_ID}}_variants_combinedCodons.vcf",
    output:
        vcf_file_variants_combinedCodons_annot = f"{sample_out_dir}/WHO_resistance/{{sample_ID}}_variants_combinedCodons.eff.vcf",
    conda:
        "./envs/variant_annotation.yaml"
    params:
        snpEff_db = config['snpEff_db'],
    shell:
        """
        snpEff eff {params.snpEff_db} -noStats -no-downstream -no-upstream -lof {input.vcf_file_variants_combinedCodons} > {output.vcf_file_variants_combinedCodons_annot}

        rm {input.vcf_file_variants_combinedCodons}
        """



rule create_WHO_catalog_variants_TSV:
    input:
        vcf_file_variants_combinedCodons_annot = f"{sample_out_dir}/WHO_resistance/{{sample_ID}}_variants_combinedCodons.eff.vcf",
    output:
        vcf_file_bgzip = temp(f"{sample_out_dir}/WHO_resistance/{{sample_ID}}_variants_combinedCodons.eff.vcf.bgz"),
        vcf_file_bgzip_tbi = temp(f"{sample_out_dir}/WHO_resistance/{{sample_ID}}_variants_combinedCodons.eff.vcf.bgz.tbi"),
        variants_file_tsv = f"{sample_out_dir}/WHO_resistance/{{sample_ID}}_variants.tsv",
    params:
        WHO_catalog_regions_BED_file = os.path.join(references_dir, "WHO_catalog_resistance", "regions.bed"),
    conda:
        "./envs/variant_annotation.yaml"
    shell:
        """
        # need to bgzip the VCF file to use bcftools view with the region argument. NEED TO PUT "" AROUND FILE NAME TO PROPERLY CONSIDER SPECIAL CHARACTERS IN FILENAME
        bgzip -c {input.vcf_file_variants_combinedCodons_annot} > {output.vcf_file_bgzip}
    
        # tabix the bgzipped file, which will create fName.bgz.tbi
        tabix -p vcf -f {output.vcf_file_bgzip}

        bcftools view -R {params.WHO_catalog_regions_BED_file} {output.vcf_file_bgzip} | SnpSift extractFields '-' POS REF ALT FILTER QUAL IMPRECISE AF DP BQ MQ IC DC ANN -e "" > {output.variants_file_tsv}
        """


rule get_WHO_catalog_resistance_predictions:
    input:
        variants_file_tsv = f"{sample_out_dir}/WHO_resistance/{{sample_ID}}_variants.tsv",
    output:
        variants_file_tsv_annot = f"{sample_out_dir}/WHO_resistance/{{sample_ID}}_variants_annot.tsv",
        predictions_fName = f"{sample_out_dir}/WHO_resistance/{{sample_ID}}_pred_AF_thresh_75.csv",
        predictions_fName_lowAF = f"{sample_out_dir}/WHO_resistance/{{sample_ID}}_pred_AF_thresh_25.csv",
    params:
        output_file_basename = f"{sample_out_dir}/WHO_resistance/{{sample_ID}}_pred",
        process_variants_WHO_catalog_script = os.path.join(scripts_dir, "process_variants_for_WHO_catalog.py"),
        get_WHO_resistance_predictions_script = os.path.join(scripts_dir, "WHO_catalog_resistance_pred.py"),
    shell:
        """
        python3 -u {params.process_variants_WHO_catalog_script} -i {input.variants_file_tsv}
        rm {input.variants_file_tsv}

        # get resistance predictions -- any Group 1 or 2 variant that passes QC leads to a prediction of R for a given drug. If not, predicted S
        python3 -u {params.get_WHO_resistance_predictions_script} -i {output.variants_file_tsv_annot} -o {params.output_file_basename}
        python3 -u {params.get_WHO_resistance_predictions_script} -i {output.variants_file_tsv_annot} -o {params.output_file_basename} --AF-thresh 0.25
        """


rule get_SNPs_for_phylogenetic_tree:
    input:
        vcf_file_gzip = f"{sample_out_dir}/pilon/{{sample_ID}}_full.vcf.gz",
    output:
        vcf_SNP_sites = temp(f"{sample_out_dir}/lineage/SNP_sites.tsv"),
        vcf_SNP_sites_gzip = f"{sample_out_dir}/lineage/SNP_sites.tsv.gz",
    conda:
        "./envs/variant_annotation.yaml"
    params:
        exclude_regions_BED_file = os.path.join(references_dir, "phylogeny", "exclude_regions.bed"),
    shell:
        """
        # exclude structural variants (SVTYPE)
        gunzip -c {input.vcf_file_gzip} | bcftools filter -e "SVTYPE == 'INS' | SVTYPE == 'DEL'" | bedtools subtract -a '-' -b {params.exclude_regions_BED_file} | SnpSift extractFields '-' POS REF ALT FILTER QUAL AF DP BQ MQ -e "" > {output.vcf_SNP_sites}

        gzip -c {output.vcf_SNP_sites} > {output.vcf_SNP_sites_gzip}
        """