import os, glob
import numpy as np
import pandas as pd

# define some paths to reduce string redundancy
sample_out_dir = f"{output_dir}/{{sample_ID}}"
run_out_dir = f"{output_dir}/{{sample_ID}}/{{run_ID}}"

scripts_dir = config["scripts_dir"]
references_dir = config["references_dir"]


rule get_input_FASTQ_files:
    output:
        fastq1 = f"{run_out_dir}/{{run_ID}}_R1.fastq.gz",
        fastq2 = f"{run_out_dir}/{{run_ID}}_R2.fastq.gz",
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
        "./envs/bioinformatics.yaml"
    params:
        min_length = config["min_length"]
    shell:
        """
        fastp -i {input.fastq1} -I {input.fastq2} -o {output.fastq1_trimmed} -O {output.fastq2_trimmed} -h {output.fastp_html} -j {output.fastp_json} --length_required {params.min_length} --dedup --thread 8

        # rm {input.fastq1} {input.fastq2}
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
        "./envs/bioinformatics.yaml"
    params:
        kraken_db = config["kraken_db"],
        output_dir = output_dir,
        classified_out_string = f"{run_out_dir}/kraken/{{run_ID}}#.kraken.filtered.fastq"
    shell:
        """
        kraken2 --db {params.kraken_db} --threads 8 --paired {input.fastq1_trimmed} {input.fastq2_trimmed} --report {output.kraken_report} --classified-out {params.classified_out_string} > {output.kraken_classifications}
        
        rm {input.fastq1_trimmed} {input.fastq2_trimmed}
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
        "./envs/bioinformatics.yaml"
    params:
        fastlin_barcodes = os.path.join(references_dir, "MTBC_barcodes.tsv"),
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
            raise ValueError(f"Halting pipeline for {wildcards.sample_ID} because the different WGS runs for it have different lineages assigned by fastlin")


rule align_reads_mark_duplicates:
    input:
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
        ref_genome = config["ref_genome"],
    conda:
        "./envs/bioinformatics.yaml"
    shell:
        """
        # index reference genome (which is required before aligning reads)
        bwa index {params.ref_genome}

        # align reads to the reference genome sequence. The RG name specifies the read group name, which is necessary if you are merging multiple WGS runs into a single BAM file
        bwa mem -M -R "@RG\\tID:{{run_ID}}\\tSM:{{run_ID}}" -t 8 {params.ref_genome} {input.fastq1_trimmed_classified} {input.fastq2_trimmed_classified} > {output.sam_file}

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


rule get_BAM_depth:
    input:
        bam_file_dedup = lambda wildcards: [f"{output_dir}/{wildcards.sample_ID}/{run_ID}/bam/{run_ID}.dedup.bam" for run_ID in sample_run_dict[wildcards.sample_ID]],
    params:
        ref_genome = config["ref_genome"],
        sample_out_dir = sample_out_dir,
    output:
        # run_IDs_file = f"{sample_out_dir}/bam/run_IDs.txt",
        depth_file = f"{sample_out_dir}/bam/{{sample_ID}}.depth.tsv",
    conda:
        "./envs/bioinformatics.yaml"
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


rule compute_BAMs_passing_QC_thresholds:
    input:
        combined_depth_file = f"{sample_out_dir}/bam/{{sample_ID}}.depth.tsv", # contains depths for all BAM files for all WGS runs
        bam_file_dedup = lambda wildcards: [f"{output_dir}/{wildcards.sample_ID}/{run_ID}/bam/{run_ID}.dedup.bam" for run_ID in sample_run_dict[wildcards.sample_ID]],
    output:
        pass_run_IDs_file = f"{sample_out_dir}/bam/pass_run_IDs.txt",
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
        python3 -u {params.BAM_depth_QC_script} -i {input.combined_depth_file} -b {input.bam_file_dedup} --median-depth {params.median_depth} --min-cov {params.min_cov} --genome-cov-prop {params.genome_cov_prop}

        # finally, gzip the depth file because it is very large
        gzip {input.combined_depth_file}
        """


rule merge_BAMs:
    input:
        pass_run_IDs_file = f"{sample_out_dir}/bam/pass_run_IDs.txt",
    output:
        merged_bam_file = f"{sample_out_dir}/bam/{{sample_ID}}.dedup.bam",
        merged_bam_index_file = f"{sample_out_dir}/bam/{{sample_ID}}.dedup.bam.bai",
    conda:
        "./envs/bioinformatics.yaml"
    params:
        ref_genome = config["ref_genome"],
        sample_out_dir = sample_out_dir,
        median_depth = config["median_depth"],
        min_cov = config["min_cov"],
        genome_cov_prop = config["genome_cov_prop"],
    shell:
        """
        num_runs_passed=$(wc -l {input.pass_run_IDs_file} | awk '{{print $1}}')

        # stop processing samples that don't pass the BAM coverage requirements
        if [ $num_runs_passed -eq 0 ]; then
            echo "No BAM files for {wildcards.sample_ID} passed the minimum coverage requirements. Halting pipeline for this sample"
            exit
            
        else 
            # if only one BAM file passed, or there is only one sequencing run for this isolate, just use that BAM file for variant calling
            echo "$num_runs_passed WGS runs for {wildcards.sample_ID} have median depth > {params.median_depth} and at least {params.genome_cov_prop} of sites covered {params.min_cov} times"

            # merge them using samtools. works because the original bam files were sorted prior to running picard and dropping duplicates (after which they remain sorted)
            samtools merge -b {input.pass_run_IDs_file} {output.merged_bam_file}

            if [ $num_runs_passed -eq 1 ]; then

                # delete the original BAM file to reduce disk space usage because it's a duplicate of the merged BAM file
                for file_path in $(cat {input.pass_run_IDs_file}); do
                    rm $file_path
                done

            fi

            # index the merged BAM file for variant calling
            samtools index {output.merged_bam_file}

        fi
        """

rule variant_calling:
    input:
        merged_bam_file = f"{sample_out_dir}/bam/{{sample_ID}}.dedup.bam",
    output:
        vcf_file = temp(f"{sample_out_dir}/pilon/{{sample_ID}}.vcf"),
        vcf_file_gzip = f"{sample_out_dir}/pilon/{{sample_ID}}_full.vcf.gz",
        vcf_file_variants_only = f"{sample_out_dir}/pilon/{{sample_ID}}_variants.vcf",
        # fasta_file = temp(f"{sample_out_dir}/pilon/{{sample_ID}}.fasta"),        
    params:
        ref_genome = config["ref_genome"],
        sample_pilon_dir = f"{sample_out_dir}/pilon",
    conda:
        "./envs/bioinformatics.yaml"
    shell:
        """
        pilon -Xmx30g --minmq 1 --genome {params.ref_genome} --bam {input.merged_bam_file} --output {wildcards.sample_ID} --outdir {params.sample_pilon_dir} --variant
            
        # then gzip the full VCF file and delete the unzipped version. Also delete the FASTA file because it's not needed
        gzip -c {output.vcf_file} > {output.vcf_file_gzip}

        # save the variants only (non-REF calls) to another VCF file
        bcftools view --types snps,indels,mnps,other {output.vcf_file_gzip} > {output.vcf_file_variants_only}
        """

rule get_sample_lineage_and_F2:
    input:
        vcf_file_gzip = f"{sample_out_dir}/pilon/{{sample_ID}}_full.vcf.gz",
        vcf_file_variants_only = f"{sample_out_dir}/pilon/{{sample_ID}}_variants.vcf",
    params:
        lineage_pos_for_F2 = os.path.join(references_dir, "Coll2014_positions_all.txt"),
        F2_metric_script = os.path.join(scripts_dir, "calculate_F2_metric.py"),
        output_dir = output_dir,
    output:
        bcf_file = temp(f"{sample_out_dir}/lineage/{{sample_ID}}.bcf"),
        bcf_index_file = temp(f"{sample_out_dir}/lineage/{{sample_ID}}.bcf.csi"),
        vcf_lineage_positions = temp(f"{sample_out_dir}/lineage/{{sample_ID}}_lineage_positions.vcf"),
        fast_lineage_caller_output = f"{sample_out_dir}/lineage/fast_lineage_caller_output.txt",
        F2_metric_output = f"{sample_out_dir}/lineage/{{sample_ID}}/F2_Coll2014.txt",
    conda:
        "./envs/bioinformatics.yaml"
    shell:
        """
        # convert the full VCF file to a BCF fileto get only the lineage-defining positions according to the Coll 2014 scheme
        bcftools view {input.vcf_file_gzip} -O b -o {output.bcf_file}

        # index bcf file
        bcftools index {output.bcf_file}

        # create VCF file of just the lineage positions, which will be used by the F2 metric script. Per the documentation, if --regions-file is a tab-delimited file, then it needs two columns (CHROM and POS), and POS is 1-indexed and inclusive
        # THIS IS DIFFERENT BEHAVIOR FROM IF IT WAS A BED FILE OR IF YOU USE BEDTOOLS. IN BOTH OF THOSE CASES, YOU NEED THREE COLUMNS (CHROM, BEG, AND END), AND THEY ARE 0-INDEXED WITH END BEING EXCLUSIVE (I.E. HALF-OPEN)
        bcftools view {output.bcf_file} --regions-file {params.lineage_pos_for_F2} -O v -o {output.vcf_lineage_positions}   
        
        fast-lineage-caller {output.vcf_file_variants_only} --pass --out {output.fast_lineage_caller_output}

        python3 -u {params.F2_metric_script} -i f"{params.output_dir}/{wildcards.sample_ID}"
        """


rule combine_codon_variants:
    input:
        vcf_file_variants_only = f"{sample_out_dir}/pilon/{{sample_ID}}_variants.vcf",
    output:
        vcf_file_variants_combinedCodons = f"{sample_out_dir}/pilon/{{sample_ID}}_variants_combinedCodons.vcf",
    params:
        combine_codon_variants_script = os.path.join(scripts_dir, "combine_codon_variants.py"),
        annot_fName = config["annot_VCF_name"],
    run:
        shell(
            """
            # this creates _variants_combined_codons.vcf in the same directory as the input VCF. You can also specify the -o flag if you want to write the output file somewhere else
            python3 -u {params.combine_codon_variants_script} -i {input.vcf_file_variants_only}            
            """
        )

        # then write the file name to the text file to keep track of files to be annotated with snpEff in the next rule
        with open({params.annot_fName}, 'w+') as file:
            for vcf_fName in {output.vcf_file_variants_combinedCodons}:
                file.write(vcf_fName + "\n")


rule annotate_variants_snpEff:
    input:
        vcf_file_variants_combinedCodons = f"{sample_out_dir}/pilon/{{sample_ID}}_variants_combinedCodons.vcf",
    output:
        vcf_file_variants_combinedCodons_annot = f"{sample_out_dir}/pilon/{{sample_ID}}_variants_combinedCodons.eff.vcf",
    conda:
        "./envs/bioinformatics.yaml"
    params:
        annot_fName = config["annot_VCF_name"],
        snpEff_db = config['snpEff_db'],
    shell:
        """
        snpEff eff {params.snpEff_db} -noStats -no-downstream -no-upstream -lof -fileList {params.annot_fName}

        # then delete the text file so that VCFs don't get re-annotated in later runs
        rm {params.annot_fName}
        """


rule write_WHO_catalog_variants_TSV:
    input:
        vcf_file_variants_combinedCodons_annot = f"{sample_out_dir}/pilon/{{sample_ID}}_variants_combinedCodons.eff.vcf",
    output:
        vcf_file_bgzip = temp(f"{sample_out_dir}/pilon/{{sample_ID}}_variants_combinedCodons.eff.vcf.bgz"),
        vcf_file_bgzip_tbi = temp(f"{sample_out_dir}/pilon/{{sample_ID}}_variants_combinedCodons.eff.vcf.bgz.tbi"),
        variants_file_tsv = f"{sample_out_dir}/pilon/{{sample_ID}}_variants.tsv",
    params:
        WHO_catalog_regions_BED_file = os.path.join(references_dir, "WHO_catalog_regions.bed"),
    conda:
        "./envs/bioinformatics.yaml"
    shell:
        """
        # need to bgzip the VCF file to use bcftools view with the region argument. NEED TO PUT "" AROUND FILE NAME TO PROPERLY CONSIDER SPECIAL CHARACTERS IN FILENAME
        bgzip -c {inpput.vcf_file_variants_combinedCodons_annot} > {output.vcf_file_bgzip}
    
        # tabix the bgzipped file, which will create fName.bgz.tbi
        tabix -0 -p vcf {output.vcf_file_bgzip} -f

        bcftools view -R {params.WHO_catalog_regions_BED_file} {output.vcf_file_bgzip} | SnpSift extractFields '-' POS REF ALT FILTER QUAL IMPRECISE AF DP BQ MQ IC DC ANN -e "" > {output.variants_file_tsv}
        """

rule get_WHO_catalog_resistance_predictions:
    input:
        variants_file_tsv = f"{sample_out_dir}/pilon/{{sample_ID}}_variants.tsv",
    output:
        variants_file_tsv_annot = f"{sample_out_dir}/pilon/{{sample_ID}}_variants_annot.tsv",
        out_fName="$out_dir/$sample_ID/WHO_resistance/${sample_ID}_pred.csv",
    params:
        process_variants_WHO_catalog_script = os.path.join(scripts_dir, "process_variants_WHO_catalog_script.py"),
        get_WHO_resistance_predictions_script = os.path.join(scripts_dir, "WHO_catalog_resistance_pred.py"),
    shell:
        """
        python3 -u {params.process_variants_WHO_catalog_script} -i {input.variants_file_tsv}

        # get resistance predictions -- any Group 1 or 2 variant that passes QC leads to a prediction of R for a given drug. If not, predicted S
        python3 -u {params.get_WHO_resistance_predictions_script} -i {output.variants_file_tsv_annot} -o $out_fName
        python3 -u {params.get_WHO_resistance_predictions_script} -i {output.variants_file_tsv_annot} -o $out_fName --AF-thresh 0.25
        """