#!/bin/bash 
#SBATCH -c 1
#SBATCH -t 0-11:59
#SBATCH -p short 
#SBATCH --mem=50G
#SBATCH -o /home/sak0914/Errors/zerrors_%j.out 
#SBATCH -e /home/sak0914/Errors/zerrors_%j.err 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skulkarni@g.harvard.edu


#################### REQUIRED COMMAND LINE ARGUMENTS (IN THIS ORDER): ####################

# 1. TSV file with 2 columns (AND NO HEADER): 
    # col1: sample_ID, i.e. "SAMEA104362063"
    # col2: comma-separated string of all sequencing runs for the sample in col1, i.e. "ERR2184222,ERR9029914"
# 2. fastq_dir: directory where paired-end Illumina sequencing reads are stored
# 3. out_dir: output directory where results should be stored. i.e. /n/data1/hms/dbmi/farhat/rollingDB/cryptic_output or /n/data1/hms/dbmi/farhat/rollingDB/genomic_data

##### NOTE: The values in col1 and col2 of the input file can be the same, i.e. for unpublished data. But for standardization, it is best to use the BioSample accession ID in col1 and sequencing run IDs in col2 #####

set -o errexit # any error will cause the shell script to exit immediately. This is not native bash behavior
source activate bioinformatics # CHANGE TO YOUR OWN ENVIRONMENT OR REMOVE IF RUNNING IN THE BASE ENV

if ! [ $# -eq 3 ]; then
    echo "Please pass in 3 command line arguments: a text file with two columns: sample_ID and sequencing runs, the FASTQ directory, and the output directory"
    exit
fi

# command line arguments
input_file=$1
fastq_dir=$2
out_dir=$3 # directory where completed variant calling results for sample_ID IDs will be stored (OUTPUT)

repair_script="/home/sak0914/anaconda3/envs/bioinformatics/bin/repair.sh" # CHANGE TO WHATEVER FILE PATH IS IN YOUR ENVIRONMENT
ref_genome="/n/data1/hms/dbmi/farhat/mtb_data/h37rv/h37rv.fna"
genome_length=$(tail -n +2 $ref_genome | tr -d '\n' | wc -c) # remove first line (FASTA header) and newline characters, then count characters to get ref genome length

min_length=50 # parameter determined by Max. Drop trimmed reads that are shorter than this length

#################### IF YOU USE THE STANDARD DATABASE, THEN YOU HAVE TO RUN EXTRACT_KRAKEN_READS.PY TO GET ONLY THOSE THAT MAP TO MTBC AND CHILDREN TAXA ####################
#################### IF YOU USE THE MTBC DATABASE, THEN YOU JUST HAVE TO EXTRACT THE CLASSIFIED READS BECAUSE THEY ARE THE ONLY ONES THAT GET CLASSIFIED TO MTBC -- DO THIS!!! ####################

# custom MTBC database
kraken_db="/n/data1/hms/dbmi/farhat/mm774/References/Kraken2_DB_Dir/Kraken2_DB"
redo_file="./redo_files.txt"
annot_fName="./fNames_for_annot.txt"

# if the file doesn't exist, create empty file, so wc -l will return 0
touch $annot_fName

# empty it (which is necessary if the file already existed before running touch)
truncate -s 0 $annot_fName

# Check if the output directory exists. If not, raise an error
if [ ! -d "$out_dir" ]; then
    echo "Output directory $out_dir doesn't exist!"
    exit 1
fi

# define function to check FASTQ file quality
function check_FASTQ_quality() {
    
    local run_ID="$1"

    # FASTQ files for the single sequencing run
    fastq1_file="$fastq_dir/$run_ID/${run_ID}_R1.fastq.gz"
    fastq2_file="$fastq_dir/$run_ID/${run_ID}_R2.fastq.gz"
            
    # unzip FASTQ files to compare
    FQ1_unzipped="${fastq1_file/.gz/}"
    FQ2_unzipped="${fastq2_file/.gz/}"

    gunzip -c $fastq1_file > $FQ1_unzipped
    gunzip -c $fastq2_file > $FQ2_unzipped

    # first check that the original FASTQ files have the same numbers of lines
    FQ1_line_count=$(wc -l $FQ1_unzipped | awk '{print $1}')
    FQ2_line_count=$(wc -l $FQ2_unzipped | awk '{print $1}')
    
    # check that neither FASTQ file has no reads
    if [ $FQ1_line_count -eq 0 ] || [ $FQ2_line_count -eq 0 ]; then
        echo "Error: $fastq1_file or $fastq2_file has no reads"
        exit 1
    # Compare the counts and raise an error if they are not equal 
    elif [ "$FQ1_line_count" -ne "$FQ2_line_count" ]; then
        echo "Error: $fastq1_file and $fastq2_file have different line counts: $FQ1_line_count and $FQ2_line_count"
        exit 1
    else
        echo "Line counts in paired-end FASTQ files for $run_ID match"
    fi

    # compare paired end read files. If they are the same, then add to error list. Suppress output with -s tag, so it doesn't print out the differences
    # If the files are identical, it returns an exit status of 0, and the condition is considered true, so an error will be returned.
    if cmp -s "$FQ1_unzipped" "$FQ2_unzipped"; then
       echo "Error: $fastq1_file and $fastq2_file are duplicates"
       exit 1
    fi

    # then delete the unzipped files
    rm $FQ1_unzipped
    rm $FQ2_unzipped
}



function make_subdirs() {

    local run_out_dir="$1"
    
    # bam will contain the deduplicated bam files, and pilon will contain the FASTA file and gzipped VCF file
    run_fastp_dir="$run_out_dir/fastp"
    run_kraken_dir="$run_out_dir/kraken"
    run_bam_dir="$run_out_dir/bam"
    
    subdirs_lst=($run_fastp_dir $run_kraken_dir $run_bam_dir)
    
    # Iterate through the array of subdirectories and create each one if it doesn't exist
    # the -p flag creates all nested directories if they don't exist
    for i in ${!subdirs_lst[@]}; do
        if [ ! -d "${subdirs_lst[$i]}" ]; then
          mkdir -p "${subdirs_lst[$i]}"
        fi
    done

}



# define function to perform the trimming, kraken-classification, and aligning steps
function repair_trim_classify_reads() {

    local run_ID="$1"
    local run_out_dir="$2"

    local run_fastp_dir="$run_out_dir/fastp"
    local run_kraken_dir="$run_out_dir/kraken"

    # FASTQ files for the single sequencing run
    fastq1_file="$fastq_dir/$run_ID/${run_ID}_R1.fastq.gz"
    fastq2_file="$fastq_dir/$run_ID/${run_ID}_R2.fastq.gz"

    # this is a rare occurrence, but is included for all samples
    # in some FASTQ files, reads are duplicated in one paired end file, so then when bwa tries to align them, it gets a read mismatch
    # it won't get fixed just by sorting the reads because there are multiple of a single read
    # use the repair script from bbmap (conda-installable)
    if [ ! -f "$run_out_dir/$run_ID.R2.fixed.fastq" ]; then
        /bin/bash $repair_script in="$fastq1_file" in2="$fastq2_file" out="$run_out_dir/$run_ID.R1.fixed.fastq" out2="$run_out_dir/$run_ID.R2.fixed.fastq"
    fi    
    
    # trim reads with the minimum allowable read length (post-trimming). This discards reads that are below the minimum length after trimming
    # default is to discard reads in the individual R1/R2 files if only one read in a pair passes QC. 
    # also perform deduplication, as a precaution. Duplicate reads can cause issues with bwa alignment later
    if [ ! -f "$run_out_dir/$run_ID.R2.fastq" ]; then
        fastp -i "$run_out_dir/$run_ID.R1.fixed.fastq" -I "$run_out_dir/$run_ID.R2.fixed.fastq" -o "$run_out_dir/$run_ID.R1.fastq" -O "$run_out_dir/$run_ID.R2.fastq" -h "$run_fastp_dir/fastp.html" -j "$run_fastp_dir/fastp.json" --length_required $min_length --dedup --thread 8
    fi

    # the classified-out flag creates output files with the format following the flag, where # is replaced with _1 or _2 for the paired end reads
    # also pipe the read classification information (which, by default, is printed to stdout) to an output file
    if [ ! -f "$run_out_dir/${run_ID}.R_2.kraken.filtered.fastq" ]; then
        kraken2 --db $kraken_db --threads 8 --paired "$run_out_dir/$run_ID.R1.fastq" "$run_out_dir/$run_ID.R2.fastq" --report "$run_kraken_dir/kraken_report_standard_db" > "$run_kraken_dir/kraken_read_classifications_standard_db"
    fi

}

# Read the TSV file line by line, skiping the header. IFS sets the field separator. Here, it is tab
while IFS=$'\t', read -r sample_ID run_IDs_string
do

    # create a directory for each sample to be processed. The sample_ID is the name of the directory 
    # in rollingDB/genomic_data or rollingDB/cryptic_output, each sample has a folder named with the sample name
    # bam dir contain the deduplicated bam files, pilon dir will contain the VCF files, and F2 dir will contain the .txt file with the F2 metric 
    sample_out_dir="$out_dir/$sample_ID"

    # # make top-level directory first
    # if [ ! -d "$sample_out_dir" ]; then
    #     mkdir "$sample_out_dir"
    # fi

    sample_bam_dir="$sample_out_dir/bam"
    sample_pilon_dir="$sample_out_dir/pilon"
    sample_lineage_dir="$sample_out_dir/lineage_calling"

    # # Iterate through the array of subdirectories and create each one if it doesn't exist
    # subdirs_lst=($sample_bam_dir $sample_pilon_dir $sample_lineage_dir)
    
    # # the -p flag creates all nested directories if they don't exist
    # for i in ${!subdirs_lst[@]}; do
    #     if [ ! -d "${subdirs_lst[$i]}" ]; then
    #       mkdir -p "${subdirs_lst[$i]}"
    #     fi
    # done

    # if they are the same, don't make an extra directory level for sequencing runs
    if [[ "$sample_ID" == "$run_IDs_string" ]]; then

        # make subdirectories for the sequencing run
        # make_subdirs $sample_out_dir

        if [ ! -f "$sample_out_dir/kraken/kraken_report_standard_db" ]; then

            # first check FASTQ quality using the function. The function will return an error if conditions are not met
            # check_FASTQ_quality $sample_ID
            repair_trim_classify_reads $sample_ID $sample_out_dir

        fi
        
        # delete all files that are not in the subdirectories for each sequencing run. these files are not critical, so it saves space
        find $sample_out_dir -maxdepth 1 -type f -delete

    else
    
        # create an array of the sequencing runs associated with a single isolate
        IFS=',' read -ra run_IDs_array <<< "$run_IDs_string"
        
        # within each sample directory, create a directory for all individual sequencing runs
        for run_ID in "${run_IDs_array[@]}"; do

            # make top-level directory first
            run_out_dir="$sample_out_dir/$run_ID"
            
            # if [ ! -d "$run_out_dir" ]; then
            #     mkdir "$run_out_dir"
            # fi

            # make subdirectories for each sequencing run
            # make_subdirs $run_out_dir

            # perform trimming, kraken-classification, and read alignment for each sequencing run
            if [ ! -f "$run_out_dir/kraken/kraken_report_standard_db" ]; then

                # first check FASTQ quality using the function. The function will return an error if conditions are not met
                # check_FASTQ_quality $run_ID
                repair_trim_classify_reads $run_ID $run_out_dir
            fi
            
            # delete all files that are not in the subdirectories for each sequencing run. these files are not critical, so it saves space
            find $run_out_dir -maxdepth 1 -type f -delete

        done
    fi

done < "$input_file"