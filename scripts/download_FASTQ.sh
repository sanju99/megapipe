set -o errexit # any error will cause the shell script to exit immediately. This is not native bash behavior

if ! [ $# -eq 2 ]; then
    echo "Please pass in 2 command line arguments: the sample output directory and the WGS run"
    exit
fi

sample_out_dir=$1
run_ID=$2
    
# Create the directory if it doesn't exist
if [ ! -d "$sample_out_dir/$run_ID" ]; then
    mkdir -p "$sample_out_dir/$run_ID"
fi

# Download the FASTQ files
fastq-dump --split-files --outdir "$sample_out_dir/$run_ID" $run_ID

FQ1_fName="$sample_out_dir/$run_ID/${run_ID}_1.fastq"
FQ2_fName="$sample_out_dir/$run_ID/${run_ID}_2.fastq"

# first check that the original FASTQ files have the same numbers of lines
FQ1_line_count=$(wc -l $FQ1_fName | awk '{print $1}')
FQ2_line_count=$(wc -l $FQ2_fName | awk '{print $1}')

# check that neither FASTQ file has no reads
if [ $FQ1_line_count -eq 0 ] || [ $FQ2_line_count -eq 0 ]; then
    echo "Error: At least one of the FASTQ files for $sample_ID/$run_ID has no reads"
    exit 1
# Compare the counts and raise an error if they are not equal 
elif [ "$FQ1_line_count" -ne "$FQ2_line_count" ]; then
    echo "Error: FASTQ files for $sample_ID/$run_ID have different line counts: $FQ1_line_count and $FQ2_line_count"
    exit 1
# else
#     echo "Line counts in paired-end FASTQ files for $sample_ID/$run_ID match"
fi

# compare paired end read files. If they are the same, then add to error list. Suppress output with -s tag, so it doesn't print out the differences
# If the files are identical, the exit status is 0, and the condition is considered true, so an error will be returned.
if cmp -s "$FQ1_fName" "$FQ2_fName"; then
   echo "Error: $FQ1_fName and $FQ2_fName are duplicates"
   exit 1
fi

# Gzip the FASTQ files
gzip -c $FQ1_fName > "$sample_out_dir/$run_ID/${run_ID}_R1.fastq.gz"
gzip -c $FQ2_fName > "$sample_out_dir/$run_ID/${run_ID}_R2.fastq.gz"

# # delete the unzipped files
# rm $FQ1_fName
# rm $FQ2_fName