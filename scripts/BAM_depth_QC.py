import numpy as np
import pandas as pd
import sys, os, argparse

parser = argparse.ArgumentParser()

# Add a required string argument for the config file
parser.add_argument("-i", dest='depths_fName', type=str, required=True)

parser.add_argument("-b", dest="BAM_files_lst", type=list, require=True)

parser.add_argument('--median-depth', dest='median_depth_thresh', default=15, type=int, help='Minimum median depth (exclusive) that an alignment must meet to be included')

parser.add_argument('--min-cov', dest='min_cov', default=20, type=int, help='Minimum number of reads (inclusive) covering a site')

parser.add_argument('--genome-cov-prop', dest='genome_cov_prop', default=0.90, type=float, help='Minimum proportion of the genome (inclusive) that must be covered')

cmd_line_args = parser.parse_args()

depths_fName = cmd_line_args.depths_fName
BAM_files_lst = cmd_line_args.BAM_files_lst
median_depth_thresh = cmd_line_args.median_depth_thresh
min_cov = cmd_line_args.min_cov
genome_cov_prop = cmd_line_args.genome_cov_prop

# sample_ID = os.path.basename(sample_dir)
depths = pd.read_csv(depths_fName, sep='\t', header=None)
print(depths.head())

# write output files to the same directory as the depth file
sample_dir = os.path.dirname(depths)
print(sample_dir)

if len(BAM_files_lst) + 2 != depths.shape[1]:
    raise ValueError(f"Number of columns in {depths_fName} is not consistent with {len(BAM_files_lst)} sequencing runs")

if depths.shape[1] < 3:
    raise ValueError(f"here should be at least 3 columns in {depths_fName}. There are only {depths.shape[1]}")

depths.columns = ['CHROM', 'POS'] + BAM_files_lst

with open(f"{sample_dir}/bam/pass_run_IDs.txt", "w+") as file:
    
    for BAM_file in BAM_files_lst:
    
        # median depth across the entire H37Rv ref genome
        median = depths[BAM_file].median()
    
        # proportion of sites with a coverage of at least 20. Round in case there are samples with 0.949 or something (saw one with 0.9498)
        prop_sites_cov_thresh = np.round(len(depths.loc[depths[BAM_file] >= min_cov]) / len(depths), 2)

        # if the BAM file passes the thresholds, write it to the pass_run_IDs file to be used in the merge BAMS rule
        if median > median_depth_thresh and prop_sites_cov_thresh >= genome_cov_prop:

            file.write(BAM_file + "\n")