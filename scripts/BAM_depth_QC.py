import numpy as np
import pandas as pd
import sys, os, argparse

parser = argparse.ArgumentParser()

# Add a required string argument for the config file
parser.add_argument("-i", dest='depths_fName', type=str, required=True, help='Input gzipped tab-separated file of read depths across an alignment')

# allow multiple sctrings for the BAM files list argument so that any number can be passed in
parser.add_argument("-b", '--bam-files', dest="BAM_files_lst", nargs='+', type=str, required=True, help='List of BAM files associated with this sample')

parser.add_argument("-o", dest='output_fName', type=str, required=True, help='Output text file to write the BAM files that pass the QC thresholds to')

parser.add_argument('--median-depth', dest='median_depth_thresh', default=15, type=int, help='Minimum median depth (exclusive) that an alignment must meet to be included')

parser.add_argument('--min-cov', dest='min_cov', default=20, type=int, help='Minimum number of reads (inclusive) covering a site')

parser.add_argument('--genome-cov-prop', dest='genome_cov_prop', default=0.90, type=float, help='Minimum proportion of the genome (inclusive) that must be covered')

cmd_line_args = parser.parse_args()

depths_fName = cmd_line_args.depths_fName
BAM_files_lst = cmd_line_args.BAM_files_lst
output_fName = cmd_line_args.output_fName
median_depth_thresh = cmd_line_args.median_depth_thresh
min_cov = cmd_line_args.min_cov
genome_cov_prop = cmd_line_args.genome_cov_prop

depths = pd.read_csv(depths_fName, sep='\t', header=None)

# write output files to the same directory as the depth file
sample_bam_dir = os.path.dirname(depths_fName)

# the first two columns are CHROM and POS, so there should be 2 + N columns in the depths dataframe, where N is the number of runs
if len(BAM_files_lst) + 2 != depths.shape[1]:
    raise ValueError(f"Number of columns in {depths_fName} is not consistent with {len(BAM_files_lst)} sequencing runs")

if depths.shape[1] < 3:
    raise ValueError(f"here should be at least 3 columns in {depths_fName}. There are only {depths.shape[1]}")

depths.columns = ['CHROM', 'POS'] + BAM_files_lst

with open(output_fName, "w+") as file:
    
    for BAM_file in BAM_files_lst:
    
        # median depth across the entire H37Rv ref genome
        median = depths[BAM_file].median()
    
        # proportion of sites with a coverage of at least 20. Round in case there are samples with 0.949 or something (saw one with 0.9498)
        prop_sites_cov_thresh = np.round(len(depths.loc[depths[BAM_file] >= min_cov]) / len(depths), 2)

        # if the BAM file passes the thresholds, write it to the pass_run_IDs file to be used in the merge BAMS rule
        if median >= median_depth_thresh and prop_sites_cov_thresh >= genome_cov_prop:

            file.write(BAM_file + "\n")