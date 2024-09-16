import pandas as pd
import numpy as np
import glob, os, warnings, argparse
warnings.filterwarnings("ignore")

import Bio.SeqUtils
import Bio.Data
from Bio import Seq, SeqIO, Entrez

# dataframe of all regions in H37Rv and their coordinates. Includes genes and non-coding regions
h37Rv_full = pd.read_csv("./references/ref_genome/mycobrowser_h37rv_v4.csv")

################################################################## NOTES ##################################################################

# initiator_codon_variant and stop_retained_variant are encoded differently in snpEff than in the catalog, but they're not Group 1-2, so ignore
# silent variants composed of MNVs are formatted like this in the catalog: katG_c.CTCinsGAT, whereas the combineCodons function outputs katG_c.CTC>GAT

# this includes LoF variants -- i.e. Rv0678_LoF is significantly associated with Bedaquiline resistance, so all component LoF mutations (frameshift, stop gained, start lost) are Group 2) Assoc w R - Interim
who_catalog = pd.read_csv("./references/WHO_catalog_resistance/V2_catalog.csv", header=[2])
who_catalog = who_catalog[['drug', 'variant', 'FINAL CONFIDENCE GRADING', 'Comment']].rename(columns={'FINAL CONFIDENCE GRADING': 'confidence'})

parser = argparse.ArgumentParser()

# Add a required string argument for the input file (annotations)
parser.add_argument("-i", "--input", dest='input_file', type=str, required=True)

# Add a required string argument for the output file (CSV of resistance predictions)
parser.add_argument("-o", "--output", dest='output_file_basename', type=str, required=True)

parser.add_argument('--AF-thresh', type=float, dest='AF_thresh', default=0.75, help='Alternative allele frequency threshold (exclusive) to consider variants present. Default = 0.75')

cmd_line_args = parser.parse_args()

input_file = cmd_line_args.input_file
output_file_basename = cmd_line_args.output_file_basename
AF_thresh = cmd_line_args.AF_thresh

# must be a float in the range [0, 1)
if AF_thresh > 1:
    AF_thresh /= 100



def get_all_drug_resistances(df):

    resistance_dict = {}
    
    for drug in np.sort(who_catalog['drug'].unique()):

        # get only Group 1-2 resistance mutations
        R_assoc_variants = who_catalog.loc[who_catalog['confidence'].str.contains('Assoc w R')].query("drug==@drug").variant.values

        # for BDQ/CFZ, LoF in mmpL5 abrogates the effect of R-assoc Rv0678 mutations
        # for AMK/KAN, LoF in eis abrogates the effect of R-assoc eis promoter mutations
        # checked that the only mutations with "Abrogates" in the comments are these mutations
        negating_muts = who_catalog.dropna(subset="Comment").query("drug==@drug & Comment.str.contains('Abrogates')").variant.values

        if len(df.query("variant in @negating_muts")) > 0:
            resistance_dict[drug] = 'S'
        else:
            if len(df.query("variant in @R_assoc_variants")) > 0:
                resistance_dict[drug] = 'R'
            else:
                resistance_dict[drug] = 'S'
    
    resistance_df = pd.DataFrame(resistance_dict, index=[0]).T.reset_index()
    resistance_df.columns = ['Drug', 'Phenotype']
    return resistance_df
    


def make_susceptible_dataframe():
    '''
    Make a dataframe with S for all 15 drugs
    '''
    
    resistance_df = pd.DataFrame(columns=['Drug', 'Phenotype'])
    resistance_df['Drug'] = np.sort(who_catalog['drug'].unique())
    resistance_df['Phenotype'] = 'S'

    return resistance_df




# add a suffix for the AF threshold
output_file = output_file_basename + f"_AF_thresh_{int(AF_thresh*100)}.csv"
# '.'.join(output_file.split(".")[:-1]) + f"_AF_thresh_{int(AF_thresh*100)}.csv"

if os.path.isfile(output_file):
    exit()

if 'annot' not in input_file:
    input_file = input_file.replace('.tsv', '_annot.tsv')

# read in processed variants file
df = pd.read_csv(input_file, sep='\t')

# the variants_annot.tsv file is created for snakemake compatibility reasons, but it may be empty, indicating that there are no non-IMPRECISE variants in the WHO catalog regions
if len(df) == 0:
    print(f"{input_file} is all drug-susceptible")
    resistance_df = make_susceptible_dataframe()
    resistance_df.to_csv(output_file, index=False)
    exit()

# some issues with indels: indels with AF in the range (0.25, 0.75] are often present at a high level, but reads that start or end near the indel position don't cover it sufficiently and won't be considered as supporting the indel. These will pull down the AF, even though the indel is present in many reads. So set the AF thresold at 0.25 for these and additionally check that at least 5 reads support the indel (IC or DC).
# these variants typically have high MQ and DP, so they are very likely real

# this is quite prevalent in ethA and to a lesser degree, other non-essential drug resistance genes like katG, pncA, Rv0678, mmpL5, and mmpS5

# check that the complement of what you want is not true because not all variant types have all fields filled in, so you will lose some if you do i.e. AF > 0.75
# for indels (including gene deletions), only check DP, AF, MQ, and FILTER because BQ and QUAL are often NaN for them
df_highConf_muts = pd.concat([df.loc[(df['REF'].str.len() > df['ALT'].str.len()) & 
                                       (~df['FILTER'].str.contains('|'.join(['Del', 'LowCov'])))
                                        ].query("~(AF <= 0.25) & ~(DP < 5) & ~(MQ < 30) & ~(DC < 5)"), # deletion
    
                                df.loc[(df['REF'].str.len() < df['ALT'].str.len()) & 
                                       (~df['FILTER'].str.contains('|'.join(['Del', 'LowCov'])))
                                        ].query("~(AF <= 0.25) & ~(DP < 5) & ~(MQ < 30) & ~(IC < 5)"), # insertion

                                # for SNPs/MNPs, check everything
                                df.loc[(df['REF'].str.len() == df['ALT'].str.len()) & 
                                       (~df['FILTER'].str.contains('|'.join(['Del', 'LowCov'])))
                                        ].query("~(QUAL < 10) & ~(AF <= @AF_thresh) & ~(DP < 5) & ~(BQ < 20) & ~(MQ < 30)") # SNP
                             ])

resistance_df = get_all_drug_resistances(df_highConf_muts)
resistance_df.to_csv(output_file, index=False)