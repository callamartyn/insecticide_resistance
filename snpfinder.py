import vcf

import sys

import pandas as pd

from Bio import SeqIO

import numpy as np

# getting a list of the vcf files to analyze
vcf_files = sys.argv[2:]
# getting the name of the reference, will be used to select relevant gene positions
ref_name = sys.argv[1].rstrip('.fasta')
# import fasta sequence and convert to string
ref_seq = str(SeqIO.read(sys.argv[1], "fasta").seq)

# importing the csv file containing all known resistance snps for all references
full_table = pd.read_csv('all_vgsg_snps.csv')
# subsetting just the snps corresponding to the reference used
snp_table = full_table[(full_table['reference'] == ref_name)]
# creating a list of tuples from the table, to be used in making a dictionary
snp_list = list(zip(snp_table.position, snp_table.name))
# adding positions for +1 and +2 nucleotides for each snp, to capture the entire codon
codon_snps = []
for snp in snp_list:
    codon_snps.append(((snp[0] + 1), snp[1]))
    codon_snps.append(((snp[0] + 2), snp[1]))
# adding the new positions to the original list
snp_list = snp_list + codon_snps
# sorting the list by nuclotide position
snp_list.sort(key=lambda tup: tup[0])
# creating a dictionary from the list of snps, key = nucleotide position value = snp name
snp_dict = dict(snp_list)

# importing translation table of possible codons and corresponding amino acids
codon_table = pd.read_csv('codon_table.csv')
# converting csv file into dictionary (key = codon, value = aa)
codon_dict = dict(zip(codon_table.codon, codon_table.amino_acid))

def compare(a, b):
    if codon_dict[a] == codon_dict[b]:
        return ['synonymous', codon_dict[a], codon_dict[b]]
    else:
        return ['non-synonymous', codon_dict[a], codon_dict[b]]

# creating an empty list for snps found in the dataset, snps from each vcf will be added here
sample_snps = []
# iterating through each vcf file and combining all of the entries into a single list
for file in vcf_files:
    vcf_reader = vcf.Reader(open(file, 'r'))
    # removing prefix and suffix from vcf file path to get the sample name
    sample = file.rstrip('.vcf')
    sample = file.lstrip('./')
    for record in vcf_reader:
        # defining the nucleotide position, reference base, and alternate base using vcf_reader
        POS = record.POS
        REF = record.REF
        ALT = str(record.ALT).strip('/[]')
        QUAL = record.QUAL
        INFO = record.INFO
        RO = INFO['RO']
        AO = INFO['AO'][0]
        # snp at codon position 1
        if len(REF) > 1:
            continue
        if POS % 3 == 1:
            start = POS - 1
            end = start + 3
            ref_codon = (ref_seq[start:end])
            alt_codon = ALT + ref_codon[1:3]
        # snp at codon position 2
        elif POS % 3 == 2:
            start = POS - 2
            end = start + 3
            ref_codon = (ref_seq[start:end])
            alt_codon = ref_codon[0] + ALT + ref_codon[2]
        # snp at codon position 3
        elif POS % 3 == 0:
            start = POS - 3
            end = start + 3
            ref_codon = (ref_seq[start:end])
            alt_codon = ref_codon[0:2] + ALT
        snp_type = compare(ref_codon, alt_codon)
        # first will check if the snp is in the list of known snps using the dictionary
        if POS in snp_dict:
            # if it is, including the name of the snp (codon) in the output
            codon = snp_dict[POS]
        else:
            # if it is not speficying that this snp is not reported as a resitance snp
            codon = 'not reported'
        sample_snps.append([sample, codon, REF, ALT, POS, QUAL, snp_type[0], snp_type[1], snp_type[2], AO, RO])

# turning the list from all the samples into a pandas dataframe
results = pd.DataFrame(sample_snps, columns = ['Sample', 'Codon', 'Ref_Base', 'Alt_Base', 'Position', 'Quality', 'Change_Type', 'Ref_AA', 'Alt_AA', 'Alt_Reads', 'Ref_reads'])
# exporting the dataframe as a csv file in the working directory
results.to_csv('results.csv')
# also printing it out for instant gratification
print(results)

# creating new dataframe with only the non-synonymous snps
non_synonymous = pd.DataFrame(data = (results[(results['Change_Type'] == 'non-synonymous')]))
non_synonymous.reset_index(drop=True, inplace=True)
dataset_snps = non_synonymous.Position.unique()
new_list = []
for snp in dataset_snps:
    sample_group = []
    quality_scores = []
    ref_reads = []
    alt_reads = []
    for i in range(len(non_synonymous.Position)):
        if(non_synonymous.Position[i] == snp):
            sample_group.append(non_synonymous.Sample[i])
            quality_scores.append(non_synonymous.Quality[i])
            alt_reads.append(non_synonymous.Alt_Reads[i])
            ref_base = non_synonymous.Ref_Base[i]
            alt_base = non_synonymous.Alt_Base[i]
            ref_aa = non_synonymous.Ref_AA[i]
            alt_aa = non_synonymous.Alt_AA[i]
            codon = non_synonymous.Codon[i]
    number = len(sample_group)
    ave_qual = np.mean(quality_scores)
    ave_reads = np.mean(alt_reads)
    snp_row = [snp, codon, number, ref_base, alt_base, ref_aa, alt_aa, ave_qual, ave_reads, sample_group]
    new_list.append(snp_row)
non_synonymous_table = pd.DataFrame(new_list, columns = ['SNP', 'Codon', 'No_Occurrences', 'Ref_Base', 'Alt_Base', 'Ref_AA', 'Alt_AA', 'Average_Quality', 'Average_no_reads', 'Samples_found_in'])
non_synonymous_table.sort_values(by='No_Occurrences', ascending=False, inplace=True)
non_synonymous_table.to_csv('non_synonymous_table.csv')
print(non_synonymous_table)
