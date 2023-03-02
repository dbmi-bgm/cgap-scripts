#!/usr/bin/env python3

################################################
#
#   Michele Berselli
#   berselli.michele@gmail.com
#
################################################

################################################
#
#   Libraries
#
################################################
import sys, argparse, os
import glob
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from granite.lib import vcf_parser



################################################
#   Functions
################################################

def get_counts(inputfile, repeat):
    """
    """
    # Open input vcf
    vcf = vcf_parser.Vcf(inputfile)
    ID_genotype = vcf.header.IDs_genotypes[0]

    # Parse and get repeat counts
    for vnt_obj in vcf.parse_variants():
        if vnt_obj.get_tag_value('REPID') == repeat:
            # Reference allele
            REF = vnt_obj.get_tag_value('REF')
            # Alternate alleles
            ALT = vnt_obj.ALT.replace('<STR', '').replace('>', '').split(',')
            # ref + alt as integer
            REF_ALT = list(map(int, [REF] + ALT))
            # Get genotype
            GT_0, GT_1 = map(int, vnt_obj.get_genotype_value(ID_genotype, 'GT').replace('|', '/').split('/'))
            # Map genotype to alleles
            #   and return repeat counts
            return int(REF), [REF_ALT[GT_0], REF_ALT[GT_1]]

def plot(x, filename, repeat, REF):
    """
    """
    # Create Histogram
    _, ax = plt.subplots(1, 1)

    # Data to numpy array
    x = np.array(x)
    unique, counts = np.unique(x, return_counts=True)
    unique_ = list(map(str, unique))
    colors = ["#8856a7" if i == REF else "#1c9099" for i in unique]

    ax.bar(unique_, counts, color=colors)

    # Axis and title
    ax.set_title(f'{repeat} Expansion')
    ax.set_ylabel('Observed Alleles in Samples')
    ax.set_xlabel('Number of Repeats (per Allele)')
    ax.set_xticks(unique_)

    # Legend
    pop_a = mpatches.Patch(color='#1c9099', label='Alternate Allele')
    pop_b = mpatches.Patch(color='#8856a7', label='Reference Allele')
    plt.legend(handles=[pop_a, pop_b])

    # Save
    plt.savefig(filename, format='png')
    # plt.show()

def main(args):

    # Data Points
    x = []
    repeat = args['repeat']

    # Initialize output file
    fo = open(args['outputname']+'.tsv', 'w')
    fo.write(f'#Repeat: {repeat}\n')

    # Loop files in inputdir
    inputdir = args['inputdir']
    for i, fn in enumerate(glob.glob(f'{inputdir}/*.vcf')):
        # Get repeat counts for file
        REF, counts = get_counts(fn, repeat)
        # Add counts and write to file
        x += counts
        fn_ = fn.split('/')[-1]
        if i == 0:
            fo.write(f'#Reference: {REF}\n')
            fo.write(f'#File\tAllele_1\tAllele_2\n')
        fo.write(f'{fn_}\t{counts[0]}\t{counts[1]}\n')

    # Close output file
    fo.close()

    if args['histogram']:
        plot(x, args['outputname']+'_histogram.png', repeat, REF)



################################################
#   Main
################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Loop through ExpansionHunter output VCF files, report and plot repeat counts for the target STR.')

    parser.add_argument('-i','--inputdir', help='path to directory containing VCF files', required=True)
    parser.add_argument('-r','--repeat', help='LocusId as defined in STR catalog', required=True)
    parser.add_argument('-o','--outputname', help='basename for output (w/o extension)', required=True)
    parser.add_argument('-p','--histogram', help='plot histogram', action='store_true', required=False)

    args = vars(parser.parse_args())

    main(args)
