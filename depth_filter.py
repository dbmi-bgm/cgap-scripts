#!/usr/bin/env python3

################################################
#
#  Script to filter on DP field of genotype(s)
#  to remove variants that do not meet a
#  certain depth threshold

################################################

################################################
#   Libraries
################################################
from granite.lib import vcf_parser
import argparse, subprocess

################################################
#   Functions
################################################

def main(args):
    in_vcf = vcf_parser.Vcf(args['inputSampleVCF'])

    min_depth_to_keep = int(args['min_depth'])

    with open(args['outputfile'], 'w') as fo:
        in_vcf.write_header(fo)
        for vnt_obj in in_vcf.parse_variants():
            sample_list = vnt_obj.IDs_genotypes
            for sample in sample_list:
                if vnt_obj.get_genotype_value(sample, "DP") != ".":
                    if int(vnt_obj.get_genotype_value(sample, "DP")) >= min_depth_to_keep:
                        in_vcf.write_variant(fo, vnt_obj)
                        break

    subprocess.run(["bgzip", args['outputfile']])
    subprocess.run(["tabix", args['outputfile']+".gz"])

################################################
#   Main
################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Filter variants based on depth (DP)')

    parser.add_argument('-i','--inputSampleVCF', help='input sample VCF file', required=True)
    parser.add_argument('-o','--outputfile', help='output file name for filtered VCF', required=True)
    parser.add_argument('-d','--min_depth', help='minimum DP (depth) value for a given variant. At least 1 sample must have DP >= min_depth to keep a variant', required=True)

    args = vars(parser.parse_args())

    main(args)
