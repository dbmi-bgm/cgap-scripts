#!/usr/bin/env python3

################################################
#
#   Add samplegeno information to variants
#       in VCF file
#
#   Michele Berselli
#   berselli.michele@gmail.com
#
################################################

################################################
#
#   Libraries
#       requires granite library
#
################################################
import sys, argparse, os
#import vcf_parser from granite
from granite.lib import vcf_parser

################################################
#
#   Functions
#
################################################
################################################
#   runner
################################################
def main(args):
    ''' '''
    # Variables
    granite_def = '##GRANITE=<ID=SAMPLEGENO>'
    samplegeno_def = '##INFO=<ID=SAMPLEGENO,Number=.,Type=String,Description="Sample genotype information. Subembedded:\'samplegeno\':Format:\'NUMGT|GT|AD|SAMPLEID\'">'

    # Buffers
    fo = open(args['outputfile'], 'w')

    # Creating Vcf object
    vcf_obj = vcf_parser.Vcf(args['inputfile'])

    # Update and write header
    vcf_obj.header.add_tag_definition(granite_def + '\n' + samplegeno_def, 'INFO')
    vcf_obj.write_header(fo)

    # Reading variants and adding samplegeno
    for vnt_obj in vcf_obj.parse_variants():
        # Init empty samplegeno
        samplegeno = []
        vnt_obj.complete_genotype()
        # Get possible alleles
        # REF is at idx 0, ALT idxs are maintained
        alleles = [vnt_obj.REF] + vnt_obj.ALT.split(',')

        # Calculate samplegeno
        for ID_genotype in vnt_obj.IDs_genotypes:
            samplegeno_ = []
            GT_0_, GT_1_ = vnt_obj.get_genotype_value(ID_genotype, 'GT').replace('|', '/').split('/')
            if GT_0_ != '.': GT_0 = alleles[int(GT_0_)]
            else: GT_0 = GT_0_
            #end if
            if GT_1_ != '.': GT_1 = alleles[int(GT_1_)]
            else: GT_1 = GT_1_
            #end if
            AD = vnt_obj.get_genotype_value(ID_genotype, 'AD').replace(',', '/')

            samplegeno_.append(GT_0_ + '/' +  GT_1_)
            samplegeno_.append(GT_0 + '/' + GT_1)
            samplegeno_.append(AD)
            samplegeno_.append(ID_genotype)

            samplegeno.append('|'.join(samplegeno_))
        #end for

        # Add samplegeno to variant INFO
        vnt_obj.add_tag_info('SAMPLEGENO={0}'.format(','.join(samplegeno)))

        # Write variant
        vcf_obj.write_variant(fo, vnt_obj)
    #end for

    fo.close()
#end def main

################################################
#
#   MAIN
#
################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Add samplegeno information to variants in VCF file. SAMPLEGENO tag added to INFO field')

    parser.add_argument('-i','--inputfile', help='input VCF file', required=True)
    parser.add_argument('-o','--outputfile', help='output VCF file', required=True)

    args = vars(parser.parse_args())

    main(args)

#end if
