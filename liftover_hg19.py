#!/usr/bin/env python3

################################################
#
#    Script to add hg19 lifover coordinates
#   into the INFO column of sample vcf files
#  Note: Should be run on filtered variant list
#
################################################

################################################
#   Libraries
################################################

from granite.lib import vcf_parser
from pyliftover import LiftOver
import sys, argparse, subprocess

################################################
#   Functions
################################################

def main(args):
    # open input vcf
    vcf = vcf_parser.Vcf(args['inputfile'])
    # add 3 new tag definitions - for hg19 liftover: chr, pos, and end
    hg19CHROM_definition = '##INFO=<ID=hg19_chr,Number=1,Type=String,Description="CHROM in hg19 using LiftOver from pyliftover">'
    hg19POS_definition = '##INFO=<ID=hg19_pos,Number=1,Type=Integer,Description="POS in hg19 using LiftOver from pyliftover (converted back to 1-based)">'
    hg19END_definition = '##INFO=<ID=hg19_end,Number=1,Type=Integer,Description="END in hg19 using LiftOver from pyliftover (converted back to 1-based)">'
    vcf.header.add_tag_definition(hg19END_definition)
    vcf.header.add_tag_definition(hg19POS_definition)
    vcf.header.add_tag_definition(hg19CHROM_definition)

    # get chain file for liftover
    lo = LiftOver(args['chainfile'])

    # write header and then loop variants, adding liftover coordiantes to INFO fields when appropriate. write all variants.
    with open(args['outputfile'], 'w') as fo:
        vcf.write_header(fo)
        for vnt_obj in vcf.parse_variants():

            # generate hg19 LO coordinates based on CHROM and POS
            hits = lo.convert_coordinate(vnt_obj.CHROM, vnt_obj.POS-1)
            if len(hits) > 0:
                #add hg19_chr
                hg19CHROM_value = 'hg19_chr='+hits[0][0].split('chr')[1]
                vnt_obj.add_tag_info(hg19CHROM_value)
                #add hg19_pos
                hg19POS_value = 'hg19_pos='+str(hits[0][1]+1)
                vnt_obj.add_tag_info(hg19POS_value)

            # also want to incorporate END position for SV and CNV
            # check if "END" exists in INFO and if it does, try a liftover
            try:
                END = int(vnt_obj.INFO.split("END=")[1].split(";")[0])
            except:
                END = ''

            if END != '':
                hits_end = lo.convert_coordinate(vnt_obj.CHROM, END-1)
                if len(hits_end) > 0:
                    try:
                        #if hg19_chr is already defined, don't add it
                        vnt_obj.get_tag_value("hg19_chr")
                        #add hg19_end
                        hg19END_value = 'hg19_end='+str(hits_end[0][1]+1)
                        vnt_obj.add_tag_info(hg19END_value)
                    except:
                        #if hg19_chr is not defined, add hg19_chr
                        hg19CHROM_value = 'hg19_chr='+hits_end[0][0].split('chr')[1]
                        vnt_obj.add_tag_info(hg19CHROM_value)
                        #add hg19_end
                        hg19END_value = 'hg19_end='+str(hits_end[0][1]+1)
                        vnt_obj.add_tag_info(hg19END_value)
            vcf.write_variant(fo, vnt_obj)

    subprocess.run(["bgzip", args['outputfile']])
    subprocess.run(["tabix",args['outputfile']+".gz"])

################################################
#   Main
################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Add hg19 liftover coordinates to INFO field for each qualified variant')

    parser.add_argument('-i','--inputfile', help='input VCF file', required=True)
    parser.add_argument('-c','--chainfile', help='input hg38-to-hg19-chain file', required=True)
    parser.add_argument('-o','--outputfile', help='output VCF file', required=True)

    args = vars(parser.parse_args())

    main(args)
