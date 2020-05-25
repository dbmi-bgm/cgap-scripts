#!/usr/bin/env python3

################################################
#
#   Remove duplicate lines from VCF file,
#    keep the longer and more informative one
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

################################################
#
#   Object
#
################################################
class Variant(object):
    ''' object to store information for variant in vcf format '''

    def __init__(self, line_strip):
        ''' initialize Variant object '''
        line_split = line_strip.split('\t')
        self.CHROM = line_split[0]
        self.POS = int(line_split[1])
        self.ID = line_split[2]
        self.REF = line_split[3]
        self.ALT = line_split[4]
        self.QUAL = line_split[5]
        self.FILTER = line_split[6]
        self.INFO = line_split[7]
        self.FORMAT = line_split[8]
        self.GENOTYPES = line_split[9:]
    #end def

    def to_string(self):
        ''' variant as string rapresentation '''
        variant_as_string = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t'.format(self.CHROM,
                                                                                    self.POS,
                                                                                    self.ID,
                                                                                    self.REF,
                                                                                    self.ALT,
                                                                                    self.QUAL,
                                                                                    self.FILTER,
                                                                                    self.INFO,
                                                                                    self.FORMAT)
        return variant_as_string + '\t'.join(self.GENOTYPES) + '\n'
    #end def

#end class Variant

################################################
#
#   Functions
#
################################################
def is_longer(vnt_1, vnt_2):
    ''' '''
    l_1, l_2 = len(vnt_1.INFO.split(';')), len(vnt_2.INFO.split(';'))
    if l_1 == l_2:
        if len(vnt_1.INFO) > len(vnt_2.INFO):
            return vnt_1
        else: return vnt_2
        #end if
    elif l_1 > l_2:
        return vnt_1
    else: return vnt_2
    #end if
#end def

################################################
#   runner
################################################
def main(args):
    ''' '''
    # Buffer output
    fo = open(args['outputfile'], 'w')

    # Remove duplicate lines
    with open(args['inputfile']) as fi:
        tmp_vnt, vnt = '', ''
        for line in fi:
            if line.startswith('#'):
                fo.write(line)
            else:
                vnt = Variant(line.rstrip())
                if not tmp_vnt:
                    tmp_vnt = vnt
                else:
                    if vnt.CHROM == tmp_vnt.CHROM \
                            and vnt.POS == tmp_vnt.POS \
                            and vnt.REF == tmp_vnt.REF \
                            and vnt.ALT == tmp_vnt.ALT:
                        tmp_vnt = is_longer(tmp_vnt, vnt)
                    else:
                        fo.write(tmp_vnt.to_string())
                        tmp_vnt = vnt
                    #end if
                #end if
            #end if
        #end for
        # write last variant
        fo.write(tmp_vnt.to_string())
        #end if
    #end with

    # Close buffer output
    fo.close()
#end def main

################################################
#
#   MAIN
#
################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Remove duplicate variants from sorted input VCF file. Keep the most informative variant')

    parser.add_argument('-i','--inputfile', help='input VCF file', required=True)
    parser.add_argument('-o','--outputfile', help='output VCF file', required=True)

    args = vars(parser.parse_args())

    main(args)

#end if
