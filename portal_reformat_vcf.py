#!/usr/bin/env python3

################################################
#
#   Reformat VCF and add custom tags and fields
#       in VCF file for portal ingestion
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
# shared_functions as *
from granite.lib.shared_functions import *
# shared_vars
from granite.lib.shared_vars import VEP_encode
from granite.lib.shared_vars import DStags

################################################
#
#   Functions
#
################################################
def get_maxds(vnt_obj, SpAItag_list, SpAI_idx_list):
    ''' '''
    # if SpliceAI is within VEP
    # fetching only the first transcript
    # expected the same scores for all transcripts
    SpAI_vals = []
    for i, SpAItag in enumerate(SpAItag_list):
        SpAI_val = get_tag_idx(vnt_obj, SpAItag, SpAI_idx_list[i])
        # if SpliceAI is with VEP and is at the end of Format
        # need to remove , that separate next transcript
        try: SpAI_vals.append(float(SpAI_val.split(',')[0]))
        except Exception:
            return None
        #end try
    #end for
    if SpAI_vals:
        return max(SpAI_vals)
    #end if
    return None
#end def

def get_worst_trscrpt(vnt_obj, VEPtag, IMPCT_encode, CNONICL_idx, IMPCT_idx):
    ''' '''
    # Get worst impact among transcripts
    # VEP_field from granite shared_functions,
    #   return a list of value at index across all transcripts
    IMPCT_list = VEP_field(vnt_obj, IMPCT_idx, VEPtag)
    IMPCT_encoded = [IMPCT_encode[IMPCT] for IMPCT in IMPCT_list]
    try:
        worst_IMPCT = min(IMPCT_encoded)
    except Exception: return None
    # Get VEP
    val_get = vnt_obj.get_tag_value(VEPtag)
    # Get worst transcripts
    worst_trscrpt_list = []
    trscrpt_list = val_get.split(',')
    # Check transcripts
    for trscrpt in trscrpt_list:
        trscrpt_impct = trscrpt.split('|')[IMPCT_idx]
        if IMPCT_encode[trscrpt_impct] == worst_IMPCT:
            # Check canonical
            trscrpt_cnonicl = trscrpt.split('|')[CNONICL_idx]
            if trscrpt_cnonicl == 'YES' or \
                trscrpt_cnonicl == '1':
                return trscrpt
            #end if
            worst_trscrpt_list.append(trscrpt)
        #end if
    #end for
    return worst_trscrpt_list[0]
#end def

def update_worst(vnt_obj, VEPtag, worst_trscrpt):
    ''' '''
    # Get VEP
    val_get = vnt_obj.get_tag_value(VEPtag)
    trscrpt_update = []
    trscrpt_list = val_get.split(',')
    # Update transcripts
    for trscrpt in trscrpt_list:
        if trscrpt == worst_trscrpt:
            trscrpt += '|1'
        else: trscrpt += '|0'
        #end if
        trscrpt_update.append(trscrpt)
    #end for
    return ','.join(trscrpt_update)
#end def

################################################
#   runner
################################################
def main(args):
    ''' '''
    # Variables
    is_verbose = args['verbose']
    VEPtag = 'CSQ'
    IMPCT_encode = {'HIGH': 1, 'MODERATE': 2, 'LOW': 3, 'MODIFIER': 4}

    # Definitions
    vep_init = '##VEP=<ID={0}>'.format(VEPtag)
    genes_init = '##CGAP=<ID=GENES>'
    spliceai_def = '##INFO=<ID=spliceaiMaxds,Number=1,Type=Float,Description="SpliceAI max delta score">'
    genes_def = '##INFO=<ID=GENES,Number=.,Type=String,Description=". Subembedded:\'genes\':Format:\'most_severe_gene|most_severe_transcript|most_severe_feature_ncbi|most_severe_hgvsc|most_severe_hgvsp|most_severe_amino_acids|most_severe_sift_score|most_severe_polyphen_score|most_severe_maxentscan_diff|most_severe_consequence\'">'
    variant_def = '##INFO=<ID=variantClass,Number=1,Type=String,Description="Variant type">'

    # Buffers
    fo = open(args['outputfile'], 'w')

    # Creating Vcf object
    vcf_obj = vcf_parser.Vcf(args['inputfile'])

    # Modify VEP definition
    vep_def = '##INFO=<ID={0},Number=.,Type=String,Description="Consequence annotations from Ensembl VEP.  Subembedded:\'transcript\':Format:\'{1}\'">'
    vep_field_list = []
    for line in vcf_obj.header.definitions.split('\n')[:-1]:
        if line.startswith('##INFO=<ID=' + VEPtag + ','): ##<tag_type>=<ID=<tag>,...
            format = line.split('Format:')[1]
            # Cleaning format
            format = format.replace(' ', '')
            format = format.replace('\'', '')
            format = format.replace('\"', '')
            format = format.replace('>', '')
            # Update definition
            vep_field_list = format.split('|')
            vep_field_list.append('most_severe')
            vep_def = vep_def.format(VEPtag, '|'.join(vep_field_list))
            break
        #end if
    #end for

    # Remove older VEP definition
    vcf_obj.header.remove_tag_definition(VEPtag)

    # Update and write custom definitions
    vcf_obj.header.add_tag_definition(vep_init + '\n' + genes_init, 'INFO')
    vcf_obj.header.add_tag_definition(spliceai_def, 'INFO')
    vcf_obj.header.add_tag_definition(genes_def, 'INFO')
    vcf_obj.header.add_tag_definition(variant_def, 'INFO')
    vcf_obj.header.add_tag_definition(vep_def, 'INFO')

    # Write header
    vcf_obj.write_header(fo)

    # Get SpliceAI ds indexes
    # DStags import from granite.shared_vars
    SpAItag_list, SpAI_idx_list = [], []
    for DStag in DStags:
        tag, idx = vcf_obj.header.check_tag_definition(DStag)
        SpAItag_list.append(tag)
        SpAI_idx_list.append(idx)
    #end for

    # Get VEP indexes
    IMPCT_idx = vcf_obj.header.get_tag_field_idx(VEPtag, 'IMPACT')
    CNONICL_idx = vcf_obj.header.get_tag_field_idx(VEPtag, 'CANONICAL')

    ENSG_idx = vcf_obj.header.get_tag_field_idx(VEPtag, 'Gene')
    ENST_idx = vcf_obj.header.get_tag_field_idx(VEPtag, 'Feature')
    MANE_idx = vcf_obj.header.get_tag_field_idx(VEPtag, 'MANE') #feature_ncbi
    HGVSC_idx = vcf_obj.header.get_tag_field_idx(VEPtag, 'HGVSc')
    HGVSP_idx = vcf_obj.header.get_tag_field_idx(VEPtag, 'HGVSp')
    AACIDS_idx = vcf_obj.header.get_tag_field_idx(VEPtag, 'Amino_acids')
    SIFT_idx = vcf_obj.header.get_tag_field_idx(VEPtag, 'SIFT')
    PPHEN_idx = vcf_obj.header.get_tag_field_idx(VEPtag, 'PolyPhen')
    MAXENTDIFF_idx = vcf_obj.header.get_tag_field_idx(VEPtag, 'MaxEntScan_diff')
    CONSEQUENCE_idx = vcf_obj.header.get_tag_field_idx(VEPtag, 'Consequence')

    # Reading variants and adding new tags
    for i, vnt_obj in enumerate(vcf_obj.parse_variants()):
        if is_verbose:
            sys.stderr.write('\r' + str(i+1))
            sys.stderr.flush()
        #end if
        # Get max SpliceAI max_ds
        maxds = get_maxds(vnt_obj, SpAItag_list, SpAI_idx_list)

        # Get most severe transcript
        worst_trscrpt = get_worst_trscrpt(vnt_obj, VEPtag, IMPCT_encode, CNONICL_idx, IMPCT_idx)

        if not worst_trscrpt: continue
        #end if

        # Get variant class
        # import from granite.shared_functions
        clss = variant_type_ext(vnt_obj.REF, vnt_obj.ALT)

        # Add MAXDS to variant INFO
        if maxds: vnt_obj.add_tag_info('spliceaiMaxds={0}'.format(maxds))
        #end if
        # Add CLASS to variant INFO
        vnt_obj.add_tag_info('variantClass={0}'.format(clss.upper()))

        # Update and replace VEP tag in variant INFO
        # Adding field most_severe (0|1) to transcripts
        VEP_update = update_worst(vnt_obj, VEPtag, worst_trscrpt)
        # Replace VEP
        vnt_obj.remove_tag_info(VEPtag)
        vnt_obj.add_tag_info('{0}={1}'.format(VEPtag, VEP_update))

        # Add GENES to variant INFO
        worst_trscrpt_ = worst_trscrpt.split('|')
        genes = '{0}|{1}|{2}|{3}|{4}|{5}|{6}|{7}|{8}|{9}'.format(
                                                        worst_trscrpt_[ENSG_idx],
                                                        worst_trscrpt_[ENST_idx],
                                                        worst_trscrpt_[MANE_idx],
                                                        worst_trscrpt_[HGVSC_idx],
                                                        worst_trscrpt_[HGVSP_idx],
                                                        worst_trscrpt_[AACIDS_idx],
                                                        worst_trscrpt_[SIFT_idx],
                                                        worst_trscrpt_[PPHEN_idx],
                                                        worst_trscrpt_[MAXENTDIFF_idx],
                                                        worst_trscrpt_[CONSEQUENCE_idx]
                                                )
        vnt_obj.add_tag_info('GENES={0}'.format(genes))

        # Write variant
        vcf_obj.write_variant(fo, vnt_obj)
    #end for

    # Close buffers
    sys.stderr.write('\n')
    fo.close()
#end def main

################################################
#
#   MAIN
#
################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Reformat VCF and add custom tags and fields in VCF file for portal ingestion')

    parser.add_argument('-i','--inputfile', help='input VCF file', required=True)
    parser.add_argument('-o','--outputfile', help='output VCF file', required=True)
    parser.add_argument('--verbose', help='verbose', action='store_true', required=False)

    args = vars(parser.parse_args())

    main(args)

#end if
