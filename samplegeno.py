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
from granite.lib import vcf_parser

################################################
#
#   Functions
#
################################################
################################################
#   runner
################################################

"""
Replaces a format tag for the provided ID genotype with the given value

:param vnt_obj: Variant to be modified
:type vnt_obj: Variant

:param format_tag: Tag to be replaced
:type format_tag: str

:param id_genotype: ID of the sample that should be modified
:type id_genotype: str

:param value: new value
:type value: str

:param sep: separator used to separate fields
:type sep: str
"""


def replace_genotype_value(vnt_obj, format_tag, id_genotype, value, sep=":"):

    format_tag_idx = vnt_obj.FORMAT.split(sep).index(format_tag)
    format_fields = vnt_obj.GENOTYPES[id_genotype].split(sep)
    format_fields[format_tag_idx] = value
    vnt_obj.GENOTYPES[id_genotype] = sep.join(format_fields)


"""
Process the GT field for the sample

:param vnt_obj: Variant to be checked
:type vnt_obj: Variant

:param ID_genotype: ID of the sample that should be checked
:type ID_genotype: str

:return: GT values 
:rtype: list
"""


def process_GT_values(vnt_obj, ID_genotype):
    GT_0_ = None
    GT_1_ = None

    # Get possible alleles
    # REF is at idx 0, ALT idxs are maintained
    alleles = [vnt_obj.REF] + vnt_obj.ALT.split(",")

    GT = vnt_obj.get_genotype_value(ID_genotype, "GT").replace("|", "/").split("/")

    # for single GT values, we replace it with two values e.g 1 -> 1/1
    if len(GT) == 1:
        replace_genotype_value(vnt_obj, "GT", ID_genotype, f"{GT[0]}/{GT[0]}")
        GT_0_ = GT[0]
        GT_1_ = GT[0]
    elif len(GT) == 2:
        GT_0_ = GT[0]
        GT_1_ = GT[1]
    else:
        raise ValueError(
            f"GT for {ID_genotype} contains unexpected number of subfields, expected 1 or 2, got {len(GT)}. Variant: {vnt_obj.tostring()}."
        )

    if GT_0_ != ".":
        GT_0 = alleles[int(GT_0_)]
    else:
        GT_0 = GT_0_
    # end if
    if GT_1_ != ".":
        GT_1 = alleles[int(GT_1_)]
    else:
        GT_1 = GT_1_

    GT_alleles = f"{GT_0}/{GT_1}"
    GT = vnt_obj.get_genotype_value(ID_genotype, "GT")
    return GT, GT_alleles


"""
Process the AD field for the sample

:param vnt_obj: Variant to be checked
:type vnt_obj: Variant

:param ID_genotype: ID of the sample that should be checked
:type ID_genotype: str

:return: AD value
:rtype: str
"""


def process_AD_values(vnt_obj, ID_genotype):
    try:
        AD = vnt_obj.get_genotype_value(ID_genotype, "AD")
        # replace dots with Os to make it compatible with the portal
        if AD == ".":
            # two zeros for ref and alt alleles
            AD = "0,0"
            replace_genotype_value(vnt_obj, "AD", ID_genotype, AD)
        if "." in AD:
            # replace all dots with zeros 
            AD = AD.replace(".", "0")
            replace_genotype_value(vnt_obj, "AD", ID_genotype, AD)


    except ValueError:

        ad_num = 1  # number of ref + alt alleles, starts from 1 for the ref
        ad_num += len(vnt_obj.ALT.split(","))

        AD = ",".join(["0"] * ad_num)

        vnt_obj.add_values_genotype(ID_genotype, AD)

    return AD.replace(",", "/")


"""
Checks if the AD tag is in the FORMAT column, if it's not then adds the AD tag

:param vnt_obj: Variant to be checked
:type vnt_obj: Variant
"""


def check_AD_field(vnt_obj):
    try:
        AD_tag = vnt_obj.get_tag_value("AD")
    except ValueError:
        vnt_obj.add_tag_format("AD")


def main(args):
    """ """
    # Variables
    granite_def = "##GRANITE=<ID=SAMPLEGENO>"
    samplegeno_def = "##INFO=<ID=SAMPLEGENO,Number=.,Type=String,Description=\"Sample genotype information. Subembedded:'samplegeno':Format:'NUMGT|GT|AD|SAMPLEID'\">"

    # Buffers
    fo = open(args["outputfile"], "w")

    # Creating Vcf object
    vcf_obj = vcf_parser.Vcf(args["inputfile"])

    # Update and write header
    vcf_obj.header.add_tag_definition(granite_def + "\n" + samplegeno_def, "INFO")
    vcf_obj.write_header(fo)

    # Reading variants and adding samplegeno
    for vnt_obj in vcf_obj.parse_variants():
        # Init empty samplegeno
        samplegeno = []
        vnt_obj.complete_genotype()
        check_AD_field(vnt_obj)
        # Calculate samplegeno
        for ID_genotype in vnt_obj.IDs_genotypes:
            samplegeno_ = []
            GT, GT_alleles = process_GT_values(vnt_obj, ID_genotype)
            AD = process_AD_values(vnt_obj, ID_genotype)

            samplegeno_.append(GT)
            samplegeno_.append(GT_alleles)
            samplegeno_.append(AD)
            samplegeno_.append(ID_genotype)

            samplegeno.append("|".join(samplegeno_))
        # end for

        # Add samplegeno to variant INFO
        vnt_obj.add_tag_info("SAMPLEGENO={0}".format(",".join(samplegeno)))

        # Write variant
        vcf_obj.write_variant(fo, vnt_obj)
    # end for

    fo.close()


# end def main

################################################
#
#   MAIN
#
################################################
if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Add samplegeno information to variants in VCF file. SAMPLEGENO tag added to INFO field"
    )

    parser.add_argument("-i", "--inputfile", help="input VCF file", required=True)
    parser.add_argument("-o", "--outputfile", help="output VCF file", required=True)

    args = vars(parser.parse_args())

    main(args)

# end if
