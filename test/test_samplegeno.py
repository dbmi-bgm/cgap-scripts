#################################################################
#   Libraries
#################################################################
import sys, os
import pytest
import filecmp
from granite.lib import vcf_parser


samplegeno = __import__("samplegeno")


def test_samplegeno(tmp_path):
    """
    Test for samplegeno
    """
    outputfile = f"{tmp_path}/output.vcf"
    # Variables and Run
    args = {"inputfile": "test/files/in_samplegeno.vcf", "outputfile": outputfile}

    samplegeno.main(args)
    result = filecmp.cmp(outputfile, "test/files/out_samplegeno.vcf")

    # Test
    assert result == True


def test_missing_format_samplegeno(tmp_path):
    """
    Test for samplegeno when AD FORMAT is missing
    """
    outputfile = f"{tmp_path}/output.vcf"
    # Variables and Run
    args = {"inputfile": "test/files/in_missingformat_samplegeno.vcf", "outputfile": outputfile}

    samplegeno.main(args)
    result = filecmp.cmp(outputfile, "test/files/out_missingformat_samplegeno.vcf")

    # Test
    assert result == True
