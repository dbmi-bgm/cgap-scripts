#################################################################
#   Libraries
#################################################################
import sys, os
import pytest
from subprocess import Popen

from liftover_hg19 import (
                            main as main_liftover_hg19
                           )


#################################################################
#   Tests
#################################################################

def test_snv():
    # Variables and Run
    args = {'inputfile': 'test/files/test_reduced_sorted','outputfile':'output.vcf','chainfile':'test/files/test.chain'}
    # Test
    main_liftover_hg19(args)
    a = os.popen('bgzip -c -d output.vcf.gz')
    b = os.popen('bgzip -c -d test/files/test_reduced_lo_sorted.vcf.gz')

    assert [row for row in a.read()] == [row for row in b.read()]
    # Clean
    os.remove('output.vcf.gz')
    os.remove('output.vcf.gz.tbi')

def test_sv():
    '''
    4 variants:
    1 that lifts over both positions,
    1 that lifts over start and fails end,
    1 that lifts over end and fails start,
    1 that fails both
    '''
    # Variables and Run
    args = {'inputfile': 'test/files/sv_liftover_in.vcf.gz','outputfile':'output.vcf','chainfile':'test/files/test.chain'}
    # Test
    main_liftover_hg19(args)
    a = os.popen('bgzip -c -d output.vcf.gz')
    b = os.popen('bgzip -c -d test/files/sv_liftover_out.vcf.gz')

    assert [row for row in a.read()] == [row for row in b.read()]
    # Clean
    os.remove('output.vcf.gz')
    os.remove('output.vcf.gz.tbi')
