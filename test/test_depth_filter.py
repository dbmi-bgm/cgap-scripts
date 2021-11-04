#################################################################
#   Libraries
#################################################################
import sys, os
import pytest
from subprocess import Popen

from depth_filter import (
                            main as main_depth_filter
                           )

#################################################################
#   Tests
#################################################################

def test_depth_filter_proband_only():
    # Variables
    args = {'inputSampleVCF': 'test/files/depth_test_file_proband.vcf.gz', 'outputfile': 'proband3.vcf', 'min_depth': '3'}
    # Run
    main_depth_filter(args)
    # Tests

    a = os.popen('bgzip -c -d proband3.vcf.gz')
    b = os.popen('bgzip -c -d test/files/depth_test_file_proband_3.vcf.gz')

    assert [row for row in a.read()] == [row for row in b.read()]

    # Clean
    os.remove('proband3.vcf.gz')
    os.remove('proband3.vcf.gz.tbi')


def test_depth_filter_trio():
    # Variables
    args = {'inputSampleVCF': 'test/files/depth_test_file_trio.vcf.gz', 'outputfile': 'trio3.vcf', 'min_depth': '3'}
    # Run
    main_depth_filter(args)

    # Tests
    a = os.popen('bgzip -c -d trio3.vcf.gz')
    b = os.popen('bgzip -c -d test/files/depth_test_file_trio_3.vcf.gz')

    assert [row for row in a.read()] == [row for row in b.read()]

    # Clean
    os.remove('trio3.vcf.gz')
    os.remove('trio3.vcf.gz.tbi')
