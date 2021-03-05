#################################################################
#   Libraries
#################################################################
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from add_altcounts_by_gene import (
                                    main as main_add_altcounts
                                  )

#################################################################
#   Tests
#################################################################
def test_run_rckTar():
    ''' '''
    # Variables
    args = {'inputfile': 'test/files/input_add_altcounts.vcf',
            'outputfile': 'test/files/test.out'}
    # Run
    main_add_altcounts(args)
    # Tests
    assert [row for row in open('test/files/input_add_altcounts.out')] == [row for row in open('test/files/test.out')]
    # Clean
    os.remove('test/files/test.out')
#end def
