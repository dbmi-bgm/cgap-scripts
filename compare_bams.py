#!/usr/bin/env python3

################################################
#
#   Compare reads line by line between two bam files
#       with the same number of reads
#
#   Michele Berselli
#   berselli.michele@gmail.com
#
################################################

################################################
#   Libraries
################################################
import sys, argparse, os

################################################
#   Functions
################################################
def to_generator(stream):
    for line in stream:
        yield line
    #end for
#end def

def main(args):

    # Variables
    threads = int(args['threads']) if args['threads'] else 2

    samfile1 = to_generator(os.popen('samtools view -@ {0} {1}'.format(threads//2, args['input_1'])))
    samfile2 = to_generator(os.popen('samtools view -@ {0} {1}'.format(threads//2, args['input_2'])))

    # Compare reads
    i, r1, r2 = 0, '', ''
    while True:
        try:
            r1 = next(samfile1)
            r2 = next(samfile2)
        except Exception:
            break
        #end try
        i += 1
        sys.stderr.write('\rread ' + str(i))
        if not r1 == r2:
            sys.stderr.write(' - mismatch\n')
            sys.stdout.write('@read {0}\n'.format(i))
            sys.stdout.write(r1)
            sys.stdout.write(r2)
        #end if
    #end while
    sys.stderr.write('\n')
#end def main

################################################
#   MAIN
################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('--input_1', help='input BAM file', required=True)
    parser.add_argument('--input_2', help='input BAM file', required=True)
    parser.add_argument('-t','--threads', help='number of threads to use for compression/decompression [1]', required=False)

    args = vars(parser.parse_args())

    main(args)
