#!/usr/bin/env python3

################################################
#
#   Add read groups to a BAM file using pysam
#    and samtools to multi-thread compression
#
#   Michele Berselli
#   berselli.michele@gmail.com
#
################################################

################################################
#   Libraries
################################################
import sys, argparse, os
import subprocess
import copy as cp
import pysam as ps

################################################
#   Functions
################################################
def header_as_str(header):
    str_header = ''
    for key in header:
        for d in header[key]:
            s = '@' + str(key) + '\t'
            for k, v in d.items():
                s += str(k) + ':' + str(v) + '\t'
            #end for
            str_header += s.rstrip() + '\n'
        #end for
    #end for
    return str_header
#end def write_header

def main(args):

    # Variables
    directory = args['directory'] if args['directory'] else '.'
    platform = args['platform'] if args['platform'] else 'illumina'
    samplename = args['samplename']
    threads = int(args['threads']) if args['threads'] else 1

    # Open bamfile using samtools and send to pipe
    pipe_in = subprocess.Popen(['samtools', 'view', '-h', '-@ {0}'.format(threads), args['inputfile']], stdout=subprocess.PIPE)

    # Data structures
    IDs = set()

    # Reading header and all possible IDs
    samfile = ps.AlignmentFile(pipe_in.stdout, 'r')

    header = cp.deepcopy(dict(samfile.header))
    header.setdefault('RG', [])

    for read in samfile:
        ID = '_'.join(read.query_name.split(':')[:4])
        IDs.add(ID)
    #end for

    # Updating header with read groups
    {header['RG'].append({'ID': ID, 'PL': platform, 'PU': ID, 'LB': ID, 'SM': samplename}) for ID in IDs}

    # Opening output file
    filename = directory + '/' + args['inputfile'].split('/')[-1].split('.')[0] + '_rg' + '.bam'
    bamfile = open(filename, 'w')

    # Buffer to stream output
    pipe_out = subprocess.Popen(['samtools', 'view', '-b', '-S', '-h', '-@ {0}'.format(threads), '-'], stdin=subprocess.PIPE, stdout=bamfile)

    # Writing header
    pipe_out.stdin.write(header_as_str(header).encode())

    # Adding read groups to alignments
    pipe_in = subprocess.Popen(['samtools', 'view', '-h', '-@ {0}'.format(threads), args['inputfile']], stdout=subprocess.PIPE)
    samfile = ps.AlignmentFile(pipe_in.stdout, 'r')
    for read in samfile:
        ID = '_'.join(read.query_name.split(':')[:4])
        # not using pysam to add read tag because somehow it is doing
        # something extremely slow and locking all the multi-threading
        read_str = read.tostring() + '\t' + 'RG:Z:{0}'.format(ID) + '\n'
        pipe_out.stdin.write(read_str.encode())
    #end for

    bamfile.close()
#end def main

################################################
#   MAIN
################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Add read groups to a BAM file')

    parser.add_argument('-i','--inputfile', help='input aligned paired-ends BAM file', required=True)
    parser.add_argument('-d','--directory', help='directory to use to write results', required=False)
    parser.add_argument('-t','--threads', help='number of threads to use for compression/decompression', required=False)
    parser.add_argument('-s','--samplename', help='name of the sample', required=True)
    parser.add_argument('-p','--platform', help='name of the sequencing platform', required=False)

    args = vars(parser.parse_args())

    main(args)

#end if
