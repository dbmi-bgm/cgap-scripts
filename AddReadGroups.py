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
        if type(header[key]) is list:
            for d in header[key]:
                s = '@' + str(key) + '\t'
                for k, v in d.items():
                    s += str(k) + ':' + str(v) + '\t'
                #end for
                str_header += s.rstrip() + '\n'
            #end for
        elif type(header[key]) is dict:
            s = '@' + str(key) + '\t'
            for k, v in header[key].items():
                s += str(k) + ':' + str(v) + '\t'
            #end for
            str_header += s.rstrip() + '\n'
        #end if
    #end for
    return str_header
#end def

def get_read_group(query_name):
    qn_as_list = query_name.split(':')
    l = len(qn_as_list) # number of fields in query_name
    if l == 7: # new style header
        return '_'.join(qn_as_list[:4])
    elif l == 5: # old style header
        return '_'.join(qn_as_list[:2])
    elif l == 1: # weird style header, just use a placeholder
        return 'RG_PLACEHOLDER'
    else:
        sys.exit('\nFORMAT ERROR: read format {0} not recognized\n'
                    .format(query_name))
    #end if
#end def

def check_EOF(filename):
    EOF_hex = b'\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00'
    size = os.path.getsize(filename)
    fb = open(filename, "rb")
    fb.seek(size - 28)
    EOF = fb.read(28)
    fb.close()
    if EOF != EOF_hex:
        sys.stderr.write('EOF is missing, adding EOF\n')
        fb = open(filename, "ab")
        fb.write(EOF_hex)
        fb.close()
    else:
        sys.stderr.write('EOF is present, skipping\n')
    #end if
#end def

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
        ID = get_read_group(read.query_name)
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
        ID = get_read_group(read.query_name)
        # not using pysam to add read tag because somehow it is doing
        # something extremely slow and locking all the multi-threading
        read_str = read.tostring() + '\t' + 'RG:Z:{0}'.format(ID) + '\n'
        pipe_out.stdin.write(read_str.encode())
    #end for

    pipe_out.communicate()
    bamfile.close()

    # Check if bamfile has EOF, if not add EOF
    check_EOF(filename)
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
