#!/usr/bin/env python3

################################################
#
#   Split a bam file in consecutive chunks,
#     each chunk contains specified number of reads
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

################################################
#   Functions
################################################
def check_EOF(filename):
    EOF_hex = b'\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00'
    size = os.path.getsize(filename)
    fb = open(filename, "rb")
    fb.seek(size - 28)
    EOF = fb.read(28)
    fb.close()
    if EOF != EOF_hex:
        sys.stderr.write(' - EOF is missing, adding EOF\n')
        fb = open(filename, "ab")
        fb.write(EOF_hex)
        fb.close()
    else:
        sys.stderr.write(' - EOF is present, skipping\n')
    #end if
#end def

def write_bam(name, directory, chunk, header, reads, threads):
    # Opening output file
    filename = directory + '/' + '{0}.{1}.chunk.bam'.format(name, chunk)
    bamfile = open(filename, 'w')

    # Buffer to stream output
    pipe_out = subprocess.Popen(['samtools', 'view', '-b', '-S', '-h', '-@ {0}'.format(threads//2), '-'], stdin=subprocess.PIPE, stdout=bamfile)

    # Writing header
    pipe_out.stdin.write(header.encode())

    # Writing reads
    for read in reads:
        pipe_out.stdin.write(read.encode())
    #end for

    # Close stream
    pipe_out.communicate()
    bamfile.close()

    # Check if bamfile has EOF, if not add EOF
    check_EOF(filename)
#end def

def main(args):

    # Variables
    directory = args['directory'] if args['directory'] else '.'
    threads = int(args['threads']) if args['threads'] else 2
    nreads = int(args['nreads']) if args['nreads'] else 5000000
    name = args['inputfile'].split('/')[-1].split('.')[0]

    # Open bamfile header
    samfile = os.popen('samtools view -H -@ {0} {1}'.format(threads, args['inputfile']))

    # Reading and saving header
    header = ''
    for line in samfile:
        header += line
    #end for

    # Open bamfile reads
    samfile = os.popen('samtools view -@ {0} {1}'.format(threads//2, args['inputfile']))

    # Reading reads
    i, c, reads = 0, 0, []
    for read in samfile:
        i += 1
        reads.append(read)
        if i == nreads:
            c += 1
            sys.stderr.write('write chunk n ' + str(c))
            write_bam(name, directory, c, header, reads, threads)
            i, reads = 0, []
        #end if
    #end for
    if i:
        c += 1
        sys.stderr.write('write chunk n ' + str(c))
        write_bam(name, directory, c, header, reads, threads)
    #end if
    sys.stderr.write('finished\n')

#end def main

################################################
#   MAIN
################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-i','--inputfile', help='input BAM file', required=True)
    parser.add_argument('-d','--directory', help='directory to use to write results [.]', required=False)
    parser.add_argument('-t','--threads', help='number of threads to use for compression/decompression [1]', required=False)
    parser.add_argument('-n','--nreads', help='number of reads to write per chunk [5000000]', required=False)

    args = vars(parser.parse_args())

    main(args)
