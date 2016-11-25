#!/usr/local/bin/python
import os
import sys
import pysam

def main():
    try:
        infile = sys.argv[1]
    except IndexError:
        sys.exit("Run using bam_to_fasta.py <bam_file>")

    bamf = pysam.AlignmentFile(infile, "rb", check_sq=False)

    # dict for mapping file descriptors to corresponding FASTA files
    fds = {}

    # creating a parent directory (if it doesn't exist) containing all FASTA files
    DIR = "fasta_data/"
    if not os.path.exists(DIR):
        os.makedirs(DIR)

    # writing sequences to corresponding FASTA files
    # naming convention for files:
    # 1. file name -> 'sequence id'
    # 2. sequence names inside files -> {qStart}_{qEnd}
    for i, line in enumerate(bamf.fetch(until_eof=True)):
        seq_name = ">"
        #print line.query_name
        q_name = line.query_name.split('/')
        seq_name = seq_name + q_name[2]
        if q_name[1] not in fds:
            file_name = DIR + q_name[1] + ".fa"
            fds[q_name[1]] = open(file_name, 'a')
        #print (type(fds[q_name[1]]))
        fds[q_name[1]].write(seq_name + "\n")
        fds[q_name[1]].write(line.query_sequence + "\n")

    # closing all 'opened' FASTA files
    for fd in fds:
        fds[fd].close()

    # closing the BAM file
    bamf.close()

if __name__ == '__main__':
    main()