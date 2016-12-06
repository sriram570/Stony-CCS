#!/usr/bin/env python

from   consensus import score_assignment, do_consensus
from   converter import reverse_complement
from   poaligner import align_sequences, convert_po_msa_to_dag

import argparse
import os
import pysam
import sys
import tempfile


# Constants
MIN_REQUIRED_SEQS = 3


def process_and_filter_seqs(cur_seq_id, seqs_well, seq_data):
    # Filter for number of sequences
    if len(seqs_well) < MIN_REQUIRED_SEQS:
        return

    # Good to add these sequences for this id to our seq_data for doing ccs
    seq_data[cur_seq_id] = {'sequences': []}
    for i, seq in enumerate(seqs_well):
        # Reverse-complement every odd-numbered sequence
        seq_string = seq.query if i % 2 == 0 else reverse_complement(seq.query)
        seq_data[cur_seq_id]['sequences'].append(seq_string)


def do_stonyccs(seqs, score_matrix_file):
    po_msa_f = tempfile.NamedTemporaryFile(delete=False)
    po_msa_f.close()

    align_sequences(seqs, score_matrix_file, po_msa_f.name)
    dag = convert_po_msa_to_dag(po_msa_f.name)
    os.unlink(po_msa_f.name)

    # Convert to final CCS
    score_assignment(dag)
    ccs = do_consensus(dag)

    return ccs


def parse_opts():
    parser = argparse.ArgumentParser()

    parser.add_argument("input_file", help="Input bam file")
    parser.add_argument("output_file_prefix", help="Prefix for output file")

    parser.add_argument("--matrix_file", help="Score matrix file for doing ccs",
                                         required=True)

    opts = parser.parse_args()
    return opts


def main():
    opts = parse_opts()
    
    inf = pysam.AlignmentFile(opts.input_file, 'rb', check_sq=False)
    
    seq_data = {}

    # Step 1: Read all sequences and choose the relevant ones to do ccs for
    cur_seq_id = -1
    seqs_well = []
    # Start reading
    for line in inf.fetch(until_eof=True):
        qname  = line.qname # Looks like "name/12345/23_34"
        seq_id = int(qname.split('/')[1])
        if seq_id != cur_seq_id and cur_seq_id != -1:
            # Process previous sequences
            process_and_filter_seqs(cur_seq_id, seqs_well, seq_data)
            seqs_well = []
            cur_seq_id = seq_id 
        seqs_well.append(line)
    process_and_filter_seqs(cur_seq_id, seqs_well, seq_data)

    # Step 2: Do ccs for chosen wells
    ccs_seqs = {}
    for well_id in seq_data:
        ccs_seq = do_stonyccs(seq_data[well_id]['sequences'], opts.matrix_file)
        ccs_seqs[well_id] = ccs_seq

    # Step 3: Write output to fasta file 
    fastaf = open(opts.output_file_prefix + '.fa', 'w')
    for well_id in sorted(ccs_seqs):
        fastaf.write('>' + str(well_id) + '/stonyccs\n')
        fastaf.write(ccs_seqs[well_id] + '\n')
    fastaf.close()

if __name__ == '__main__':
    main()
