#!/usr/bin/env python

from   consensus import score_assignment, do_consensus
from   converter import reverse_complement
from   poaligner import align_sequences, convert_po_msa_to_dag, get_best_score

import argparse
import os
import pysam
import sys
import tempfile


# http://stackoverflow.com/questions/24101524/finding-median-of-list-in-python
def find_median(arr):
    arr = sorted(arr)
    left_mid = (len(arr) - 1) / 2
    if len(arr) % 2 == 0:
        return float(arr[left_mid])
    else:
        return float((arr[left_mid] + arr[left_mid+1]) / 2.0) 


# Various tweaks and knobs
# These exist as global variables as we don't want to pass them around needlessly
MIN_REQUIRED_SEQS = 3
MIN_READ_QUALITY  = 0.8
MIN_SNR           = 3.75
MIN_READ_LENGTH   = 10
MAX_READ_LENGTH   = 7000

DISABLED_FEATURES = []

# Logger
LOG_FH = None
def log_info(message):
    LOG_FH.write('INFO: ' + message + '\n')

# ============================ Rejection checks ================================

def enough_required_sequences(seqs):
    return bool(len(seqs) >= MIN_REQUIRED_SEQS) 

# ================================= Filters ====================================

def do_adapter_filter(seqs_well):
    out_seqs_well = []
    for seq in seqs_well:
        for tag in seq.tags:
            if tag[0] == 'cx':
                break
        # 3 means both adapter_start and adapter_end are present
        # http://pacbiofileformats.readthedocs.io/en/3.0/BAM.html
        if int(tag[1]) == 3:
            out_seqs_well.append(seq)
    return out_seqs_well

def do_read_quality_filter(seqs_well):
    out_seqs_well = []
    for seq in seqs_well:
        for tag in seq.tags:
            if tag[0] == 'rq':
                break
        if float(tag[1]) >= MIN_READ_QUALITY:
            out_seqs_well.append(seq)
    return out_seqs_well

def do_snr_filter(seqs_well):
    out_seqs_well = []
    for seq in seqs_well:
        for tag in seq.tags:
            if tag[0] == 'sn':
                break
        if all(float(v) >= MIN_SNR for v in tag[1]):
            out_seqs_well.append(seq)
    return out_seqs_well

def do_length_filter(seqs_well):
    return [s for s in seqs_well if MIN_READ_LENGTH <= len(s.query) <= MAX_READ_LENGTH]

def do_median_filter(seqs_well):
    seq_lengths = [len(seq.query) for seq in seqs_well]
    median = find_median(seq_lengths)
    out_seqs_well = []
    for seq, seq_length in zip(seqs_well, seq_lengths):
        if (seq_length / 2.0) < median < (seq_length * 2.0):
            out_seqs_well.append(seq)
    return out_seqs_well


# ========================== Ordering Heuristics ===============================

def star_algorithm_ordering(sequences, score_matrix_file, only_forward=False):
    """
    Order the sequences using the STAR alogirthm
    
    If only_forward is 
    The problem here is we don't know if the strands are in the forward/reverse
    direction. So, we perform a slight tweak on the original STAR algorithm.

    For each strand, we take the forward and reverse versions and find the alignment score
    with the forward and reverse versions of all other strands. We then compute the
    overall distance for each forward and reverse version of each strand with all
    other strands in both forward and reverse directions. Note that for a strand
    S_c, we will choose the maximum among the forward and reverse alignment scores
    with another strand S_i (Same will be done for all other strands)
    """
    print("Doing STAR algorithm ordering for %s seqs..." % len(sequences))
    seq_data = {i: {'fw': seq,
                    'rv': reverse_complement(seq)} for i, seq in enumerate(sequences)}

    neg_inf  = -float("inf")

    scores = {}
    fw, rv = 'fw', 'rv'
    # Assign scores
    for i in range(len(seq_data)):
        for j in range(i+1, len(seq_data)):
            i_fw_seq, j_fw_seq = seq_data[i]['fw'], seq_data[j]['fw']
            scores.setdefault((i,j), {})
            scores.setdefault((j,i), {})
            # Only forward scores first
            fw_fw_score = get_best_score([i_fw_seq, j_fw_seq], score_matrix_file)
            scores[(i,j)][(fw,fw)] = fw_fw_score
            scores[(j,i)][(fw,fw)] = fw_fw_score
            if not only_forward:
                i_rv_seq, j_rv_seq = seq_data[i]['rv'], seq_data[j]['rv']
                # Reverse-reverse scores
                rv_rv_score = get_best_score([i_rv_seq, j_rv_seq], score_matrix_file)
                scores[(i,j)][(rv,rv)] = rv_rv_score
                scores[(j,i)][(rv,rv)] = rv_rv_score
                # Cross pairing
                i_fw_j_rv_score = get_best_score([i_fw_seq, j_rv_seq], score_matrix_file)
                i_rv_j_fw_score = get_best_score([i_rv_seq, j_fw_seq], score_matrix_file)
                scores[(i,j)][(fw,rv)] = i_fw_j_rv_score
                scores[(i,j)][(rv,fw)] = i_rv_j_fw_score
                scores[(j,i)][(fw,rv)] = i_rv_j_fw_score
                scores[(j,i)][(rv,fw)] = i_fw_j_rv_score

    # Choose best
    best_scores = {}
    for i in range(len(seq_data)):
        for orntn in ('fw', 'rv'):
            if only_forward and orntn == 'rv': continue
            score = 0
            for j in range(len(seq_data)):
                if i == j: continue
                score += max(scores[(i,j)][(orntn,'fw')], scores[(i,j)].get((orntn,'rv'),neg_inf))
            if score not in best_scores:
                best_scores[score] = (i, orntn)

    best_i, orntn = best_scores[max(best_scores)]
    ordered_scores = []
    for j in range(len(seq_data)):
        if j == best_i: continue
        best_orntn = 'fw' if scores[(best_i,j)][(orntn,'fw')] >= scores[(best_i,j)].get((orntn,'rv'),neg_inf) else 'rv'
        best_score = scores[(best_i,j)][(orntn,best_orntn)]
        ordered_scores.append([j, best_orntn, best_score])
    ordered_scores = sorted(ordered_scores, key=lambda s: s[2])

    # Order and return
    final_ordered_sequences = [seq_data[best_i][orntn]]
    for el in ordered_scores:
        final_ordered_sequences.append(seq_data[el[0]][el[1]])

    return final_ordered_sequences


# ============================= Main Read Sanitizer ============================

def process_and_filter_seqs(cur_seq_id, seqs_well, seq_data):
    # Rejection checks - 
    # (We check this later on as well but we have a lot of single-read wells
    #  and we want to reject them early for efficiency)
    if not enough_required_sequences(seqs_well):
        return

    # Save read order for later access if necessary
    read_order = {s.qname: i for (i, s) in enumerate(seqs_well)}

    # Filters
    # A lot of redundant looping here. Optimize?
    seqs_well = do_adapter_filter(seqs_well)
    seqs_well = do_read_quality_filter(seqs_well)
    seqs_well = do_snr_filter(seqs_well)
    seqs_well = do_length_filter(seqs_well)
    seqs_well = do_median_filter(seqs_well)

    if not enough_required_sequences(seqs_well):
        return

    log_info('Adding %s sequences with id %s for ccs' % (len(seqs_well), cur_seq_id))
    # Good to add these sequences for this id to our seq_data for doing ccs
    sequences = [s.query for s in seqs_well]
    seq_data[cur_seq_id] = {'sequences': sequences}
    #for i, seq in enumerate(seqs_well):
        # Reverse-complement every odd-numbered sequence
        # seq_string = seq.query if i % 2 == 0 else reverse_complement(seq.query)
        #seq_data[cur_seq_id]['sequences'].append(seq_string)



# =============================== Main CCS =====================================

def do_stonyccs(seqs, score_matrix_file):
    po_msa_f = tempfile.NamedTemporaryFile(delete=False)
    po_msa_f.close()

    ordered_seqs = star_algorithm_ordering(seqs, score_matrix_file)

    align_sequences(ordered_seqs, score_matrix_file, po_msa_f.name, do_progressive=False)
    dag = convert_po_msa_to_dag(po_msa_f.name)
    os.unlink(po_msa_f.name)

    # Convert to final CCS
    score_assignment(dag)
    ccs = do_consensus(dag)

    return ccs


# ==============================================================================

def parse_opts():
    global MIN_REQUIRED_SEQS, MIN_READ_QUALITY, MIN_SNR, MIN_READ_LENGTH, \
           MAX_READ_LENGTH, LOG_FH

    parser = argparse.ArgumentParser()

    parser.add_argument("input_file", help="Input bam file")
    parser.add_argument("output_file_prefix", help="Prefix for output file")

    parser.add_argument("--matrix_file", help="Score matrix file for doing ccs",
                                         required=True)

    parser.add_argument("--min_required_sequences", type=int,
        help="Min. no. of sequences required to do ccs (default %s)" % MIN_REQUIRED_SEQS)
    parser.add_argument("--min_read_quality", type=float,
        help="Min. read quality required to do ccs (default %s)" % MIN_READ_QUALITY)
    parser.add_argument("--min_snr", type=float,
        help="Min. snr on all 4 channels required to do ccs (default %s)" % MIN_SNR)
    parser.add_argument("--min_read_length", type=int,
        help="Min. read length required to do ccs (default %s)" % MIN_READ_LENGTH)
    parser.add_argument("--max_read_length", type=int,
        help="Max. read length limit for ccs (default %s)" % MAX_READ_LENGTH)

    parser.add_argument("--disable_features", help="Comma-separated feature names to disable")

    parser.add_argument("--log_file", type=str,
        help="Log file to write logs to. Defaults to stonyccs_report.txt in cwd")

    opts = parser.parse_args()

    if opts.min_required_sequences:
        MIN_REQUIRED_SEQS = opts.min_required_sequences
    if opts.min_read_quality:
        MIN_READ_QUALITY = opts.min_read_quality
    if opts.min_snr:
        MIN_SNR = opts.min_snr
    if opts.min_read_length:
        MIN_READ_LENGTH = opts.min_read_length
    if opts.max_read_length:
        MAX_READ_LENGTH = opts.max_read_length

    if not opts.log_file:
        opts.log_file = os.path.join(os.getcwd(), 'stonyccs_report.txt')
    LOG_FH = open(opts.log_file, 'w')

    return opts


def main():
    opts = parse_opts()

    inf = pysam.AlignmentFile(opts.input_file, 'rb', check_sq=False)
    
    seq_data = {}

    # Step 1: Read all sequences and choose the relevant ones to do ccs for
    cur_seq_id = -1
    seqs_well = []
    # Start reading
    for i, line in enumerate(inf.fetch(until_eof=True)):
        qname  = line.qname # Looks like "name/12345/23_34"
        seq_id = int(qname.split('/')[1])
        if i == 0:
            cur_seq_id = seq_id
        if seq_id != cur_seq_id:
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

    # Close the log file
    LOG_FH.close()

if __name__ == '__main__':
    main()
