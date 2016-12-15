#!/usr/bin/env python
"""
Script to run generate consensus sequences given a set of input reads

Do 'stonycss --help' for detailed information
"""

from   consensus import (scoring_function, do_consensus,
                         SCORING_FUNCTIONS, TRAVERSAL_ALGOS)
from   converter import  reverse_complement
from   poaligner import  align_sequences, convert_po_msa_to_dag, get_best_score

import argparse
import os
import pysam
import sys
import tempfile


PROG_DESC = """
This is a ccs (consensus calling) tool that, given a set of reads, generates a
consensus sequence.

It accepts as input a sorted-by-qname .bam file containing reads & a scoring
matrix file and it generates consensus sequences for each well among the reads.
The .bam file has to be in pacbio format. The tool allows configuring various
settings as seen below.
"""

# http://stackoverflow.com/questions/24101524/finding-median-of-list-in-python
def find_median(arr):
    if len(arr) == 1:
        return float(arr[0])
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
MEDIAN_DIFFER_ALLOWANCE = 5.0

DO_FILTERING      = True
ORDERING_ALGOS    = ["star_only_forward",
                     "star_forward_reverse",
                     "no_star_iterative",
                     "no_star_progressive",
                     "no_star_alternate_reversed_progressive"]
MY_ORDERING_ALGO  = "star_forward_reverse"
MY_SCORING_FUNC   = "edge_weight_based_score"
MY_TRAVERSAL_ALGO = "max_score"


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
    if len(seqs_well) == 0:
        return seqs_well
    seq_lengths = [len(seq.query) for seq in seqs_well]
    median = find_median(seq_lengths)
    out_seqs_well = []
    for seq, seq_length in zip(seqs_well, seq_lengths):
        if (median / MEDIAN_DIFFER_ALLOWANCE) < seq_length < (median * MEDIAN_DIFFER_ALLOWANCE):
            out_seqs_well.append(seq)
    return out_seqs_well


# ========================== Ordering Heuristics ===============================

def star_algorithm_ordering(sequences, score_matrix_file, only_forward=False):
    """
    Order the sequences using the STAR alogirthm
    
    If only_forward is True, all sequences are assumed to be in the forward
    direction. If it is set to False, the algorithm will choose the best among
    the forward and reverse orientations. The problem here is we don't know if
    the strands are in the forward/reverse direction. So, we perform a slight
    tweak on the original STAR algorithm.

    For each strand, we take the forward and reverse versions and find the alignment score
    with the forward and reverse versions of all other strands. We then compute the
    overall distance for each forward and reverse version of each strand with all
    other strands in both forward and reverse directions. Note that for a strand
    S_c, we will choose the maximum among the forward and reverse alignment scores
    with another strand S_i (Same will be done for all other strands)
    """
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
    read_order = {s.query: i for (i, s) in enumerate(seqs_well)}

    # Filters
    # A lot of redundant looping here. Optimize?
    if DO_FILTERING:
        seqs_well = do_adapter_filter(seqs_well)
        seqs_well = do_read_quality_filter(seqs_well)
        seqs_well = do_snr_filter(seqs_well)
        seqs_well = do_length_filter(seqs_well)
        seqs_well = do_median_filter(seqs_well)

    if not enough_required_sequences(seqs_well):
        return

    log_info('Adding %s sequences with id %s for ccs (filtered?=%s)' % (len(seqs_well), cur_seq_id, DO_FILTERING))
    # Good to add these sequences for this id to our seq_data for doing ccs
    sequences = [s.query for s in seqs_well]
    seq_data[cur_seq_id] = {'sequences': sequences}
    if MY_ORDERING_ALGO == 'no_star_alternate_reversed_progressive':
        for i, seq in enumerate(sequences):
            # Reverse-complement every odd-numbered sequence
            real_index = read_order[seq]
            if real_index % 2 != 0:
                seq_data[cur_seq_id]['sequences'][i] = reverse_complement(seq)


# =============================== Main CCS =====================================

def do_stonyccs(well_id, seqs, score_matrix_file):
    po_msa_f = tempfile.NamedTemporaryFile(delete=False)
    po_msa_f.close()

    if MY_ORDERING_ALGO == 'star_only_forward':
        ordered_seqs = star_algorithm_ordering(seqs, score_matrix_file, only_forward=True)
    elif MY_ORDERING_ALGO == 'star_forward_reverse':
        ordered_seqs = star_algorithm_ordering(seqs, score_matrix_file, only_forward=False)
    else:
        ordered_seqs = seqs

    do_progressive = False
    if MY_ORDERING_ALGO in ('no_star_progressive', 'no_star_alternate_reversed_progressive'):
        do_progressive = True

    log_info("Doing ccs for id %s" % well_id)
    align_sequences(ordered_seqs, score_matrix_file, po_msa_f.name, do_progressive=do_progressive)
    dag = convert_po_msa_to_dag(po_msa_f.name)
    os.unlink(po_msa_f.name)

    # convert to final CCS
    # assignment of scoring function
    scoring_function(dag, scoring_func=MY_SCORING_FUNC)
    # generating consensus based on given traversal algorithm
    ccs = do_consensus(dag, traversal_algo=MY_TRAVERSAL_ALGO)

    return ccs


# ==============================================================================

def parse_opts():
    global MIN_REQUIRED_SEQS, MIN_READ_QUALITY, MIN_SNR, MIN_READ_LENGTH, \
           MEDIAN_DIFFER_ALLOWANCE, MAX_READ_LENGTH, LOG_FH, MY_ORDERING_ALGO, \
           MY_SCORING_FUNC, MY_TRAVERSAL_ALGO, DO_FILTERING

    parser = argparse.ArgumentParser(description=PROG_DESC)

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
    parser.add_argument("--median_differ_allowance", type=int,
        help="Only allow reads whose length is greater than median/allowance and less than"
             "median*allowance (default %s)" % MEDIAN_DIFFER_ALLOWANCE)

    parser.add_argument("--disable_filters", action="store_true",
        help="If this flag is set, no filtering will be done (but wells with \n"
             "lesser than MIN_REQUIRED_SEQS sequences will still be rejected)")
    parser.add_argument("--ordering_algo", type=str,
        help=("Ordering algorithm to use. Specify one in \n" + 
              "[" + ", ".join(ORDERING_ALGOS) + "]\n"
              "(default %s)" % MY_ORDERING_ALGO))
    parser.add_argument("--scoring_func", type=str,
        help=("Graph scoring function to use. Specify one in \n" + 
              "[" + ", ".join(SCORING_FUNCTIONS) + "]\n"
              "(default %s)" % MY_SCORING_FUNC))
    parser.add_argument("--traversal_algo", type=str,
        help=("Graph traversal algorithm to use. Specify one in \n" + 
              "[" + ", ".join(TRAVERSAL_ALGOS) + "]\n"
              "(default %s)" % MY_TRAVERSAL_ALGO))

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
    if opts.median_differ_allowance:
        MEDIAN_DIFFER_ALLOWANCE = opts.median_differ_allowance

    DO_FILTERING = not(opts.disable_filters)
    if opts.ordering_algo:
        if opts.ordering_algo not in ORDERING_ALGOS:
            raise ValueError("Invalid ordering algo: " + opts.ordering_algo)
        MY_ORDERING_ALGO = opts.ordering_algo
    if opts.scoring_func:
        if opts.scoring_func not in SCORING_FUNCTIONS:
            raise ValueError("Invalid scoring function: " + opts.scoring_func)
        MY_SCORING_FUNC = opts.scoring_func
    if opts.traversal_algo:
        if opts.traversal_algo not in TRAVERSAL_ALGOS:
            raise ValueError("Invalid traversal algo: " + opts.traversal_algo)
        MY_TRAVERSAL_ALGO = opts.traversal_algo

    if not opts.log_file:
        opts.log_file = os.path.join(os.getcwd(), 'stonyccs_report.txt')
    LOG_FH = open(opts.log_file, 'w')

    message = "Command being run: " + \
              "python stonyccs.py {0} {1} ".format(opts.input_file, opts.output_file_prefix) + \
              "--matrix_file {0} ".format(opts.matrix_file) + \
              "--ordering_algo {0} ".format(MY_ORDERING_ALGO) + \
              "--scoring_func {0} ".format(MY_SCORING_FUNC) + \
              "--traversal_algo {0} ".format(MY_TRAVERSAL_ALGO)
    if not DO_FILTERING:
        message += "--disable_filters "
    log_info(message)

    return opts


def main():
    opts = parse_opts()

    print("All logs go to %s..." % LOG_FH.name)

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
        if seq_id > cur_seq_id:
            # Process previous sequences
            process_and_filter_seqs(cur_seq_id, seqs_well, seq_data)
            seqs_well = []
            cur_seq_id = seq_id
        elif seq_id < cur_seq_id:
            raise ValueError('This program expects a sorted .bam file')
        seqs_well.append(line)
    process_and_filter_seqs(cur_seq_id, seqs_well, seq_data)
    total_seqs_read = i + 1
    
    print("\nRead %s sequences\n" % total_seqs_read)
    print("Doing ccs for %s wells...\n" % len(seq_data))

    # Step 2: Do ccs for chosen wells
    log_info("Configuration for ccs - "
             "Ordering_algo: {0}, Scoring_function: {1}, Traversal_algo: {2}, "
             "Doing Filtering?: {3}".format(MY_ORDERING_ALGO, MY_SCORING_FUNC,
                                            MY_TRAVERSAL_ALGO, DO_FILTERING))
    ccs_seqs = {}
    seqs_used_for_ccs = 0
    for well_id in sorted(seq_data):
        seqs_used_for_ccs += len(seq_data[well_id]['sequences']) 
        ccs_seq = do_stonyccs(well_id, seq_data[well_id]['sequences'], opts.matrix_file)
        ccs_seqs[well_id] = ccs_seq

    # Step 3: Write output to fasta file 
    if seqs_used_for_ccs > 0:
        fastaf = open(opts.output_file_prefix + '.fa', 'w')
        for well_id in sorted(ccs_seqs):
            fastaf.write('>' + str(well_id) + '/stonyccs\n')
            fastaf.write(ccs_seqs[well_id] + '\n')
        fastaf.close()
        log_info("Generated consensus file - %s" % fastaf.name)

    log_info("Read total {0} reads from input file, did ccs on total {1} reads "
             "({2} separate wells). Used {3:.2f} % of the input reads".format(
                total_seqs_read, seqs_used_for_ccs, len(seq_data),
                (seqs_used_for_ccs / float(total_seqs_read))*100.0))

    # Close the log file
    LOG_FH.close()

if __name__ == '__main__':
    main()
