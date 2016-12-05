#!/usr/bin/env python
"""
Module to convert sequence data between various formats
"""
import pysam


_SEQ_COMPLEMENTS = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def complement(sequence):
    """
    Returns a complemented version of a given DNA sequence string
    """
    return ''.join(map(lambda s: _SEQ_COMPLEMENTS[s], sequence))


def reverse_complement(sequence):
    """
    Returns a reversed complemented string of a given DNA sequence string
    """
    return complement(sequence)[::-1]


def sam_to_bam(sam_file, bam_file, check_sq=False):
    """
    Convert sam to bam file

    @sam_file: Input sam filename
    @bam_file: Output bam filename
    """
    in_f    = pysam.AlignmentFile(sam_file, 'r', check_sq=check_sq)
    in_segs = [seg for seg in in_f.fetch(until_eof=True)]

    out_f = pysam.AlignmentFile(bam_file, 'wb', header=in_f.header)
    for seg in in_segs:
        a = pysam.AlignedSegment()
        a.query_name           = seg.query_name
        a.query_sequence       = seg.query_sequence
        a.flag                 = seg.flag
        a.reference_id         = seg.reference_id
        a.reference_start      = seg.reference_start
        a.mapping_quality      = seg.mapping_quality
        a.cigar                = seg.cigar
        a.next_reference_id    = seg.next_reference_id
        a.next_reference_start = seg.next_reference_start
        a.template_length      = seg.template_length
        a.query_qualities      = seg.query_qualities
        a.tags                 = seg.tags
        outf.write(a)

    in_f.close()
    out_f.close()

