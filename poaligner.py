#!/usr/bin/env python
"""
Module to perform Partial Order Alignments
"""
from __future__ import print_function

import os
import re
import subprocess
import tempfile

# See if poa is accessible when module loads
with open(os.devnull, 'w') as dnull:
    exit_code = subprocess.call('poa', shell=True, stdout=dnull, stderr=dnull)
    if exit_code == 127:
        raise EnvironmentError("'poa' not found in PATH. Add it and retry")


def _align(input_files_command,
           score_matrix_file,
           po_out_file,
           pir_out_file=None,
           clustal_out_file=None,
           do_global=False,
           do_progressive=True):
    """
    The common aligner used by other methods. Constructs a command and call's
    the C POA program. 
    """
    poa_command  = 'poa'
    poa_command += ' ' + input_files_command
    poa_command += ' ' + score_matrix_file

    poa_command += ' -po ' + po_out_file
    if pir_out_file:
        poa_command += ' -pir ' + pir_out_file
    if clustal_out_file:
        poa_command += ' -clustal ' + clustal_out_file

    if do_global:
        poa_command += ' -do_global'
    if do_progressive:
        poa_command += ' -do_progressive'

    print('Triggering Partial Order Alignment...')
    print('Running "%s"' % poa_command)

    subprocess.check_call(poa_command, shell=True)


def align_sequences_from_fasta_file(fasta_file,
                    score_matrix_file,
                    po_out_file,
                    pir_out_file=None,
                    clustal_out_file=None,
                    do_global=False,
                    do_progressive=True):
    """
    Align a set of sequences present in a fasta file into a po_msa
    """
    input_files_command = '-read_fasta ' + fasta_file
    _align(input_files_command, score_matrix_file, po_out_file, 
            pir_out_file, clustal_out_file, do_global, do_progressive)


def align_sequences(sequences,
                    score_matrix_file,
                    po_out_file,
                    pir_out_file=None,
                    clustal_out_file=None,
                    do_global=False,
                    do_progressive=True):
    """
    Align a list of sequence strings into a po_msa
    """
    temp_fasta_file = tempfile.NamedTemporaryFile(delete=False)

    # Converting sequences to fasta format
    for i, sequence in enumerate(sequences):
        temp_fasta_file.write('>' + 'Sequence_%s' % i + '\n')
        temp_fasta_file.write(sequence + '\n')
    temp_fasta_file.close()

    align_sequences_from_fasta_file(temp_fasta_file.name, score_matrix_file,
                                    po_out_file, pir_out_file, clustal_out_file,
                                    do_global, do_progressive)

    os.unlink(temp_fasta_file.name)


def align_po_msas(po_msa_files,
                    score_matrix_file,
                    po_out_file,
                    pir_out_file=None,
                    clustal_out_file=None,
                    do_global=False,
                    do_progressive=True):
    """
    Align a list of po_msa files to a consensus po_msa
    """
    temp_list_file = tempfile.NamedTemporaryFile(delete=False)

    for po_msa_file in po_msa_files:
        temp_list_file.write(po_msa_file + '\n')
    temp_list_file.close()

    input_files_command = '-read_msa_list ' + temp_list_file.name
    _align(input_files_command, score_matrix_file, po_out_file, 
            pir_out_file, clustal_out_file, do_global, do_progressive)

    os.unlink(temp_list_file.name)


def convert_po_msa_to_dag(po_msa_file):
    """
    Convert a po_msa_file to a Directed Acyclic Graph

    The output is a list with each element being a dict containing keys:
        index(vertex id),
        incoming(list of incoming vertex ids)
        sequences(list of sequence ids contained in this vertex)
        charecter(the actual character here)
    """
    dag = []
    
    f = open(po_msa_file, 'r')
    i = 0
    for line in f:
        if '=' in line:
            continue
        char, rest = line.split(':')
        node_data = re.findall('([LS]\d+)', rest)
        node = {'index': i,
                'character': char,
                'incoming': [],
                'sequences': set()}
        for data in node_data:
            identifier, number = data[0], int(data[1:])
            if identifier == 'L':
                node['incoming'].append(number)
            elif identifier == 'S':
                node['sequences'].add(number)
        dag.append(node)
        i += 1

    f.close()
    
    return dag




