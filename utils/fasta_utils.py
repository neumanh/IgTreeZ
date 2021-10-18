#!/usr/bin/python3
"""
Contains functions for AIRR scheme json file processing
"""

# Info

__author__ = 'Hadas Neuman'
__version__ = '1.0'
__date__ = '19/10/20'

import multiprocessing as mp
# Imports
import os
from itertools import repeat  # For starmap

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from ete3 import PhyloTree, Tree
from objects import Logger
from utils import general_utils, tree_utils


def find_gl_seq_in_fasta(sequence_file: str, gl: str):
    """
    Find the germline sequencce in a Fasta file
    :param sequence_file: The file with the GL sequence
    :param gl: The GL name
    :return: The germline sequece
    """
    gl_seq = None
    for record in SeqIO.parse(sequence_file, "fasta"):
        if (record.name == gl) and (not gl_seq):
            gl_seq = str(record.seq)

    return gl_seq


def attach_gl_seqs(t: Tree, sequences_file: str, gl):
    """
    Attaches the germline sequences to nodes
    :param t: a Tree
    :param sequences_file: The sequences file
    :param gl: The GL name
    :return: The updated tree
    """
    for node in t.traverse("preorder"):
        if not hasattr(node, 'sequence'):  # The node is not linked to any sequence
            if hasattr(node, 'name') and (node.name == gl):
                gl_seq = find_gl_seq_in_fasta(sequences_file, gl)
                if gl_seq:
                    node.sequence = gl_seq
                    # print('Attaching gl to the node', node.name) # DEBUG

    return t


def process_sequences_names(sequence_file: str):
    """
    Processes the nodes names, for IgTree input files
    :param sequence_file: The sequences file
    :return: The output file names
    """
    out_file = os.path.basename(sequence_file) + "_names-updated.fasta"
    updt_records = []
    for record in SeqIO.parse(sequence_file, "fasta"):
        updated_record_name = record.id.replace(":", "-")
        updated_record_name = updated_record_name.replace(";", "-")
        record.id = updated_record_name
        new_record = SeqRecord(record.seq, updated_record_name, '', '')
        updt_records.append(new_record)

    # Writing the records to file
    SeqIO.write(updt_records, out_file, "fasta")
    return out_file


def get_seqs_dict(fasta_file: str):
    """
    Creates a dictionary of seq_id: sequence
    :param fasta_file: A fasta file.
    :return: A dictionary of sequence ids as keys, and sequences as values
    """
    seq_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_dict[record.name] = str(record.seq)
    return seq_dict


def get_seq_from_dict(node_name: str, seq_dict: dict):
    """
    Gets the sequence from the seq dictionary, if exists.
    :param node_name: The node name.
    :param seq_dict: The sequences dictionary
    :return: The corresponding sequence, or None if the sequence was not found.
    """
    seq = None
    for seq_name in seq_dict:
        if (seq_name in node_name) and (not seq):
            seq = seq_dict[seq_name]
    return seq


def link_alignment_to_tree(fasta_file: str, t: Tree, gl_name: str, log_object: Logger = None, illumina: str = False):
    """
    Links sequences to tree
    :param fasta_file: the IgTree input txt file
    :param t: An ETE3 Tree object
    :param gl_name: the GL name
    :param log_object: The Logger object
    :param illumina: Whether to replace colons with under-line
    :return: a tree that contains the aligned sequences
    """

    if not log_object:
        log_object = general_utils.create_general_log_object()

    # Saving the tree id
    tree_id = t.id

    if illumina:
        # Processing the fasta file
        temp_file = process_sequences_names(fasta_file)
    else:
        temp_file = fasta_file

    # Load the tree
    try:
        t = PhyloTree(t.write(format=1), alignment=temp_file, format=1)
        alright_flag = True
        t.id = tree_id
    except Exception as e:
        t = None
        alright_flag = False
        general_utils.handle_errors(f'Could not create a PhyloTree object for {tree_id}: {e}', log_object)

    if t:
        # Attaching the clone id to the tree
        t = general_utils.replace_new_line_in_ig_tree(t)

        # Attaching the GL sequence to the tree
        gl_seq = find_gl_seq_in_fasta(temp_file, gl_name)
        if gl_seq:
            t.sequence = gl_seq
        else:
            general_utils.handle_errors(f'Could not find a germline sequence named {gl_name} in the file {fasta_file}.',
                                        log_object)
            alright_flag = False

        # t = attach_gl_seqs(t, fasta_file, gl_name)

        # For backing up the linking
        seq_dict = get_seqs_dict(temp_file)
        print(seq_dict)  # TEMP

        if hasattr(t, 'sequence'):
            seq_num = 0  # The sequences counter
            for node in t.traverse('preorder'):  # Going from root to leaves to fill all the unknown sequences
                if hasattr(node, 'sequence'):
                    node.sequence = node.sequence.replace('.', '-').upper()
                    seq_num += 1  # A sequence

                # If the node is identical to its father
                elif (node.get_distance(node.up) == 0) and (hasattr(node.up, 'sequence')):
                    node.sequence = node.up.sequence.replace('.', '-').upper()

                elif hasattr(node, 'name'):
                    seq_form_dict = get_seq_from_dict(node.name, seq_dict)
                    if seq_form_dict:
                        node.sequence = seq_form_dict.replace('.', '-').upper()
                        seq_num += 1  # A sequence

            for node in t.traverse('postorder'):  # Going from leaves to root to generate missing sequences
                if not hasattr(node, 'sequence'):  # this in a hypothetical node w/o sequence
                    generated_seq = tree_utils.generate_missing_sequence(node.children, node.up, log_object)
                    if generated_seq:
                        node.sequence = generated_seq

                    else:
                        alright_flag = False

            t.seq_num = seq_num

    # Removing teh temporary fasta file
    if temp_file != fasta_file:
        # Removing the temporary file
        os.remove(temp_file)

    if not alright_flag:
        t = None

    return t


def attach_trees_using_fasta(args, log_object: Logger = None):
    """
    Attahces sequences to trees.
    :param args: The input arguments.
    :param log_object: The Logger object
    :return: A list of ETE trees
    """
    if not log_object:
        log_object = Logger(args.name, args.silent)

    # Collecting the fasta files
    if os.path.isdir(args.fasta[0]):
        seq_list = general_utils.get_files_from_dir(args.fasta[0])
    else:
        seq_list = args.fasta

    tree_list = general_utils.collect_trees(args, log_object)

    if len(tree_list) != len(seq_list):
        general_utils.handle_errors('The number of trees and fasta files is not equal.', log_object, exit_stat=True)

    # Linking the trees to sequences
    linked_trees = []

    number_of_workers = min(args.ncors, len(tree_list))
    if len(tree_list) > 0:
        with mp.Pool(number_of_workers) as pool:
            # Linking the trees with tha fasta files in prallel
            linked_trees = pool.starmap(link_alignment_to_tree,
                                        zip(seq_list, tree_list, repeat(args.gl_name), repeat(log_object),
                                            repeat(args.illumina)))

    return linked_trees
