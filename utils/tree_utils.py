#!/usr/bin/python3
"""
Functions for linking tree nodes and sequence alignments
"""

# Info

__author__ = 'Hadas Neuman'
__version__ = '1.0'
__date__ = '13/10/20'

import pandas as pd
# Imports
from utils import general_utils


def generate_missing_sequence(sons: list, father=None, log_object=None):
    """
    Generates sequence by finding the consensus for the father and sons sequences
    :param father: The parent of the node
    :param sons: The sons of the node
    :param log_object: The log object
    :return:
    """

    if not log_object:
        log_object = general_utils.create_general_log_object()

    # Collecting the references list
    son_seqs = []
    for son in sons:
        if hasattr(son, 'sequence'):
            son_seqs.append(son.sequence)

    if father and hasattr(father, 'sequence'):
        father_seq = father.sequence.upper()
        son_seqs.append(father_seq)

    cons_seq = None

    # Preform only if all the sequences have the same length:
    if len(son_seqs) > 0:
        # Generation the sequence
        cons_seq = ""
        if general_utils.is_equal_len(son_seqs):
            # for i, nuc in enumerate(son_seqs[0]):
            for i in range(0, len(son_seqs[0])):  # Going over the nucleotides in all the sequences
                # Finds all nucleotides in that position
                nucs_in_position = []
                for son in son_seqs:
                    nucs_in_position.append(son[i])

                cons_nuc = find_cons(nucs_in_position)
                cons_seq = cons_seq + cons_nuc
        else:
            general_utils.handle_errors("The sequences have different lengths. The tree could not be analysed.",
                                        log_object=log_object, exit_stat=False)

    return cons_seq


def find_cons(nucs: list):
    """
    Find the consensus in a single position
    :param nucs: The nucleotides list
    :return:
    """

    nuc_dic = {
        "-": 0,
        "N": 0,
        "A": 0,
        "C": 0,
        "G": 0,
        "T": 0,
    }
    for nuc in nucs:
        if nuc == "A":
            nuc_dic["A"] += 1
        elif nuc == "C":
            nuc_dic["C"] += 1
        elif nuc == "G":
            nuc_dic["G"] += 1
        elif nuc == "T":
            nuc_dic["T"] += 1
        elif nuc == "N":
            nuc_dic["N"] += 1
        elif nuc == "-":
            nuc_dic["-"] += 1

    # Return the nucleotide with the maximum value
    max_nuc = max(nuc_dic, key=nuc_dic.get)

    return max_nuc


def clean_db_from_none_vals(df: pd.DataFrame):
    """
    Cleans the mutations database from empty values
    :param df: A datafrmae
    :return: A clean dataframe
    """
    clean_df = pd.DataFrame()
    if len(df.columns) > 2:  # Not just the 'sample' and the 'tree_id' columns
        clean_df = df.dropna(thresh=3)  # Drop all the rows that do not have at least three values that are not NaN
    return clean_df
