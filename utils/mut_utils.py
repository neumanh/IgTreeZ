#!/usr/bin/python3
"""
Functions for mutations characterization
"""

# Info
__author__ = 'Hadas Neuman'
__version__ = '1.2'
__date__ = '30/9/20'

# Imports
import re

from ete3 import PhyloTree
from objects import Mutation, Logger
from utils import general_utils

gls = ["GL", "G.L", 'G.L.', "Germline", "Germline_UCA", 'EX;G.L', 'EX_G.L']


def replace_new_line_in_ig_tree(t: PhyloTree):
    """
    Removes everything after the \n charecter is sequence names
    :param t: A tree
    :return: The updated tree
    """
    for node in t.traverse("postorder"):  # Going from leaves to root to generate missing sequences
        if hasattr(node, 'name'):
            seq_name = node.name
            if '\\n' in seq_name:
                seq_name = re.sub(r'\\n.*', '', seq_name)
                node.name = seq_name
    return t


def find_first_non_gap_pos(seq):
    """
    Find the first position with a nucleotide, for finding the correct reading frame
    :param seq: The sequence
    :return: The first non-gap position
    """
    i = 0
    nuc = seq[i]

    while ((nuc == "-") or (nuc == ".")) and (i < len(seq)):
        i = i + 1
        nuc = seq[i]

    return i


def is_regular_nuc(nuc: str):
    """
    Check whether the input string is a single regular nucleotide
    :param nuc: a nucleotide string
    :return: True if it is A/C/G/T or False otherwise
    """
    nuc = nuc.upper()
    if (nuc == 'A') or (nuc == 'C') or (nuc == 'G') or (nuc == 'T'):
        real_nuc = True
    else:
        real_nuc = False
    return real_nuc


def get_rf_start(gl_start: int, seq_start: int):
    """
    Finds the reading frame start position
    :param gl_start: The germline first non-gap position
    :param seq_start: The sequence first non-gap position
    :return: The reading frame start position
    """
    if gl_start >= seq_start:
        rfs = gl_start
    else:
        while ((seq_start - gl_start) % 3) > 0:
            seq_start += 1
        rfs = seq_start

    return rfs


def count_observed_in_tree(t, no_trunk: bool = False, no_cdr3: bool = False, log_object: Logger = None):
    """
    Counts the number of observed mutations in a tree
    :param t: The tree
    :param no_trunk: Indication whether to analyze root's children
    :param no_cdr3: Indication whether to analyze the CDR3 region
    :param log_object: The Logger object
    :return: The mutation dictionary of one tree
    """
    dist_flag = True

    if not log_object:
        log_object = general_utils.create_general_log_object()

    mut_dict = {}
    if t:
        tree_dist = 0
        clone_id = t.id
        gl_seq = t.sequence
        gl_start_pos = find_first_non_gap_pos(gl_seq)  # For the correct reading frame

        root = t

        # For mutation analysis
        new_mut_list = []

        # Going through the tree from root to leaves
        for node in t.traverse("preorder"):
            if node != root:  # This is not the root
                if dist_flag:
                    tree_dist += node.get_distance(node.up)
                # If the no_trunk flag is on, analyze nodes which are not children fo the root
                if (not no_trunk) or (node.up != root):
                    seq = node.sequence.upper()
                    parent_seq = node.up.sequence.upper()
                    if parent_seq != seq:
                        # print('\nFrom node', node.up.name, 'to', node.name, ':') # TEMP
                        regions = get_regions_by_tree(t)
                        temp_mut_list = count_observed_in_node(parent_seq, seq, regions, gl_start_pos, no_cdr3,
                                                               log_object=log_object)
                        if temp_mut_list:
                            new_mut_list += temp_mut_list
                            # print('Between node', node.name, 'and node', node.up.name, 'found', len(temp_mut_list),
                            #       'mutations')  # TEMP
        mut_dict = sum_mut_for_one_tree(new_mut_list, clone_id)
        mut_dict['clonal_germline'] = gl_seq
        if dist_flag:
            mut_dict['all_muts_by_tree_topology'] = tree_dist

        # Add the node count
        mut_dict['nodes'] = general_utils.count_tree_nodes(t)
        mut_dict['sequences'] = t.seq_num
        if hasattr(t, 'cdr3_end') and (not no_cdr3):
            mut_dict['cdr3_end'] = t.cdr3_end

    return mut_dict


def count_observed_in_node(gl: str, seq: str, regions: dict, gl_start_pos=0, nocdr3: bool = False,
                           log_object: Logger = None):
    """
    Counts the number of observed mutations. Returns dictionary of position:R/S
    :param gl: The parent sequence
    :param seq: The son sequence
    :param regions: The FWR-CDR starting posiotions
    :param gl_start_pos: The start position of the GL reading frame
    :param nocdr3: Avoid from defining the CDR3 region
    :param log_object: The logger object
    :return: A dictionary of position:R/S
    """

    if not log_object:
        log_object = general_utils.create_general_log_object()

    mut_arr = []
    if not general_utils.is_equal_len([gl, seq]):
        general_utils.handle_errors("The sequences must have the same lengths", log_object)

    else:
        start_pos = max(get_rf_start(gl_start_pos, find_first_non_gap_pos(seq)),
                        get_rf_start(gl_start_pos, find_first_non_gap_pos(gl)))
        non_gapped_i = start_pos
        gaps = 0
        for i in range(start_pos, len(gl)):

            if gl[i] != seq[i]:
                if (is_regular_nuc(gl[i])) and (is_regular_nuc(seq[i])):  # A mutation

                    # Getting the gl codon
                    start_i = i - ((i - gaps - start_pos) % 3)
                    codon_gl = gl[start_i:(start_i + 3)]

                    # Getting the sequence codon with only one change
                    if (start_i + 2) < len(seq):
                        if (is_regular_nuc(seq[start_i])) and \
                                (is_regular_nuc(seq[start_i + 1])) and \
                                (is_regular_nuc(seq[start_i + 2])):  # A regular codon
                            temp_gl = gl[:i] + seq[i] + gl[i + 1:]
                            codon_seq = temp_gl[start_i:(start_i + 3)]
                        else:
                            # Don't replace with GL nucleotide
                            codon_seq = seq[start_i] + seq[start_i + 1] + seq[start_i + 2]
                    else:  # End of string
                        codon_seq, codon_gl = None, None

                    # print('codon_gl:', codon_gl, '\tcodon_seq:', codon_seq)  # TEMP
                    mut = Mutation()
                    mut.describe_mut_by_nuc(gl[i], seq[i])
                    if codon_gl and codon_seq:
                        mut.describe_mut_by_codon(codon_gl, codon_seq)
                    mut.define_mut_pos(i, regions, nocdr3)

                    # print(mut) # TEMP

                    if mut:
                        mut_arr.append(mut)
            if gl[i] != '-':  # Not a gap
                non_gapped_i += 1
            if gl[i] == '-':  # Not a gap
                gaps += 1

    return mut_arr


def sum_mut_for_one_tree(mut_array: list, clone_id: str):
    """
    Sums the mutation for one tree
    :param clone_id: The tree id (a sting)
    :param mut_array: A list of mutation object
    :return: A summary dictionary
    """

    dic = {
        general_utils.id_str: clone_id,
        'clone_id': general_utils.get_clone_number(clone_id),
        'CDR': 0,
        'FWR': 0,
        'hydro_change': 0,
        'charge_change': 0,
        'hydro_keep': 0,
        'charge_keep': 0,
        'source_neutral': 0,
        'dest_neutral': 0,
        'mu_count_cdr_r': 0,
        'mu_count_cdr_s': 0,
        'mu_count_fwr_r': 0,
        'mu_count_fwr_s': 0,
        'all_mutations': len(mut_array)
    }
    for mut in mut_array:  # For each mutaion
        temp_mut_dict = mut.get_dict()
        for key in temp_mut_dict:  # For each attribute
            if temp_mut_dict[key]:  # If the attribute is True
                if key in dic:
                    dic[key] += 1
                else:
                    dic[key] = 1
        if (mut.source_hydrophilic == mut.dest_hydrophilic) and (mut.source_hydrophobic == mut.dest_hydrophobic):
            dic['hydro_keep'] += 1
        else:
            dic['hydro_change'] += 1
        if (mut.source_positive == mut.dest_positive) and (mut.source_negative == mut.source_negative):
            dic['charge_keep'] += 1
        else:
            dic['charge_change'] += 1
        if (not mut.source_positive) and (not mut.source_negative):
            dic['source_neutral'] += 1
        if (not mut.dest_positive) and (not mut.dest_negative):
            dic['dest_neutral'] += 1

        if mut.cdr1 or mut.cdr2 or mut.cdr3:
            dic['CDR'] += 1
        elif mut.fwr1 or mut.fwr2 or mut.fwr3:
            dic['FWR'] += 1

        # Selection fields
        if mut.r and (mut.cdr1 or mut.cdr2 or mut.cdr3):  # Including the CDR3 region and later compatibility to shazam
            dic['mu_count_cdr_r'] += 1
        if mut.s and (mut.cdr1 or mut.cdr2):
            dic['mu_count_cdr_s'] += 1
        if mut.r and (mut.fwr1 or mut.fwr2 or mut.fwr3):
            dic['mu_count_fwr_r'] += 1
        if mut.s and (mut.fwr1 or mut.fwr2 or mut.fwr3):
            dic['mu_count_fwr_s'] += 1

    return dic


def get_next_codon(seq: str, start_idx: int):
    """
    Gets the next three nucleotides, ignoring gaps.
    :param seq: The nucleotide sequence
    :param start_idx: The start index of the codon
    :return: The codon and the next start-index
    """
    if (start_idx + 2) < len(seq):
        new_index = start_idx + 3
        codon = seq[start_idx:(start_idx + 3)]  # Get the first codon
        gap_index = codon.find('-')  # Find the first gap
        i = start_idx
        while (gap_index >= 0) and ((i + 3) < len(seq)):
            gi = gap_index + i  # Get the adjusted gap index
            codon = seq[i:gi] + seq[(gi + 1):(i + 4)]
            gap_index = codon.find('-')
            i += 1
            new_index = i + 3

    else:
        new_index = -1
        codon = None
    return codon, new_index


def get_regions_by_tree(t):
    """
    Updates the region dictionary by the tree data.
    :param t: An ETE3 tree
    :return: The updated region locations as a dictionary
    """
    temp_regions = general_utils.imgt_regions

    if hasattr(t, 'fwr1'):
        temp_regions['fwr1'] = t.fwr1
    if hasattr(t, 'cdr1'):
        temp_regions['cdr1'] = t.cdr1
    if hasattr(t, 'fwr2'):
        temp_regions['fwr2'] = t.fwr2
    if hasattr(t, 'cdr2'):
        temp_regions['cdr2'] = t.cdr2
    if hasattr(t, 'fwr3'):
        temp_regions['fwr3'] = t.fwr3
    if hasattr(t, 'cdr3'):
        temp_regions['cdr3'] = t.cdr3
    if hasattr(t, 'cdr3_end'):
        temp_regions['cdr3_end'] = t.cdr3_end
    elif hasattr(t, 'cdr3_length'):
        temp_regions['cdr3_end'] = temp_regions['cdr3'] + temp_regions['cdr3_length']

    return temp_regions
