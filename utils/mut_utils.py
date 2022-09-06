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
import numpy as np

from ete3 import PhyloTree
from objects import Mutation, Logger
from utils import general_utils

gls = ["GL", "G.L", 'G.L.', "Germline", "Germline_UCA", 'EX;G.L', 'EX_G.L']
all_mut_field = 'all_mutations'
mut_per_seq_field = 'mutations_per_sequence'
seq_field = 'sequences'
all_fields = [general_utils.sample_str, general_utils.id_str, 'nodes', seq_field,
              'source_a', 'source_c', 'source_g', 'source_t', 'transition', 'transversion',
              'replacement', 'silent',
              'cdr1', 'cdr2', 'cdr3', 'CDR', 'fwr1', 'fwr2', 'fwr3', 'FWR',
              'positive_source', 'negative_source', 'neutral_charge_source',
              'positive_target', 'negative_target', 'neutral_charge_target', 'charge_keep', 'charge_change',
              'hydrophobic_source', 'hydrophilic_source', 'neutral_hydro_source',
              'hydrophobic_target', 'hydrophilic_target', 'neutral_hydro_target', 'hydro_keep', 'hydro_change',
              'vs_volume_source', 'small_volume_source', 'medium_volume_source', 'large_volume_source',
              'vl_volume_source', 'vs_volume_target', 'small_volume_target', 'medium_volume_target',
              'large_volume_target', 'vl_volume_target', 'volume_decrease', 'volume_increase', 'volume_keep',
              'amide_source', 'acidic_source', 'basic_source', 'hydroxyl_source', 'sulfur_source', 'aromatic_source',
              'aliphatic_source', 'amide_target', 'acidic_target', 'basic_target', 'hydroxyl_target', 'sulfur_target',
              'aromatic_target', 'aliphatic_target', 'chemical_keep', 'chemical_change',
              'hydro_donor_source', 'hydro_acceptor_source', 'hydro_donor_acceptor_source', 'hydro_da_none_source',
              'hydro_donor_target', 'hydro_acceptor_target', 'hydro_donor_acceptor_target', 'hydro_da_none_target',
              'hydro_da_keep', 'hydro_da_change', 'polar_source', 'non_polar_source', 'polar_target',
              'non_polar_target', 'polarity_keep', 'polarity_change',
              'average_aa_distance', 'median_aa_distance', 'overall_aa_distance',
              'average_replacement_distance', 'median_replacement_distance', 'overall_replacement_distance',
              'mutations_per_sequence', 'all_muts_by_tree_topology', 'cdr3_length', all_mut_field]


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
    dist_flag = True  # Count the mutations by tree topology

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
        new_dist_list = []

        regions = get_regions_by_tree(t)

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
                        temp_mut_list, distance_list = count_observed_in_node(parent_seq, seq, regions, gl_start_pos,
                                                                              no_cdr3, log_object=log_object)
                        if temp_mut_list:
                            new_mut_list += temp_mut_list
                        if distance_list:
                            new_dist_list += distance_list
                        if general_utils.debug_mode > 1:
                            print(temp_mut_list)
                        
        mut_dict = sum_mut_for_one_tree(new_mut_list, clone_id)
        if general_utils.debug_mode:
            print(mut_dict)
        if dist_flag:
            mut_dict['all_muts_by_tree_topology'] = tree_dist

        # Add more counts
        mut_dict['clonal_germline'] = gl_seq
        mut_dict['nodes'] = general_utils.count_tree_nodes(t)
        mut_dict['sequences'] = t.seq_num

        if general_utils.debug_mode:
            print(mut_dict)

        # Add distance counts
        if len(new_dist_list) > 0:
            mut_dict['average_aa_distance'] = np.round(np.average(new_dist_list), 2)
            mut_dict['median_aa_distance'] = np.median(new_dist_list)
            mut_dict['overall_aa_distance'] = np.sum(new_dist_list)

        # Only for replacement mutation
        positive_dist_list = [i for i in new_dist_list if i != 0]
        if len(positive_dist_list) > 0:
            mut_dict['average_replacement_distance'] = np.round(np.average(positive_dist_list), 2)
            mut_dict['median_replacement_distance'] = np.median(positive_dist_list)
            mut_dict['overall_replacement_distance'] = np.sum(positive_dist_list)

        # Add CDR3 length
        if hasattr(t, 'cdr3_end') and (not no_cdr3):
            mut_dict['cdr3_end'] = t.cdr3_end
            mut_dict['cdr3_length'] = t.cdr3_end - regions['cdr3']

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
    dist_arr = []
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
                    if general_utils.debug_mode > 2:
                        print('codon_gl:', codon_gl, '\tcodon_seq:', codon_seq)  # TEMP
                    mut = Mutation(codon_gl, codon_seq, i, regions, nocdr3)

                    if mut:
                        mut_arr.append(mut.get_one_mut_list())
                        dist = mut.define_distance()
                        # TODO: should I add zero distance?
                        if dist >= 0:  # Else, no distance was found
                            dist_arr.append(dist)

            if gl[i] != '-':  # Not a gap
                non_gapped_i += 1
            if gl[i] == '-':  # Not a gap
                gaps += 1

    if general_utils.debug_mode > 1:
        print(mut_arr, dist_arr)

    return mut_arr, dist_arr


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
        'mu_count_cdr_r': 0,  # Selection fields
        'mu_count_cdr_s': 0,
        'mu_count_fwr_r': 0,
        'mu_count_fwr_s': 0
        # 'transition_rep': 0,
        # 'transversion_rep': 0
    }

    # Initializing all mutation fields
    for field in all_fields:
        if field != general_utils.id_str:
            dic[field] = 0

    for mut in mut_array:  # For each mutation
        dic[all_mut_field] += 1
        for field in all_fields:
            if field in mut:
                dic[field] += 1
        if ('CDR' in mut) and ('replacement' in mut):
            dic['mu_count_cdr_r'] += 1
        elif ('CDR' in mut) and ('silent' in mut):
            dic['mu_count_cdr_s'] += 1
        elif ('FWR' in mut) and ('replacement' in mut):
            dic['mu_count_fwr_r'] += 1
        elif ('FWR' in mut) and ('silent' in mut):
            dic['mu_count_fwr_s'] += 1

        # Transition/transversion - replacement
        # if ('transition' in mut) and ('replacement' in mut):
        #     dic['transition_rep'] += 1
        # elif ('transversion' in mut) and ('replacement' in mut):
        #     dic['transversion_rep'] += 1

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
