#!/usr/bin/python3
"""
Contains functions for AIRR scheme json file processing
"""

# Info
__author__ = 'Hadas Neuman'
__version__ = '1.0'
__date__ = '19/10/20'

# Imports
import json

from objects import Logger
from utils import general_utils, tree_utils


def get_tree(clone_json_object: dict, log_object: Logger):
    """
    Gets a tree from a clone json object
    :param clone_json_object: The dictionary of the clone json object
    :param log_object: The Logger object
    :return: An ETE tree
    """
    clone_id = clone_json_object['clone_id']
    trees = clone_json_object['trees']
    tree_str = trees[0]['newick']  # Looking at the first tree # TODO: Choose the correct tree
    t = general_utils.create_tree(tree_str, log_object)
    t.id = clone_id

    return t


def collect_trees_from_scheme(json_file, log_object: Logger = None):
    """
    Collects the Newick trees from the scheme file and returns a list of trees.
    :param json_file: An AIRR scheme file
    :param log_object: The Logger object
    :return: A list of trees. Each tree hols the clone id in as t.id
    """

    if not log_object:
        log_object = general_utils.create_general_log_object()

    with open(json_file) as f:
        json_object = json.load(f)

    tree_list = []
    # Going over the input dictionary
    for clone in json_object['clones']:
        t = get_tree(clone, log_object)
        if t:
            tree_list.append(t)
    return tree_list


def get_nodes_dict(nodes_in_tree_obj):
    """
    Created a dictionary of node_id: seq
    :param nodes_in_tree_obj: A dictionary of the tree object
    :return: A dictionary of nodes an seqs
    """
    nodes_dict = {}
    for js_node in nodes_in_tree_obj['nodes']:
        node_id = js_node  # The primary key
        seq = nodes_in_tree_obj['nodes'][js_node]['sequence_alignment']
        nodes_dict[node_id] = seq
    return nodes_dict


def get_tree_and_seqs(clone_json_object: dict, log_object: Logger):
    """
    Gets a tree from a clone json object
    :param clone_json_object: A dictionary of one clone.
    :param log_object: The Logger object.
    :return: An ETE tree
    """
    clone_id = clone_json_object['clone_id']
    gl = clone_json_object['germline_alignment']
    alright_flag = True

    trees = clone_json_object['trees']

    tree_json_obj = trees[0]  # Looking at the first tree # TODO: Choose the correct tree

    # Creating a tree object
    # print(clone_id, gl, tree_json_obj['newick'], '\n\n\n')

    t = general_utils.create_tree(tree_json_obj['newick'], log_object)

    if t:

        t.id = clone_id
        t.sequence = gl

        # Attaching the sequences
        nodes = get_nodes_dict(tree_json_obj)

        seq_num = 0  # SThe sequences counter
        # Going from up to down
        for node in t.traverse("preorder"):
            if hasattr(node, 'name') and (node.name in nodes):
                node.sequence = nodes[node.name].replace('.', '-').upper()
                seq_num += 1
            # If the node is identical to its father
            elif node != t:  # This is not the root
                if (node.get_distance(node.up) == 0) and (hasattr(node.up, 'sequence')):
                    node.sequence = node.up.sequence.replace('.', '-').upper()

        # Generating the missing sequences
        for node in t.traverse("postorder"):  # Going from leaves to root
            if not hasattr(node, 'sequence'):
                generated_seq = tree_utils.generate_missing_sequence(node.children, node.up, log_object)
                if generated_seq:
                    node.sequence = generated_seq
                else:
                    log_object.warning(f'Could not generate sequence for {node.name}')
                    alright_flag = False
        t.seq_num = seq_num

    else:
        log_object.error(f'Could not create the clone {clone_id}')
        alright_flag = False

    if not alright_flag:
        t = None
    return t


def collect_trees_and_seqs_from_scheme(args, log_object: Logger):
    """
    Collects the Newick trees and corresponding sequences from the scheme file and returns a list of trees
    (linked to sequences).
    :param args: The input arguments, including the AIRR scheme file
    :param log_object: The Logger object.
    :return: A list of trees. Each tree hols the clone id in as t.id
    """
    json_file = args.json

    if not log_object:
        log_object = general_utils.create_general_log_object()

    if args.illumina or args.igtree:
        log_object.warning('The arguments \'--illumina\' and \'--igtree\' are ignored for AIRR scheme analysis')

    with open(json_file) as f:
        json_object = json.load(f)

    tree_list = []
    # Going over the input dictionary
    for clone in json_object['clones']:
        t = get_tree_and_seqs(clone, log_object)
        if t:
            tree_list.append(t)
    return tree_list
