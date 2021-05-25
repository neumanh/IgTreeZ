#!/usr/bin/python3

"""
Describes the numerical attributes of lineage trees in Newick format
"""

# Info
__author__ = 'Hadas Neuman'
__version__ = '1.5'
__date__ = '11/10/2020'

import multiprocessing as mp
from argparse import ArgumentParser

import numpy as np
import pandas as pd
# Project imports
import utils.general_utils as general_utils
from ete3 import Tree
from objects.logger import Logger
from objects.treat import Treat  # A tree-attribute object


def get_od_array(t: Tree):
    """
    Gets the out degree (the number of sons) of all nodes but the trunk (root's sons)
    :param t: The tree to analyse
    :return: The out degree list
    """

    leaves = t.get_leaves()

    # Initiate variables
    od_array = []
    dasn = []
    drsn = []
    nodes_num = 0

    for node in t.traverse():
        # -- For OD array --
        if node not in leaves:

            od_array.append(len(node.get_children()))  # The number of sons

            # TODO: Check if true
            # Ignore the trunk's length as the GL (root) is a speculation
            if not node.is_root() and (not node.up.is_root()):
                dist = int(node.get_distance(node.up))  # The distance to the father
                if dist > 1:
                    # temp_arr = np.repeat(1, (dist - 1))  # Create an array of ones in the length of the distance-1
                    temp_arr = [1 for _ in range(dist - 1)]  # Create an array of ones in the length of the distance-1
                    od_array += temp_arr

        if not node.is_root():
            # -- For nodes' number
            if hasattr(node, 'name'):
                nodes_num += 1

            # For DASN array -  the distances between adjacent split nodes
            dist = node.get_distance(node.up)
            dasn.append(dist)

            # For DRSN array - the distances from the root to any split node/fork
            num_sons = len(node.get_children())
            if num_sons > 1:
                dist = node.get_distance(t)
                drsn.append(dist)

    return nodes_num, np.array(od_array), np.array(dasn), np.array(drsn)


def get_trunk_len(t: Tree):
    """
    Finds the length of tree trunk from root to the first split node, that
    is the number of mutations shared by all leaves sequences found.
    :param t: An ETE tree
    :return: The trunk's length
    """

    first_split_node = None
    node = t

    while not first_split_node:
        num_sons = len(node.get_children())
        if num_sons == 1:  # There is only one son
            node = node.get_children()[0]  # Look for the next node
        else:
            first_split_node = node  # Stop the search in you got to a split node or to a leaf

    if first_split_node.is_root():  # The root splits
        dist_arr = []
        for child in first_split_node.get_children():
            dist_to_root = child.get_distance(t)
            # dist_arr.append(dist_to_root)
            dist_arr.append(dist_to_root)
        if dist_arr:
            dist_arr = np.array(dist_arr)
            trunk = np.mean(dist_arr)
        else:
            trunk = 0  # The root is all the tree
    else:
        trunk = first_split_node.get_distance(t)

    return trunk, first_split_node


def get_nodes_num(t):
    """
    Finds the number of nodes
    :param t: An ETE tree
    :return: The number of nodes
    """
    nodes_num = 0
    for node in t.traverse():
        if not node.is_root():
            if hasattr(node, 'name'):
                nodes_num += 1
    return nodes_num


def get_leaves_num(t):
    """
    Finds the number of leaves in the tree
    :param t: An ETE tree
    :return: The number of leaves
    """
    leaves = t.get_leaves()
    num_leaves = len(leaves)
    return num_leaves


def get_root_degree(t):
    """
    Finds the root degree
    :param t: An ETE tree
    :return: The root degree
    """
    rootd = len(t.get_children())
    return rootd


def get_pl_arrays(t: Tree, first_split_node: Tree):
    """
    Finde the PL's - path to leaves -of the tree
    :param t: An ETE tree
    :param first_split_node: The first split node
    :return: The PL array
    """
    leaves = t.get_leaves()

    pls_root = []
    pls_splitnode = []

    for leaf in leaves:
        dist_to_root = leaf.get_distance(t)  # Find the distance to the root
        dist_to_splitnode = leaf.get_distance(first_split_node)  # Find the distance to the the first split node
        pls_root.append(dist_to_root)
        pls_splitnode.append(dist_to_splitnode)
    return np.array(pls_root), np.array(pls_splitnode)


def initiate_results_dict(parameters: list):
    """
    Initiates the results dictionary.
    :param parameters: The keys for the results dictionary
    :return: The dictionary
    """
    res_dict = {}
    for key in parameters:
        res_dict[key] = []

    return res_dict


def get_one_tree_attr(t: Tree):
    """
    Gets one tree's attributes
    :param t: An ete tree
    :return: A dictionary based on a treat object
    """

    clone_id = general_utils.get_clone_number(general_utils.get_clone_name(t.id))

    if t and len(list(t)) > 0:  # If there is more than one node
        tr = Treat()  # A tree attributes object
        leaves = get_leaves_num(t)
        nodes, od_array, dasn_array, drsn_array = get_od_array(t)  # Half of the time
        trunk, first_split_node = get_trunk_len(t)
        pl_array, dlfsn_array = get_pl_arrays(t, first_split_node)  # Half of the time

        # Saving the relevant results
        if drsn_array.any():
            drsn_min = np.min(drsn_array)
        else:
            drsn_min = None

        # Update the treat object
        tr.id = clone_id
        tr.nodes = nodes
        tr.leaves = leaves
        tr.od_avg = np.mean(od_array)
        tr.rootd = get_root_degree(t)
        tr.drsn_min = drsn_min
        tr.dasn_min = np.min(dasn_array)
        tr.trunk = trunk
        tr.dlfsn_avg = np.mean(dlfsn_array)
        tr.pl_min = np.min(pl_array)

        tr = tr.get_dict()

    else:
        tr = None

    return tr


def mtree(args):
    """
    The main function.
    :param args: The input arguments
    :return:
    """

    log_object = Logger(args.name, args.silent, 'mtree')

    # Getting the input arguments
    tree_list = general_utils.collect_trees(args, log_object)
    log_object.info(f'Finish collecting {len(tree_list)} trees')

    args.ncors = int(args.ncors)
    number_of_workers = min(args.ncors, len(tree_list))
    with mp.Pool(number_of_workers) as pool:
        # Running poptree in parallel
        ret_list = pool.map(get_one_tree_attr, tree_list)

    log_object.info(f'Finish going over the trees')

    # Save to pandas dataframe
    if ret_list:
        df = pd.DataFrame(ret_list)
        df.insert(0, general_utils.sample_str, args.name)
        dir_name = general_utils.create_output_dir(args.name, sub_dir_name=None, log_object=log_object)
        table_file_name = f'{dir_name}/{args.name}_mtree.csv'
        df.to_csv(table_file_name, index=False)
        log_object.info(f'The output was saved to: {table_file_name}')

    log_object.done()


def get_arg_parser():
    """
    Defines the input and the output field help message
    :return: The parser
    """

    desc = "Collects lineage tree attributes"

    # Define argument parser
    local_parser = ArgumentParser(description=desc)
    local_parser = general_utils.add_basic_arguments(local_parser)
    local_parser.add_argument('--version', action='version', version='%(prog)s:' + ' %s:%s' % (__version__, __date__))
    return local_parser


if __name__ == "__main__":
    # Parse command line arguments
    parser = get_arg_parser()
    global_args = parser.parse_args()

    mtree(global_args)
