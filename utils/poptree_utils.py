#!/usr/bin/python3
"""
Counts the mutations in Newick trees.
"""

# Info
__author__ = 'Hadas Neuman'
__version__ = '1.0'
__date__ = '13/10/20'

import re

# Imports
from ete3 import Tree
from objects.node import Node


def sort_pops(pop_dict: dict, input_pops: list):
    """
    Sorts the populations dictionary by the input order.
    :param pop_dict: A dictionary to sort
    :param input_pops: AThe input population names
    :return: The sorted dictionary
    """
    sorted_dict = {}
    for key in input_pops:
        if key in pop_dict:
            sorted_dict[key] = pop_dict[key]
    return sorted_dict


def sort_trans(trans_dict: dict, input_pops: list):
    """
    Sorts the transitions dictionary by the input order.
    :param trans_dict: A dictionary to sort
    :param input_pops: AThe input population names
    :return: The sorted dictionary
    """
    sorted_dict = {}
    for key1 in input_pops:
        for key2 in input_pops:
            if (key1 != key2) and ((key1, key2) in trans_dict):
                sorted_dict[(key1, key2)] = trans_dict[(key1, key2)]
    return sorted_dict


def replace_under_line_in_trans(trans_dict: dict):
    """
    Replaces the prefix under line, if exists.
    :param trans_dict: The transitions dictionary
    :return: The updated transitions dictionary
    """
    nu_dict = {}
    for two_pops in trans_dict:
        pop1 = two_pops[0]
        pop2 = two_pops[1]
        pop1_nu = re.sub(r'^_', '', pop1)
        pop2_nu = re.sub(r'^_', '', pop2)
        nu_dict[(pop1_nu, pop2_nu)] = trans_dict[(pop1, pop2)]
    return nu_dict


def replace_under_line_in_pops(pop_dict: dict):
    """
    Replaces the prefix under line, if exists.
    :param pop_dict: The populations dictionary
    :return: The updated transitions dictionary
    """
    nu_dict = {}
    for pop in pop_dict:
        pop_nu = re.sub(r'^_', '', pop)
        nu_dict[pop_nu] = pop_dict[pop]
    return nu_dict


def sum_dict(dic: dict):
    """
    Sums the values of a dictionary
    :param dic: A dictionary to summarize
    :return: The sum
    """
    count = 0
    # print(dic)
    for key in dic:
        count += (dic[key].flatten().sum())
    return count


def count_dict(dic: dict):
    """
    Sums the values of a dictionary
    :param dic: A dictionary to summarize
    :return: The sum
    """
    count = 0
    # print(dic)
    for key in dic:
        count += len(dic[key])
    return count


def count_pops_num_in_tree(pop_df):
    """
    Counts the number of populations in the trees
    :param pop_df: The populations-in-tree dataframe
    :return: A list of number of populations in trees
    """
    pop_count = []

    # Collecting the number of different populations

    for index, row in pop_df.iterrows():
        pop_count.append(len(row.dropna()) - 1)

    return pop_count


def normalize_by_the_source(trans_dict: dict, pops_dict: dict):
    """
    Normalizes the transition dictionary by the population dictionary
    :param trans_dict: The transition dictionary
    :param pops_dict: The population dictionary
    :return: A normalized transition dictionary
    """
    norm_dict = {}
    all_trans_count = count_dict(trans_dict)
    # Going over the transitions dictionary
    for trans_key in trans_dict:
        source = trans_key[0]  # The source of the transition
        if source in pops_dict:
            pop_size = len(pops_dict[source])  # The corresponding population size
            norm_dict[trans_key] = (len(trans_dict[trans_key]) / pop_size)
            # norm_dict[trans_key] *= all_trans_size # To work with larger numbers
            norm_dict[trans_key] *= (all_trans_count*10)  # To work with larger numbers

    return norm_dict


def normalize_by_the_dest(trans_dict: dict, pops_dict: dict):
    """
    Normalizes the transition dictionary by the population dictionary
    :param trans_dict: The transition dictionary
    :param pops_dict: The population dictionary
    :return: A normalized transition dictionary
    """
    norm_dict = {}
    all_trans_count = count_dict(trans_dict)
    # Going over the transitions dictionary
    for trans_key in trans_dict:
        dest = trans_key[1]  # The source of the transition
        if dest in pops_dict:
            pop_size = len(pops_dict[dest])  # The corresponding population size
            norm_dict[trans_key] = (len(trans_dict[trans_key]) / pop_size)
            # norm_dict[trans_key] *= all_trans_size # To work with larger numbers
            norm_dict[trans_key] *= (all_trans_count*10)  # To work with larger numbers

    return norm_dict


def get_path_list(tree: Tree):
    """
    Creates a list of all paths in the tree - from leaf to root
    :param tree: a Tree object
    :return: List of paths
    """

    # t1 is always the bigger tree
    rootname = tree.name
    all_pathes = []
    for leaf in tree.get_leaves():
        temp_path = []
        node = leaf
        while node.name != rootname:
            # node_object = Node(node.name, node.get_distance(node.up.name))
            node_object = node
            temp_path = temp_path + [node_object]
            node = node.up
        # Add root
        temp_path = temp_path + [Node(rootname, 0)]
        all_pathes.append(temp_path)
    return all_pathes


def sum_path(path: list, index: int = 0):
    """
    Sums the path's distances
    :param path: A list of nodes
    :param index: An index
    :return: The sum of the path's distance
    """
    # Initializing the counters
    sum_dist = 0
    idx = len(path) - 1
    # Summing the distances
    while idx >= index:
        sum_dist = sum_dist + path[idx].dist
        idx = idx - 1
    return sum_dist


def multiple_pops_in_name(name: str, pops: list):
    """
    Checks if more than on populations are in the tree name
    :param name: The node name
    :param pops: Thr populations
    :return: True if more than one populations are in the node name, and False otherwise.
    """
    first_pop, second_pop = False, False

    for pop in pops:
        if pop in name:
            if not first_pop:
                first_pop = True
            elif first_pop and (not second_pop):
                second_pop = True
    return second_pop
