#!/usr/bin/python3

"""
Filters the trees by their populations
"""

# Info
import multiprocessing as mp
from argparse import ArgumentParser

from ete3 import Tree
from objects.logger import Logger
from utils import general_utils

__author__ = 'Hadas Neuman'
__version__ = '1.1'
__date__ = '2/11/20'


def get_pops_of_one_tree(t: Tree, pops: list):
    """
    Searches one tree for it's populations.
    :param t: An ETE3 tree
    :param pops: The population to search for
    :return: A set of the populations in the tree
    """
    pop_set = set([])
    if t:

        for node in t.traverse():
            if hasattr(node, 'name'):  # The the node has a name
                node_name = node.name
                for pop in pops:  # Look for all the populations in the node name
                    if pop in node_name:
                        pop_set.add(pop)
    return pop_set


def search_by_or(tree_dicts: dict, or_pops: list):
    """
    Collects the trees with either one of the given populations.
    :param tree_dicts: The tree-name: populations dictionary
    :param or_pops: The populations to search for
    :return: The relevant tree file list
    """
    filtered_trees = []
    for tree_file in tree_dicts:
        tree_pops = tree_dicts[tree_file]
        any_pops_flag = False
        for pop in or_pops:
            if pop in tree_pops:
                any_pops_flag = True  # Atleast one population presents
        if any_pops_flag:
            filtered_trees.append(tree_file)
    return filtered_trees


def search_by_not(tree_dicts: dict, not_pops: list):
    """
    Collects the trees with all of the given populations.
    :param tree_dicts: The tree-name: populations dictionary
    :param not_pops: The populations to search for
    :return: The relevant tree file list
    """
    filtered_trees = []
    for tree_file in tree_dicts:
        tree_pops = tree_dicts[tree_file]
        none_pops_flag = True
        for pop in not_pops:
            if pop in tree_pops:
                none_pops_flag = False  # At lest One pop is in the tree
        if none_pops_flag:
            filtered_trees.append(tree_file)

    return filtered_trees


def search_by_and(tree_dicts: dict, and_pops: list):
    """
    Collects the trees with all of the given populations.
    :param tree_dicts: The tree-name: populations dictionary
    :param and_pops: The populations to search for
    :return: The relevant tree file list
    """
    filtered_trees = []
    for tree_file in tree_dicts:
        tree_pops = tree_dicts[tree_file]
        all_pops_flag = True
        for pop in and_pops:
            if pop not in tree_pops:
                all_pops_flag = False  # At least One pop is not in the tree
        if all_pops_flag:
            filtered_trees.append(tree_file)

    return filtered_trees


def generate_out_file_name(pop_list: list, delimited: str):
    """
    Generates output file name.
    :param pop_list: The populations list.
    :param delimited: The delimited to add between the populations.
    :return: The output file name.
    """

    if delimited == 'nodes':
        filename = f'more_than_{pop_list[0]}_nodes'
        if len(pop_list) > 1:
            filename += f'_and_less_than_{pop_list[1]}'
        filename += '_nodes'
    elif delimited == 'leaves':
        filename = f'more_than_{pop_list[0]}_leaves'
        if len(pop_list) > 1:
            filename += f'_and_less_than_{pop_list[1]}'
        filename += '_nodes'
    else:  # or/and/not
        if delimited == 'not':
            filename = 'not_'
        else:  # or/and
            filename = ''
        pop_len = len(pop_list)
        idx = 0
        for pop in pop_list:
            filename += pop  # Append the populatio name to the file name
            idx += 1
            if idx < pop_len:
                filename += f'_{delimited}_'  # If there are more population to add, add the delimited
    return filename


def save_trunkless_trees_to_file(filtered_tree_list: list, output_dir_name: str = '.',
                                 log_object: Logger = None):
    """
    Saves the file list to the given output file and copy the tree files if needed.
    :param filtered_tree_list: The tree list to save
    :param output_dir_name: The output directory
    :param log_object: The Logger object
    :return: The output file name and the output directory (if the copy_flag was used)
    """

    if not log_object:
        log_object = general_utils.create_general_log_object()

    general_utils.create_dir_if_needed(output_dir_name)
    for t in filtered_tree_list:
        if t:
            tree_file_name = f'{output_dir_name}/{t.id}.nw'
            try:
                t.write(format=1, outfile=tree_file_name)
            except IOError:
                general_utils.handle_errors(
                    f'Could not write the tree to the file {tree_file_name} in the directory {output_dir_name}',
                    exit_stat=True)
    return output_dir_name


def remove_trunk_from_tree(t):
    """
    Removes all the node until the first split node. If there is no trunk - returns None
    :param tree: An ETE tree
    :return: The tree w/o trunk or Null
    """
    split_node = None

    if len(t.children) < 2:
        for n in t.traverse('levelorder'):  # By levels
            #if hasattr(n,'name'):
                # print(f'node {n.name} has {len(n.children)} sons')  # TEMP
            if (not split_node) and (len(n.children) > 1):
                split_node = n

        #print('The new root should be', split_node.name)  # TEMP

    #print(split_node)  # TEMP

    if split_node:
        split_node.id = t.id

    return split_node


def filter_rep(args):
    """
    The main function - filter the repetoire according to given filter parameters
    :param args: The input arguments
    :return: None
    """

    log_object = Logger(args.name, args.silent, 'trunk-removing')
    # Collecting the trees
    tree_list = general_utils.collect_trees(args, log_object)
    args.ncors = int(args.ncors)
    number_of_workers = min(args.ncors, len(tree_list))
    with mp.Pool(number_of_workers) as pool:
        # Running poptree in parallel
        ret_list = pool.starmap(remove_trunk_from_tree, zip(tree_list))
    ret_list = [t for t in ret_list if t]  # Removing the None objects
    # _=[print(t.id) for t in ret_list]  # TEMP

    log_object.info(f'Finish going over {len(tree_list)} trees')

    # Saving the files
    output_dir = general_utils.create_output_dir(args.name, 'Trunkless_trees', log_object)

    out_dir = save_trunkless_trees_to_file(ret_list, output_dir, log_object)

    # Saying goodbye
    if out_dir:
        log_object.info(f'The {len(ret_list)} trunkless trees were copied to: {out_dir}')
    else:
        log_object.warning('No trees were found.')
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
    # local_parser.add_argument('-AND', action='store', nargs='+', help='Collect trees with all the given populations.')
    # local_parser.add_argument('-OR', action='store', nargs='+',
    # help='Collect trees with one of the given populations.')
    # local_parser.add_argument('-NOT', action='store', nargs='+',
    #                           help='Collect trees without the given populations.')
    # local_parser.add_argument('--copy', action='store_true',
    #                           help='If given, copies the collected trees in addition to listing them.')

    local_parser.add_argument('--version', action='version', version='%(prog)s:' + ' %s:%s' % (__version__, __date__))
    return local_parser


if __name__ == "__main__":
    # Parse command line arguments
    parser = get_arg_parser()
    global_args = parser.parse_args()

    filter_rep(global_args)
