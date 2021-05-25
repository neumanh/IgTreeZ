#!/usr/bin/python3

"""
Filters the trees by their populations
"""

# Info
import multiprocessing as mp
from argparse import ArgumentParser
from itertools import repeat  # For starmap

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


def attach_result_to_tree(tree_list, pop_results):
    """
    Attaches trees to starmap results
    :param tree_list: The tree list
    :param pop_results: The strmap results
    :return: A dictionary of tree_id: strmap result
    """
    res_dic = {}
    for idx, t in enumerate(tree_list):
        res_dic[t.id] = pop_results[idx]
    return res_dic


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


def save_filtered_trees_to_file(filtered_tree_list: list, pop_list: list, delimiter: str, tree_dict: dict,
                                output_dir: str = '.', copy_flag: bool = False, log_object: Logger = None):
    """
    Saves the file list to the given output file and copy the tree files if needed.
    :param filtered_tree_list: The tree list to save
    :param pop_list: The output file
    :param delimiter: The delimiter to add to the output name
    :param tree_dict: The tree_id: tree dictionary
    :param output_dir: The output directory
    :param copy_flag: Whether to copy the filtered files
    :param log_object: The Logger object
    :return: The output file name and the output directory (if the copy_flag was used)
    """

    if not log_object:
        log_object = general_utils.create_general_log_object()

    base_out_name = generate_out_file_name(pop_list, delimiter)
    out_file = f'{output_dir}/{base_out_name}.csv'

    try:
        with open(out_file, 'w') as f:
            f.write('tree_id\n')
            for tree_id in filtered_tree_list:
                f.write(f'{tree_id}\n')
    except IOError:
        general_utils.handle_errors(f'Could not write to the output file {out_file}', log_object, exit_stat=True)

    output_dir_name = None  # The copied files directory
    if copy_flag:
        output_dir_name = f'{output_dir}/{base_out_name}'
        general_utils.create_dir_if_needed(output_dir_name)
        for tree_id in filtered_tree_list:
            t = tree_dict[tree_id]
            tree_file_name = f'{output_dir_name}/{tree_id}.nw'
            try:
                t.write(format=1, outfile=tree_file_name)
            except IOError:
                general_utils.handle_errors(
                    f'Could not write the tree to the file {tree_file_name} in the directory {output_dir_name}',
                    exit_stat=True)
    return out_file, output_dir_name


def check_input_args(args):
    """
    Check the input arguments
    :param args: The input arguments
    :return: True if arguments are valid or false otherwise
    """
    valid = True
    if not (args.NOT or args.AND or args.OR or args.nodes or args.leaves):
        print('One of these arguments must be given: AND OR NOT nodes leaves')
        valid = False
    if args.nodes:
        if len(args.nodes) > 2:
            print(
                'The \'nodes\' argument must receive one or two arguments - '
                'for minimum size or for size range (including)')
            valid = False
        for size in args.nodes:
            if (not size.isdigit()) or (int(size) < 0):
                print('The \'nodes\' argument must receive one or two positive numbers')
                valid = False

    if args.leaves:
        if len(args.leaves) > 2:
            print(
                'The \'leaves\' argument must receive one or two arguments - '
                'for minimum size or for size range (including)')
            valid = False
        for size in args.leaves:
            if (not size.isdigit()) or (int(size) < 0):
                print('The \'leaves\' argument must receive one or two positive numbers')
                valid = False

    return valid


def get_pops_from_args(args):
    """
    Gets all the populations from the input arguments
    :param args: Input arguments
    :return: The population list
    """
    pop_list = []
    if args.OR:
        pop_list += args.OR
    if args.AND:
        pop_list += args.AND
    if args.NOT:
        pop_list += args.NOT

    pop_set = set([])
    for pop in pop_list:
        pop_set.add(pop)

    return pop_set


def search_by_size(tree_dict: dict, min_size: int, max_size: int):
    """
    Filters trees by the tree size.
    :param tree_dict: The tree dictionary of the format: tree_id: tree_size
    :param min_size: The minimum tree size to accept.
    :param max_size: The maximum tree size to accept.
    :return: The filtered tree list
    """
    filtered_trees = []
    for t_id in tree_dict:
        tree_size = tree_dict[t_id]
        if (tree_size >= min_size) and (tree_size <= max_size):
            filtered_trees.append(t_id)
    return filtered_trees


def filter_by_logic(tree_pops_dict: dict, or_pops, and_pops, not_pops):
    """
    Calls the filtering functions.
    :param tree_pops_dict:
    :param or_pops: The populations to filter by the 'or' rule
    :param and_pops: The populations to filter by the 'and' rule
    :param not_pops: The populations to filter by the 'not' rule
    :return: The three filtered lists
    """
    or_trees, and_trees, not_trees = [], [], []
    if or_pops:
        or_trees = search_by_or(tree_pops_dict, or_pops)
    if and_pops:
        and_trees = search_by_and(tree_pops_dict, and_pops)
    if not_pops:
        not_trees = search_by_not(tree_pops_dict, not_pops)
    return or_trees, and_trees, not_trees


def filer_by_size(tree_size_dict: dict, sizes: list):
    """
    Calls the filtering function.
    :param tree_size_dict: A dictionary in the format tree_id: size
    :param sizes: The input size argument
    :return: The filtered tree list.
    """
    size_trees = []
    if sizes:
        min_size = int(sizes[0])
        if len(sizes) < 2:  # Only the first argumnet was given
            max_size = 9999999999  # Arbitrary
        else:
            max_size = int(sizes[1])
        size_trees = search_by_size(tree_size_dict, min_size, max_size)
    return size_trees


def filter_rep(args):
    """
    The main function - filter the repetoire according to given filter parameters
    :param args: The input arguments
    :return: None
    """

    log_object = Logger(args.name, args.silent, 'filter')

    # Checking input sanity
    if check_input_args(args):

        # Collecting the trees
        tree_list = general_utils.collect_trees(args, log_object)

        log_object.info(f'Finish collecting {len(tree_list)} trees')

        pops = get_pops_from_args(args)

        args.ncors = int(args.ncors)
        number_of_workers = min(args.ncors, len(tree_list))
        tree_dict = general_utils.tree_list_to_tree_dict(tree_list)  # Creating a dictionary of tree_id: ETE tree
        tree_pops_dict, tree_node_dict, tree_leaves_dict = {}, {}, {}

        if args.OR or args.AND or args.NOT:
            with mp.Pool(number_of_workers) as pool:
                # Running poptree in parallel
                ret_list = pool.starmap(get_pops_of_one_tree, zip(tree_list, repeat(pops)))
            tree_pops_dict = attach_result_to_tree(tree_list, ret_list)
        if args.nodes:
            with mp.Pool(number_of_workers) as pool:
                ret_size_list = pool.map(general_utils.count_tree_nodes, tree_list)
            tree_node_dict = attach_result_to_tree(tree_list, ret_size_list)

        if args.leaves:
            with mp.Pool(number_of_workers) as pool:
                ret_size_list = pool.map(general_utils.count_tree_leaves, tree_list)
            tree_leaves_dict = attach_result_to_tree(tree_list, ret_size_list)

        log_object.info(f'Finish going over the trees')

        or_trees, and_trees, not_trees = filter_by_logic(tree_pops_dict, args.OR, args.AND, args.NOT)

        node_size_trees = filer_by_size(tree_node_dict, args.nodes)
        leaf_size_trees = filer_by_size(tree_leaves_dict, args.leaves)

        # Saving the files
        output_dir = general_utils.create_output_dir(args.name, 'Filtered_files', log_object)

        out_files, out_dirs = [], []
        if or_trees:
            out_file, out_dir = save_filtered_trees_to_file(or_trees, args.OR, 'or', tree_dict, output_dir, args.copy,
                                                            log_object)
            out_files.append(out_file)
            out_dirs.append(out_dir)
        if and_trees:
            out_file, out_dir = save_filtered_trees_to_file(and_trees, args.AND, 'and', tree_dict, output_dir,
                                                            args.copy, log_object)
            out_files.append(out_file)
            out_dirs.append(out_dir)
        if not_trees:
            out_file, out_dir = save_filtered_trees_to_file(not_trees, args.NOT, 'not', tree_dict, output_dir,
                                                            args.copy, log_object)
            out_files.append(out_file)
            out_dirs.append(out_dir)

        if node_size_trees:
            out_file, out_dir = save_filtered_trees_to_file(node_size_trees, args.nodes, 'nodes', tree_dict, output_dir,
                                                            args.copy, log_object)
            out_files.append(out_file)
            out_dirs.append(out_dir)

        if leaf_size_trees:
            out_file, out_dir = save_filtered_trees_to_file(leaf_size_trees, args.leaves, 'leaves', tree_dict, output_dir,
                                                            args.copy, log_object)
            out_files.append(out_file)
            out_dirs.append(out_dir)

        # Saying goodbye
        if out_files:
            log_object.info(f'The filtered trees are listed in: {" ".join(out_files)}')
            if args.copy:
                if out_files:
                    log_object.info(f'The filtered trees were copied to: {" ".join(out_dirs)}')
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
