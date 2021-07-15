#!/usr/bin/python3
"""
Counts the mutations in Newick trees.
"""

# Info
__author__ = 'Hadas Neuman'
__version__ = '1.0'
__date__ = '3/8/20'

import glob
import multiprocessing as mp
import os
import re
import shutil
import subprocess
import sys

from ete3 import Tree
from objects import Logger

# Global variables
sample_str = 'sample'
id_str = 'tree_id'
parent_dir_name = 'IgTreeZ_output'
# defined by IMGT http://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html
imgt_regions = {'fwr1': 0, 'cdr1': 78, 'fwr2': 114, 'cdr2': 165, 'fwr3': 195, 'cdr3': 313, 'cdr3_end': 312}

curr_log_object = None


def get_node_names(t: Tree):
    """
    Gets the names of the nodes
    :param t: An ETE tree
    :return: A list of node names
    """

    names = []

    for node in t.traverse():
        if hasattr(node, 'name'):
            names.append(node.name)
    return names


def is_equal_len(seqs):
    """
    Test if whether all sequences have the same length
    :param seqs: The sequences to compare
    :return: True if the length are the same, and False otherwise
    """
    same_len = True
    if len(seqs) > 1:
        seq_len = len(seqs[0])
        for seq in seqs:
            if len(seq) != seq_len:
                same_len = False
    return same_len


def handle_errors(error_str: str, log_object: Logger = None, exit_stat=False):
    """
    Handles errors.
    :param error_str: The error explanation
    :param log_object: The log object
    :param exit_stat: If True - exit program
    :return: None
    """
    if log_object:
        log_object.error(error_str)
    else:
        print(f'Error: {error_str}')
    if exit_stat:
        sys.exit()


def replace_new_line_in_ig_tree(t: Tree):
    """
    Removes everything after the \n character is sequence names
    :param t: A tree
    :return: The updated tree
    """
    if t:
        for node in t.traverse("postorder"):  # Going from leaves to root to generate missing sequences
            if hasattr(node, 'name'):
                seq_name = node.name
                if '\\n' in seq_name:
                    seq_name = re.sub(r'\\n.*', '', seq_name)
                    node.name = seq_name
    return t


def get_clone_name(filename):
    """
    Gets the basic name of a given clone
    :param filename: The filename
    :return: The clone name
    """
    basname = os.path.basename(filename)  # Deleting the path
    basename_no_ext = os.path.splitext(basname)[0]  # Deleting the extension
    return basename_no_ext


def get_clone_number(clone_name: str):
    """
    Finds the clone number from the clone name
    :param clone_name: The clone name
    :return: The number associated with this clone
    """
    number = clone_name
    number_regex = r'(\d+)'
    m = re.search(number_regex, clone_name)
    if m:
        number = int(m.group(1))

    return number


def get_files_from_dir(tree_dir: str):
    """
    Gets tree files list - for processing large number of input files
    :param tree_dir: The directory of the Newick trees
    :return: The sequences file list
    """
    tree_list = []
    for filename in glob.glob(tree_dir + '/*.*'):
        tree_list.append(filename)
    return tree_list


def get_tree_files(tree_dir: str):
    """
    Gets tree files list - for processing large number of input files
    :param tree_dir: The directory of the Newick trees
    :return: The sequences file list
    """
    tree_list = []
    for filename in glob.glob(tree_dir + '/*nw'):
        tree_list.append(filename)
    return tree_list


def test_file_existence(file_list):
    """
    Tests file existence
    :param file_list: Input files to test
    :return: The name of the file which does not exist
    """
    not_exs = None
    for file in file_list:
        if (not os.path.isfile(file)) and (not not_exs):
            not_exs = file
    return not_exs


def print_list(list_to_print: list):
    """
    Prints list, usually of output files
    :param list_to_print: The list to print
    :return: None
    """
    for item in list_to_print:
        print(item)


def create_dir_if_needed(dir_name):
    """
    Check if the directory exists, and create it if not.
    :param dir_name: The directory name
    :return: None
    """
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)


def move_files_to_dir(files: list, dir_name: str, log_object: Logger = curr_log_object):
    """
    Moves all given files to a given directory. If the file does not exist - handles the error
    :param files: Files to move
    :param dir_name: The directory name
    :param log_object: The Logger object
    :return: None
    """
    if not log_object:
        log_object = create_general_log_object()

    # Creates the output directory only if needed
    create_dir_if_needed(dir_name)

    for file in files:
        if file:
            if not os.path.isfile(file):
                handle_errors(f'The file {file} was not created', log_object)
            else:
                shutil.move(file, dir_name + '/' + file)


def remove_files(files: list):
    """
    Removes files. Check if the file exists to prevent an error.
    :param files: List of files to delete.
    :return: The names of the missing files.
    """

    missing_files = []

    for file in files:
        if file:
            if os.path.isfile(file):
                os.remove(file)
            else:
                missing_files.append(file)

    return missing_files


def sort_dict(dic: dict):
    """
    Sorts a dictionary
    :param dic: A dictionary to sort
    :return: The sorted dictionary
    """
    sorted_dict = {}
    for key in sorted(dic):
        sorted_dict[key] = dic[key]
    return sorted_dict


def call_command(command: list, log_object: Logger = curr_log_object, out_file: str = None):
    """
    Calls OS command using subprocess.Popopen
    :param command: The command to run
    :param log_object: The logger object
    :param out_file: The output file to write to
    :return: None
    """
    if not log_object:
        log_object = create_general_log_object()

    e = 0  # The error code to return

    log_object.info(f'Running: {" ".join(command)}', command=True)

    try:
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()

        stdout_str = stdout.decode("utf-8")
        stderr_str = stderr.decode("utf-8")

        if stdout:
            if out_file:
                with open(out_file, "w") as f:
                    f.write(stdout_str + "\n")
            else:
                log_object.info(stdout_str)

        if stderr:
            log_object.warning(stderr_str)

            print(stderr_str)
            print('An error accrued. See the log file for details.\n\n')

    except OSError as e:
        error_str = f'Could not run the command:{command}: {e}\n'
        log_object.error(error_str)
    except Exception as e:
        error_str = f'Could not run the command:{command}: {e}\n'
        log_object.error(error_str)
    return e


def create_tree(tree_file: str, log_object: Logger = curr_log_object, frmt=1):
    """
    Creates an ETE tree from a given tree file.
    :param tree_file: Newick format tree file.
    :param log_object: The Logger object
    :param frmt: The ETE format
    :return: The ETE tree.
    """
    if not log_object:
        log_object = create_general_log_object()

    try:
        t = Tree(tree_file, format=frmt)
    except Exception:
        t = None
        err_str = f'Could not open the file {tree_file} as a Newick tree'
        handle_errors(err_str, log_object)
    return t


def create_general_log_object():
    """
    Creates a general log object.
    :return: The genreal log object.
    """
    log_object = Logger('IgTreeZ', False)
    return log_object


def create_output_dir(name: str, sub_dir_name=None, log_object: Logger = None):
    """
    Creates an output directory and the initial log file
    :param name: The analysis name
    :param sub_dir_name: The sub-directory
    :param log_object: The log object
    :return: The directory name
    """
    global curr_log_object

    create_dir_if_needed(parent_dir_name)
    dir_name = f'{parent_dir_name}/{name}'
    create_dir_if_needed(dir_name)

    if log_object:
        curr_log_object = log_object
    else:
        curr_log_object = create_general_log_object()

    curr_log_object.move_log(dir_name)

    if sub_dir_name:
        dir_name = f'{parent_dir_name}/{name}/{sub_dir_name}'
        create_dir_if_needed(dir_name)

    return dir_name


def collect_trees(args, log_object: Logger = curr_log_object):
    """
    Collects the trees according to the input arguments.
    :param args: The input arguments
    :param log_object: The log object
    :return: The tree list
    """
    if not log_object:
        log_object = create_general_log_object()

    tree_list = []
    if args.json:
        from utils import json_utils
        tree_list = json_utils.collect_trees_from_scheme(args.json, log_object)
    elif args.tree:
        tree_file_list = []
        for item in args.tree:
            if os.path.isdir(item):
                tree_file_list += get_files_from_dir(item)
            else:  # A file
                tree_file_list.append(item)

        # Check if the input files exist
        if test_file_existence(tree_file_list):
            error_str = f'The file {test_file_existence(tree_file_list)} does not exist. Exiting.'
            handle_errors(error_str, log_object, exit_stat=True)
        else:  # All files exist
            for tree_file in tree_file_list:
                t = create_tree(tree_file, log_object)
                # t.id = get_clone_number(get_clone_name(tree_file))
                if t:
                    t.id = get_clone_name(tree_file)
                    t.file = os.path.abspath(tree_file)
                    tree_list.append(t)

    # Sample the tree if needed
    if args.sample:
        tree_list = sample_tree_list(tree_list, args.sample)

    if len(tree_list) < 1:
        handle_errors('Could not find trees', log_object, exit_stat=True)

    return tree_list


def is_real_node(node: Tree):
    """
    Checks if the node is a real node. Based on IgTree, GCtree and ChangeO inference methods.
    :param node: A tree node
    :return: True if this is a real node.
    """
    true_node = True
    if hasattr(node, 'name'):
        if ('nferred' in node.name) or ('Germline' == node.name) or (re.match(r'^\d{1,3}$', node.name)):
            true_node = False
    else:
        true_node = False
    return true_node


def count_tree_nodes(t: Tree):
    """
    Counts the number of real nodes in tree
    :param t: An ETE tree object
    :return: The number tree nodes with names
    """
    count = -1  # Ignoring the root
    for _ in t.traverse():
        # if is_real_node(node):
        count += 1
    return count


def count_tree_leaves(t: Tree):
    """
    Counts the number of leavess in tree
    :param t: An ETE tree object
    :return: The number tree nodes with names
    """
    leaves = t.get_leaves()
    leaves_num = 0
    for leaf in leaves:
        if 'Germline' not in leaf.name:
            leaves_num += 1

    return leaves_num


def tree_list_to_tree_dict(tree_list: list):
    """
    Creates a dictionary in the format tree_id: tree
    :param tree_list: The tree list
    :return: The tree dictionary
    """
    tree_dict = {}
    for t in tree_list:
        tree_dict[t.id] = t
    return tree_dict


def sample_tree_list(population: list, k: str, log_object: Logger = curr_log_object):
    """
    Samples the input tree list by k
    :param population: The tree list
    :param k: The sample size
    :param log_object: The Logger object
    :return: The samples tree list
    """
    import random
    sampled_list = []

    if not log_object:
        log_object = create_general_log_object()

    # TODO: move the check to input validation
    if k.isnumeric():
        k = int(k)
    else:
        handle_errors(f'Sample size {k} is not a valid number', log_object, exit_stat=True)

    if k <= 0:
        handle_errors(f'Sample size should be greater than zero', log_object, exit_stat=True)
    elif k > len(population):
        log_object.warning(
            f'Sample size {k} is bigger than the list size {len(population)}. Taking the whole tree list.')
        sampled_list = population
    else:
        sampled_list = random.sample(population, k)

    return sampled_list


def add_basic_arguments(parser):
    """
    Adds some common basic arguments to an ArgumentParser object
    :param parser: ArgumentParser object
    :return: The updated ArgumentParser object
    """
    accessible_cores = mp.cpu_count()  # IgBlast uses 5 parallel processes

    # Define argument parser
    parser.add_argument('-n', '--name', action='store', help='The analysis name', required=True)

    parser.add_argument('-t', '--tree', action='store', nargs='+', help='The trees in Newick format',
                        required=True)
    parser.add_argument('-j', '--json', action='store',
                        help='A jason file in the AIRR \'Clone and Lineage Tree Schema\' format,'
                             'as described in '
                             'https://docs.airr-community.org/en/stable/datarep/clone.html')

    parser.add_argument('--ncors', help='The number of CPU cores. default ' + str(accessible_cores),
                        default=accessible_cores)

    parser.add_argument('--sample', help='Samples the input tree list by the given number')

    parser.add_argument('--silent', action='store_true', help='Avoid output')

    return parser
