#!/usr/bin/python3

"""
Converts Newick trees to dot for better drawing
"""

# Info
__author__ = 'Hadas Neuman'
__version__ = '1.5'
__date__ = '11/10/2020'

import multiprocessing as mp
import os
import re
import shutil
from argparse import ArgumentParser
from itertools import repeat  # For starmap

import utils.draw_utils as draw_utils
# Project imports
from utils import general_utils
from objects import Logger


# Global variables
# log_object = None


def adjust_line_length(name: str, line_len: int = 30):
    """
    Adjusts the line length by adding \n in the middle
    :param name: The node namae
    :param line_len: The maximun line length
    :return: The edited line
    """
    line_len = int(line_len)
    name = name.replace('-', '_')
    if len(name) > line_len:
        i = 0
        adj_name = ''
        while i < len(name):  # Go over the string
            x = name.rfind('_', i, i + line_len)  # Find the last '_'
            if x > 0:
                adj_name = adj_name + name[i:x] + '\\n'  # Add \n after the _
                i = x + 1
            else:
                adj_name = adj_name + name[i:i + line_len] + '\\n'  # Add \n in the end on the line
                i += line_len
    else:
        adj_name = name

    return adj_name


def new_to_dot(t, fontsize: int, line_width: int, output_file=None):
    """
    Converts the Newick to dot file
    :param fontsize: The font size in the dot file
    :param line_width: The line width in the node name
    :param t: An ETE tree object
    :param output_file: The output dot file name
    :return:
    """
    # print(t.get_ascii(show_internal=True))

    if t:
        tree_id = t.id
        dot_str = "digraph tree {\nedge [style=bold]\n"
        counter = 0
        root = t

        for node in t.traverse():  # Traversing by level order (default)
            if node == root:
                node.name = 'GL'
            else:
                if not hasattr(node, 'name'):  # An inferred node
                    node.name = ''
                elif 'Germline' in node.name: # Non root germline
                    continue # TEMP

            node.nid = str(counter)

            # Add node data
            node_label = adjust_line_length(node.name, line_width)

            node_str = f'{node.nid} [style=filled,fontcolor=black, fontsize={fontsize}, label="{node_label}"];'
            dot_str = dot_str + '\n' + node_str

            # Add the edge data for all nodes but the root
            if node != root:
                dist = round(node.get_distance(node.up))
                if dist != 1:
                    dist_label = f' [fontsize={int(int(fontsize) / 1.5)}, label="{dist}"]'
                else:
                    dist_label = ''
                edge_str = f'{node.up.nid} -> {node.nid}{dist_label};'
                dot_str = dot_str + '\n' + edge_str

            counter += 1
        dot_str = dot_str + '\n}'

        if not output_file:
            output_file = tree_id + '.dot'

        with open(output_file, 'w') as f:
            f.write(dot_str)
    else:
        output_file = None
    return output_file


def test_input(args):
    """
    Tests if the input is valid
    :param args: The input argument.
    :return: The argument error type. If no error was found - return None
    """
    arg_error = None
    valid_formats = ['bmp', 'canon', 'cmap', 'cmapx', 'cmapx_np', 'dot', 'eps', 'fig', 'gtk', 'gv', 'ico', 'imap',
                     'imap_np', 'ismap', 'jpe', 'jpeg', 'jpg', 'pdf', 'pic', 'plain', 'plain-ext', 'png', 'pov', 'ps',
                     'ps2', 'svg', 'svgz', 'tif', 'tiff', 'tk', 'vml', 'vmlz', 'x11', 'xdot', 'xlib']
    # Check if the 'dot' program exists
    if not shutil.which(args.dot):
        arg_error = f'The program {args.dot} does not exist.'

    # Check if the input files exist
    elif (general_utils.test_file_existence(args.tree)) and (not os.path.isdir(args.tree[0])):
        arg_error = f'The file {general_utils.test_file_existence(args.tree)} does not exist'

    # Check if the output format is supported by dot
    elif args.format not in valid_formats:
        arg_error = f'The format {args.format} is not supported by the dot program. ' \
                    f'The supported formats are {valid_formats}'

    # Check if the font size is a number
    elif not re.search(r'\d+', args.fontsize):
        arg_error = f'The fontsize parameter must except a number'

    return arg_error


def draw_using_dot(dot_file: str, dot_path: str, log_object, out_format='png'):
    """
    Draws the tree using the dot program
    :param dot_file: The dot file
    :param dot_path: The dot-program path
    :param log_object: The log object
    :param out_format: The output format
    :return: None
    """
    output_draw = f'{dot_file}.{out_format}'
    command = [dot_path, '-T', out_format, '-O', dot_file]
    try:
        general_utils.call_command(command, log_object)
    except OSError:
        error_str = f'Could not run the program {dot_path}'
        general_utils.handle_errors(error_str, log_object)
        output_draw = None
    return output_draw


def draw_trees(args):
    """
    The main function.
    :param args: The input arguments
    :return:
    """

    # global log_object

    log_object = Logger(args.name, args.silent, 'draw')
    # log_object.start()

    input_test_res = test_input(args)  # Check input arguments
    if input_test_res:  # There is a problem
        general_utils.handle_errors(input_test_res, log_object=log_object, exit_stat=True)  # Exit

    # Getting the input arguments
    tree_list = general_utils.collect_trees(args, log_object)

    intermediate_files = []

    log_object.info(f'Finish collecting {len(tree_list)} trees')

    # The current directory
    start_dir = os.getcwd()

    # Defining the output directory name
    output_dir = general_utils.create_output_dir(args.name, 'Drawn_trees', log_object)

    # Moving to the drawn trees directory
    os.chdir(output_dir)

    args.ncors = int(args.ncors)
    number_of_workers = min(args.ncors, len(tree_list))

    with mp.Pool(number_of_workers) as pool:

        # Converting all files to a dot format in parallel
        dot_files = pool.starmap(new_to_dot, zip(tree_list, repeat(args.fontsize), repeat(args.linewidth)))
        intermediate_files += dot_files

        if args.pops:
            color_dict = draw_utils.create_color_dict(args.pops)  # Get the color dictionary

            # Color all trees in parallel
            color_dot_trees = pool.starmap(draw_utils.color_trees, zip(dot_files, repeat(color_dict)))

            # Create a legend file
            legend_file = draw_utils.create_legend(color_dict)

            # Overwrite the dot_files variable
            dot_files = color_dot_trees + [legend_file]

            # For later deleting or moving
            intermediate_files += dot_files

        log_object.info(f'Finish going over the trees')

        # Drawing trees using the dot program in parallel
        pool.starmap(draw_using_dot,
                                   zip(dot_files, repeat(args.dot), repeat(log_object), repeat(args.format)))

        if not args.keep:
            # Delete the intermediate files
            general_utils.remove_files(intermediate_files)

    # Going back to the previous directory
    os.chdir(start_dir)

    # Print output files
    log_object.info(f'The output files were saved in: {output_dir}')
    log_object.done()


def get_arg_parser():
    """
    Defines the input and the output field help message
    :return: The parser
    """

    desc = "Counts observed mutations"

    # Define argument parser
    local_parser = ArgumentParser(description=desc)
    local_parser = general_utils.add_basic_arguments(local_parser)
    local_parser.add_argument('--dot', help='The dot program path. default: dot', default='dot')
    local_parser.add_argument('--format', help='The output format. default: png', default='png')
    local_parser.add_argument('--keep', action='store_true', help='Keep the intermediate dot-format files')
    local_parser.add_argument('-p', '--pops', action='store', nargs='+',
                              help='The population to search for in the trees')

    local_parser.add_argument('-s', '--fontsize', action='store', help='The font size in the output figure. '
                                                                       'Default - 40', default='40')
    local_parser.add_argument('--version', action='version', version='%(prog)s:' + ' %s:%s' % (__version__, __date__))

    return local_parser


if __name__ == "__main__":
    # Parse command line arguments
    parser = get_arg_parser()
    global_args = parser.parse_args()

    input_error = test_input(global_args)

    if input_error:
        print(input_error)
    else:
        draw_trees(global_args)
