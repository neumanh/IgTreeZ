#!/usr/bin/python3

"""
Utility functions for drawing
"""

# Info
__author__ = 'Hadas Neuman'
__version__ = '1.1'
__date__ = '22/9/2020'

import os
import re
import subprocess
from os.path import basename

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

# Global variables
keep_silence = False
delete_all_lables = False


def multiple_pops(node_name: str, pops: dict):
    """
    Find if there are more than one population in node's name
    :param node_name: The node name
    :param pops: The population names
    :return: True if there are more than one population in node's name
    """
    first_pop = False
    second_pop = False
    for pop in pops:
        if pop in node_name:
            if first_pop and (not second_pop):
                second_pop = True
            else:
                first_pop = True

    return second_pop


def from_two_list_to_dict(keys: list, values: list):
    """
    Creates a dictionary from two list, while the first list represents the keys, and the second list - the values.
    :param keys: The keys list
    :param values: The values list
    :return:
    """
    dic = {}
    for idx, key in enumerate(keys):
        dic[key] = values[idx]
    return dic


def color_multiple_pops(color_dic: dict, line: str):
    """
    Colors the node's by population
    :param color_dic: The population-color dictionary
    :param line: The line in the dot file to update
    :return: The updated line
    """
    # color_dic = create_color_dict(pops, colors) # Convert to dictionary
    label = find_label_in_node(line)  # Find the label

    pops_in_label = []

    # Collect the relevant populations
    for pop in color_dic:
        if pop in label:
            pops_in_label.append(pop)  # Add to the list

    # Create the color string
    colors_string = ""
    for pop in pops_in_label:
        colors_string = colors_string + color_dic[pop] + ":"

    # Removing the last ":"
    regex = r':$'
    m = re.match(regex, colors_string)
    if m:
        colors_string = re.sub(regex, '', colors_string)

    full_color_string = f', style = \"wedged\", fillcolor = \"{colors_string}\"]'
    updated_line = line.replace(']', full_color_string)
    return updated_line


def find_label_in_node(node_name: str):
    """
    Finds the label in the node's name
    :param node_name: The node's name
    :return: The label
    """
    label = -1

    regex = r'^\w+\s+\[.*label\s*=\s*\"(.*?)\"'  # A node description

    m = re.match(regex, node_name)

    if m:
        label = m.group(1)
    return label


def color_trees(file: str, color_dict: dict, out_dir: str = "."):
    """
    Colors a given dot tree by a given pattern and deletes all labels
    :param color_dict: The population-colors dictionary
    :param out_dir: The output directory
    :param file: The dot file
    :return: The colored dot file
    """
    colored_tree_file = None

    if file:
        file_ext = os.path.splitext(file)[1]
        temp_file_name = basename(file)
        colored_tree_file = f'{out_dir}/{temp_file_name.replace(file_ext, "_colored" + file_ext)}'
        with open(file, "rt") as fin:  # The non colored file
            with open(colored_tree_file, "wt") as fout:  # The output, colored file
                for line in fin:
                    line_to_write = line
                    # for i in range(0, len(populations)):
                    for pop in color_dict:
                        color = color_dict[pop]
                        node_label = find_label_in_node(line)
                        if node_label != -1:  # The is a label
                            # If the label was found to be '' or contains I/inferred- a hypothetical node
                            if (node_label == '') or ('nferred' in node_label) or re.match(r'^\d{1,3}$', node_label):
                                if 'fillcolor=white' not in line_to_write:  # The color was not added previously
                                    line_to_write = line_to_write.replace(']', f', fillcolor=white]')
                            elif pop in node_label:
                                if multiple_pops(node_label, color_dict):
                                    line_to_write = color_multiple_pops(color_dict, line)
                                else:
                                    # Color the node
                                    line_to_write = line_to_write.replace(']', f', fillcolor=\"{color}\"]')

                    # Delete the labels of the nodes
                    regexp = re.compile(r'^[0-9]+ \[')  # a node
                    if (regexp.search(line_to_write)) or delete_all_lables:  # This is a node and not an edge
                        line_to_write = re.sub(r'label=\".*?\"', 'label=\"\"', line_to_write)
                    fout.write(line_to_write)

    return colored_tree_file


def create_legend(color_dict: dict, legend_file_name: str = "legend.dot", directory: str = "."):
    """
    Creates a legend base on populations and colors
    :param color_dict: The population-color dictionary
    :param directory: The output directory
    :param legend_file_name: The chosen output legend name
    :return: legend_file_name
    """
    legend_file_name = f'{directory}/{legend_file_name}'
    hd = "\"Hypothetical\\nsplit node\""
    with open(legend_file_name, "w") as fout:
        # Writing the header and the hypothetical node
        fout.write('digraph legend {\nsubgraph cluster_legend {\nlabel = \"Legend\"\ncolor=black\n\n' + hd + '[]')
        arrows_string = ""
        for pop0 in color_dict:
            # pop = pop0.replace("-", "_")  # For later drawing by dot
            pop = re.sub(r'^_', '', pop0)  # Removing the _ in the beginning
            fout.write(f'\"{pop}\" [fillcolor=\"{color_dict[pop0]}\" style=filled];\n')

            arrows_string = arrows_string + "\"" + pop + "\""
            arrows_string = arrows_string + " -> "
        arrows_string = arrows_string + hd
        arrows_string = arrows_string + ' [style=\"invis\"]\n'
        fout.write(arrows_string)

        fout.write("}\n}")  # Closing the graph

        return legend_file_name


def draw_tree(tree_file: str, out_format: str):
    """
    Draws tree using the dot program, and saves them by the given output format
    :param out_format: The dot output format
    :param tree_file:
    :return: The drawn tree name
    """
    format_str = f'-T{out_format}'  # the output format string
    command = ["dot", format_str, tree_file, "-O"]  # The dot command
    subprocess.run(command)  # rRunning the dot command
    out_file = f'{tree_file}.{out_format}'  # The output file
    if not os.path.exists(out_file):  # If the file was created - returned the file name. if not - return None
        out_file = None
        print("Error: no figure was created by dot. Command:\n", " ".join(command))
    return out_file


def initiate_vars(args):
    """
    Initiate variables.
    :param args: The input arguments
    :return: None
    """
    global keep_silence, delete_all_lables

    if args.dir:
        directory = args.dir
        if not (os.path.exists(directory)):  # Create a dir if doesn't exist
            os.makedirs(directory)

    if args.silent:
        keep_silence = True

    if args.delete_all:
        delete_all_lables = True

    populations = args.populations

    # Input colors
    if args.colors:
        colors = args.colors
        if len(colors) < len(populations):
            print("Error: the number of colors must be equal or greater than the number of populations")
            exit()

        color_dic = from_two_list_to_dict(populations, colors)
    else:
        color_dic = create_color_dict(populations)
    return color_dic


def create_color_dict(pops):
    """
    Creates a dictionary of "population: color"
    :param pops: The populations
    :return: The color dictionary
    """
    delta = 1 / len(pops)
    col = 0.0 + (delta / 2)
    color_palt = plt.cm.Spectral
    color_dict = {}

    # Initial value
    hex_col = mcolors.to_hex(color_palt(col))

    for pop in pops:
        color_dict[pop] = hex_col
        col = col + delta
        hex_col = mcolors.to_hex(color_palt(col))
    return color_dict


def create_color_dict_from_colors(pops: list, colors: list):
    """
    Creates a dictionary of "population: color"
    :param pops: The populations
    :param colors: The colors
    :return: The color dictionary
    """
    color_dict = {}

    i = 0
    for pop in pops:
        color_dict[pop] = colors[i]
        i += 1
    return color_dict
