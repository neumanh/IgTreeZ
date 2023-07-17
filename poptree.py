#!/usr/bin/python3

"""
Counts the number of transitions in repertoire's trees and saves the data on the population levels
"""

# Info
__author__ = 'Hadas Neuman'
__version__ = '1.6'
__date__ = '11/10/2020'

import datetime  # TEMP
import multiprocessing as mp
from argparse import ArgumentParser
from itertools import repeat  # For starmap

import numpy as np
import pandas as pd
from objects import Logger
from utils import general_utils, plot_utils, poptree_utils


def plot(trans_dict: dict, norm_dict_by_source: dict, norm_dict_by_dest: dict, pop_dict: dict, pop_count: np.array,
         name: str,
         outdir: str = '.', log_object: Logger = None):
    """
    Plots the poptree graphical outputs
    :param trans_dict: The transition dictionary
    :param norm_dict_by_source: The transition dictionary, normalized by source population
    :param norm_dict_by_dest: The transition dictionary, normalized by destination population
    :param pop_dict: The population dictionary
    :param pop_count: The number of populations per tree count
    :param name: The analysis name
    :param outdir: The output directory
    :param log_object: The Logger object
    :return: None
    """

    # Remove the prefix underline
    pop_dict = poptree_utils.replace_under_line_in_pops(pop_dict)
    trans_dict = poptree_utils.replace_under_line_in_trans(trans_dict)
    norm_dict_by_source = poptree_utils.replace_under_line_in_trans(norm_dict_by_source)
    norm_dict_by_dest = poptree_utils.replace_under_line_in_trans(norm_dict_by_dest)

    if not log_object:
        log_object = general_utils.create_general_log_object()

    if pop_dict:
        plot_utils.box_plot(pop_dict, fig_name=f'{outdir}/{name}_population_levels.pdf', xlab='Populations')
        plot_utils.bar_plot(pop_dict, fig_name=f'{outdir}/{name}_population_count.pdf')
        plot_utils.basic_pie_plot(pop_dict, fig_name=f'{outdir}/{name}_populations_pie.pdf')
    else:
        log_object.warning('No population detected for plotting')

    if trans_dict:
        plot_utils.box_plot(trans_dict, fig_name=f'{outdir}/{name}_transition_distance.pdf')
        plot_utils.bar_plot(trans_dict, fig_name=f'{outdir}/{name}_transition_count.pdf')
        plot_utils.sunburst_pie_plot(trans_dict, pop_dict, fig_name=f'{outdir}/{name}_transitions_pie.pdf')
    else:
        log_object.warning('No transition detected for plotting')

    # Plotting the normalized_by_source dictionary
    if norm_dict_by_source:
        plot_utils.bar_plot(norm_dict_by_source, xlab='Normalized transitions',
                            fig_name=f'{outdir}/{name}_transition_count_normalized_by_source.pdf', norm=True)
        plot_utils.sunburst_pie_plot(norm_dict_by_source, pop_dict,
                                     fig_name=f'{outdir}/{name}_transition_pie_normalized_by_source.pdf')

    # Plotting the normalized_by_destination dictionary
    if norm_dict_by_dest:
        plot_utils.bar_plot(norm_dict_by_dest, xlab='Normalized transitions',
                            fig_name=f'{outdir}/{name}_transition_count_normalized_by_destination.pdf', norm=True)
        plot_utils.sunburst_pie_plot(norm_dict_by_dest, pop_dict,
                                     fig_name=f'{outdir}/{name}_transition_pie_normalized_by_destination.pdf')

    if len(pop_count) > 0:
        plot_utils.plot_hist(pop_count, fig_name=f'{outdir}/{name}_population_count_frequency.pdf')


def save_output_tables(trans_dict: dict, norm_dict_by_source: dict, norm_dict_by_dest: dict, pop_dict: dict, trans_df,
                       trans_df_dist, pop_df, name: str, out_dir: str):
    """
    Saves the output tables
    :param trans_dict: The transitions dictionary
    :param norm_dict_by_source: The transitions dictionary, normalized by source population
    :param norm_dict_by_dest: The transitions dictionary, normalized by destination population
    :param pop_dict: The populations dictionary
    :param trans_df: The transitions dataframe (with tree id)
    :param trans_df_dist: The transition distances dataframe (with tree id)
    :param pop_df: The populations dataframe (with tree id)
    :param name: The analysis name
    :param out_dir: The output directory
    :return: None
    """
    # Defining output names
    transitions_by_tree_name = f'{out_dir}/{name}_transition_counts_by_tree.csv'
    populations_by_tree_name = f'{out_dir}/{name}_population_counts_by_tree.csv'
    transitions_by_tree_name_with_dist = f'{out_dir}/{name}_transition_distances_by_tree.csv'
    transitions_csv_name = f'{out_dir}/{name}_transition_distances.csv'
    transitions_norm_by_source_csv_name = f'{out_dir}/{name}_transition_counts_normalized_by_source.csv'
    transitions_norm_by_dest_csv_name = f'{out_dir}/{name}_transition_counts_normalized_by_destination.csv'
    populations_csv_name = f'{out_dir}/{name}_population_levels.csv'
    transitions_stats_name = f'{out_dir}/{name}_transition_distances_summary.csv'
    transitions_nor_source_stats_name = f'{out_dir}/{name}_transition_summary_normalized_by_source.csv'
    transitions_nor_dest_stats_name = f'{out_dir}/{name}_transition_summary_normalized_by_destination.csv'
    populations_stats_name = f'{out_dir}/{name}_population_levels_summary.csv'

    # Saving the outputs
    if len(trans_dict) > 0:
        save_to_file(trans_dict, transitions_csv_name, transitions_stats_name, name)
    if len(norm_dict_by_source) > 0:
        save_to_file(norm_dict_by_source, transitions_norm_by_source_csv_name, transitions_nor_source_stats_name, name)
    if len(norm_dict_by_dest) > 0:
        save_to_file(norm_dict_by_dest, transitions_norm_by_dest_csv_name, transitions_nor_dest_stats_name, name)
    if len(pop_dict) > 0:
        save_to_file(pop_dict, populations_csv_name, populations_stats_name, name)
    if trans_df.any:
        # Insert a sample name in the beginning of the name
        trans_df = poptree_utils.add_sample_column(trans_df, name)
        trans_df.to_csv(transitions_by_tree_name, index=False)
    if pop_df.any:
        # Insert a sample name in the beginning of the name
        pop_df = poptree_utils.add_sample_column(pop_df, name)
        pop_df.to_csv(populations_by_tree_name, index=False)
    if trans_df_dist.any:
        # Insert a sample name in the beginning of the name
        trans_df_dist = poptree_utils.add_sample_column(trans_df_dist, name)
        trans_df_dist.to_csv(transitions_by_tree_name_with_dist, index=False)


def poptree(args):
    """
    The main function.
    :param args: The input arguments
    :return:
    """
    log_object = Logger(args.name, args.silent, 'poptree')

    # Getting the input arguments
    pops = args.pops
    log_object.info(f'Looking for the populations:{pops}')

    tree_list = general_utils.collect_trees(args, log_object)

    log_object.info(f'Finish collecting {len(tree_list)} trees')

    args.ncors = int(args.ncors)
    number_of_workers = min(args.ncors, len(tree_list))
    with mp.Pool(number_of_workers) as pool:
        # Running poptree in parallel
        result = pool.starmap(create_dicts, zip(tree_list, repeat(pops), repeat(args.zero), repeat(args.max_dist)))
    log_object.info(f'Finish going over the trees')

    # Collecting the output
    trans_dict, pop_dict, trans_list, pop_list, trans_list_dist = unite_starmap_poptree_results(result)

    norm_dict_by_source = poptree_utils.normalize_by_the_source(trans_dict, pop_dict)
    norm_dict_by_dest = poptree_utils.normalize_by_the_dest(trans_dict, pop_dict)

    # Sorting and updating the dictionaries
    pop_dict = poptree_utils.sort_pops(pop_dict, args.pops)
    trans_dict = poptree_utils.sort_trans(trans_dict, args.pops)
    norm_dict_by_source = poptree_utils.sort_trans(norm_dict_by_source, args.pops)
    norm_dict_by_dest = poptree_utils.sort_trans(norm_dict_by_dest, args.pops)

    output_dir = general_utils.create_output_dir(args.name, 'PopTree_results', log_object)

    trans_df = pd.DataFrame(trans_list)
    trans_df_dist = pd.DataFrame(trans_list_dist)
    pop_df = pd.DataFrame(pop_list)

    pop_count = np.array(poptree_utils.count_pops_num_in_tree(pop_df))

    if args.plot:
        sub_output_dir = f'{output_dir}/plots'
        general_utils.create_dir_if_needed(sub_output_dir)
        plot(trans_dict, norm_dict_by_source, norm_dict_by_dest, pop_dict, pop_count, args.name, outdir=sub_output_dir,
             log_object=log_object)

    # From the (A, B) form to the A_to_B form
    trans_dict = edit_transitions_dict(trans_dict)
    norm_dict_by_source = edit_transitions_dict(norm_dict_by_source)
    norm_dict_by_dest = edit_transitions_dict(norm_dict_by_dest)

    # Saving the analysis results to files
    save_output_tables(trans_dict, norm_dict_by_source, norm_dict_by_dest, pop_dict, trans_df, trans_df_dist, pop_df,
                       args.name, output_dir)

    log_object.info(f'The output files were saved in: {output_dir}')
    log_object.done()


def zero_count_in_node(node_name: str, pops: list):
    """
    Tests weather there are two populaition in the tree node.
    :param node_name: The node name
    :param pops: The populations
    :return: A list of tuples: [(A,B) (B,C).. ]
    """
    ret_list = []
    for idx1, pop1 in enumerate(pops):
        for idx2, pop2 in enumerate(pops):

            if idx1 < idx2:  # Also prevent pop1 = pop2
                if (pop1 in node_name) and (pop2 in node_name):
                    # Add the populations as a tuple
                    ret_list.append((pop1, pop2))
    return ret_list


def add_dist_to_dict(key, dist: float, dic: dict):
    """
    Adds an item to the dictionary. If the key exists in the dictionary -
    adds the value to the dictionary to the list of the key
    :param key: the key of the dictionary
    :param dist: the distance to add
    :param dic: the dictionary to add to
    :return: the updated dictionary
    """

    if key not in dic:  # The tuple is the key in the transition count dictionary
        dic[key] = np.array([dist])
    else:
        dic[key] = np.append(dic[key], np.array([dist]))

    return dic


def unite_two_dicts(dict1: dict, dict2: dict):
    """
    Adds the items to a dictionary to the other. If the key exists in the dictionary -
    adds the value to the dictionary to the list of the other dictionary key's
    :param dict1: The first dictionary
    :param dict2: The second dictionary
    :return: A united dictionary
    """
    ret_dict = dict1
    for key2 in dict2:
        if key2 not in ret_dict:
            ret_dict[key2] = dict2[key2]
        else:
            ret_dict[key2] = np.append(ret_dict[key2], dict2[key2])

    return ret_dict


def sum_dict(dic: dict):
    """
    Sums the number of elements in each dictionary item (a list)
    :param dic: A dictionary of lists
    :return: A summed-up dictionary
    """
    sum_dic = {}
    for key in dic:
        sum_dic[key] = len(dic[key])
    return sum_dic


def str_dict(dic: dict):
    """
    Returns a string of elements in each dictionary item (a list) in the format '1 2 2 1 3'
    :param dic: A dictionary of lists
    :return: A summed-up dictionary
    """
    sum_dic = {}
    for key in dic:
        temp_dist_list = []
        for dist in dic[key]:
            temp_dist_list.append(str(int(dist)))
        sum_dic[key] = ' '.join(temp_dist_list)
    return sum_dic


def unite_starmap_poptree_results(res: list):
    """
    Unites the poptree-starmap results
    :param res: A list of tuples in the format - transition-dictionary, population-dictionary
    :return: Two dictionaries - the transition dictionary and the population dictionary
    """
    transitions_dict = {}
    pops_dict = {}
    transitions_list = []
    pops_list = []
    transitions_list_dist = []

    for three_dicts in res:
        # Update the dictionaries for plotting
        curr_transitions_dict = three_dicts[0]
        curr_pops_dict = three_dicts[1]
        curr_transitions_id = three_dicts[2]
        # Unite the dictionaries
        transitions_dict = unite_two_dicts(transitions_dict, curr_transitions_dict)
        pops_dict = unite_two_dicts(pops_dict, curr_pops_dict)

        id_item = {general_utils.id_str: curr_transitions_id}

        # Add to the lists of dictionaries
        edited_curr_transitions_dict = edit_transitions_dict(curr_transitions_dict)
        transitions_list.append({**id_item, **sum_dict(edited_curr_transitions_dict)})
        transitions_list_dist.append({**id_item, **str_dict(edited_curr_transitions_dict)})
        pops_list.append({**id_item, **sum_dict(curr_pops_dict)})

    return transitions_dict, pops_dict, transitions_list, pops_list, transitions_list_dist


def edit_transitions_dict(trans_dict: dict):
    """
    Replace the keys from the (A, B) form to A_to_B
    :param trans_dict: The transitions dictionary
    :return: The edited transitions dictionary
    """
    edited_dict = {}
    for key in trans_dict:
        new_key = f'{key[0]}_to_{key[1]}'
        edited_dict[new_key] = trans_dict[key]
    return edited_dict


def save_to_file(data_dict: dict, table_file_name: str, stat_file_name: str, sample_name: str):
    """
    Saves poptree's results
    :param data_dict: The transition or population dictionary
    :param table_file_name: The output data file name
    :param stat_file_name: The output summary file name
    :param sample_name: The analysis name
    :return: None
    """
    df = pd.DataFrame.from_dict(data_dict, orient='index')
    df = df.transpose()

    # Insert a sample name in the beginning of the dataframe
    df = poptree_utils.add_sample_column(df, sample_name)

    # Saving to CSV
    df.to_csv(table_file_name, index=False, index_label=False)

    # Summarize the results and saves to CSV
    df_sum = df.describe()
    df_sum.to_csv(stat_file_name)


def create_dicts(t, pops: list, zero: bool = False, max_dist = None):
    """
    Creates the transition and population dictionaries
    :param t: An ETE tree
    :param pops: Populations to search for in the tree
    :param zero: Whether to count transitions in the same node
    :return: transitions dictionary, populations dictionary, and tree-id
    """
    trans_dict, pop_dict = {}, {}
    tree_id = None

    # print(t.get_ascii(show_internal=True))
    if t:
        tree_id = t.id

        # Getting the path list
        path_list = poptree_utils.get_path_list(t)
        trans_set = set([])  # For the unique recognition based on node's name
        pop_set = set([])  # For the unique recognition based on node's name
        node_set = set([])

        for path in path_list:
            # Initializing the variables
            node1, node2, node0 = None, None, None
            pop1, pop2, pop0 = None, None, None
            idx1, idx2, idx0 = 0, 0, 0
            dist_idx = 0

            # Going over all the paths
            for idx, node_object in enumerate(path):
                # print('########### idx', idx )
                dist_idx = dist_idx + node_object.dist
                node_name = node_object.name

                for pop in pops:
                    if pop in node_name:

                        # Saving the nodes' names
                        if not poptree_utils.multiple_pops_in_name(node_name, pops):  # To prevent unknown transition
                            if not node2:
                                node2 = node_name
                                pop2 = pop
                                idx2 = idx
                                # print('pop2', pop2, 'idx2', idx2)

                            elif not node1:
                                node1 = node_name
                                pop1 = pop
                                idx1 = idx

                        node0 = node_name
                        idx0 = idx

                # Saving the transition data
                if pop1 == pop2:
                    node2, pop2 = node1, pop1  # Move to the next transition
                    node1, pop1 = None, None

                elif node1 and node2:  # And not pop1 == pop2
                    dist = poptree_utils.sum_path(path, idx2) - poptree_utils.sum_path(path, idx1)
                    if ((node1, node2) not in trans_set) and ((max_dist is None) or (dist <= max_dist)):  # A new transition
                        trans_set.add(
                            (node1, node2))  # Add the tuple to the set to make share averey edge is handled once
                        # Update the transition dictionary,
                        # while the tuple is the key in the transition count dictionary
                        trans_dict = add_dist_to_dict(key=(pop1, pop2), dist=dist, dic=trans_dict)

                    node2, pop2, idx2 = node1, pop1, idx1  # Move to the next transition
                    node1, pop1 = None, None

                # Count zero-length transitions
                if (zero and node0) and (node0 not in node_set):
                    node_set.add(node0)  # Each node is handled once
                    zero_trans = zero_count_in_node(node0, pops)
                    for zt in zero_trans:
                        trans_dict = add_dist_to_dict(key=zt, dist=0, dic=trans_dict)

                # Saving the population statistics data
                if node0:
                    # Looking for all the populations in the node's name
                    for pop in pops:
                        if pop in node0:
                            pop0 = pop
                            node_tup = (node0, pop)  # Validating that each node+pop will considered only once
                            if node_tup not in pop_set:  # A new node
                                # dist = len(path) - idx0 - 1
                                dist = poptree_utils.sum_path(path, idx0)
                                pop_set.add(node_tup)  # Add to the set
                                # Update the population dictionary,
                                # while the population is the key in the transition count dictionary
                                pop_dict = add_dist_to_dict(key=pop0, dist=dist, dic=pop_dict)

    return trans_dict, pop_dict, tree_id


def get_arg_parser():
    """
    Defines the input and the output field help message
    :return: The parser
    """

    desc = "Counts observed mutations"

    # Define argument parser
    local_parser = ArgumentParser(description=desc)
    local_parser = general_utils.add_basic_arguments(local_parser)
    local_parser.add_argument('-p', '--pops', action='store', nargs='+',
                              help='The population to search for in the trees',
                              required=True)

    local_parser.add_argument('-z', '--zero', action='store_true', help='Analyse zero distance transitions as well')
    local_parser.add_argument('--version', action='version', version='%(prog)s:' + ' %s:%s' % (__version__, __date__))

    return local_parser


if __name__ == "__main__":
    # Parse command line arguments
    parser = get_arg_parser()
    global_args = parser.parse_args()

    poptree(global_args)
    print("Done:", datetime.datetime.now())
