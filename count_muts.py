#!/usr/bin/python3
"""
Counts the mutations in Newick trees.

version 1.6 - uses generators for memory saving
"""

# Info
__author__ = 'Hadas Neuman'
__version__ = '1.6'
__date__ = '1/10/20'

# Imports
import multiprocessing as mp
from argparse import ArgumentParser
from itertools import repeat  # For starmap

import pandas as pd
from objects import Logger
from utils import fasta_utils, general_utils, mut_utils, db_utils, tree_utils

# Global variables
SAMPLE_FIELD = ""

gls = ["GL", "G.L", 'G.L.', "Germline", "Germline_UCA", 'EX;G.L', 'EX_G.L']


def order_mut_df_columns(df: pd.DataFrame, log_object: Logger, sample_name: str):
    """
    Orders the mutations dataframe columns.
    :param df: The mutations dataframe.
    :param log_object: The logger object.
    :param sample_name: The sample name.
    :return: The general ordered dataframe and the selection ordered dataframe.
    """
    new_columns = []
    new_shazam_columns = []

    # Inserting the sample column in the beginning of the dataframe
    df[general_utils.sample_str] = sample_name

    all_shazam_columns = [general_utils.sample_str, general_utils.id_str, 'clone_id', 'clonal_germline',
                          'mu_count_cdr_r', 'mu_count_cdr_s', 'mu_count_fwr_r', 'mu_count_fwr_s',
                          'cdr3_end', 'sequences']

    df_columns = df.columns.tolist()

    # Collecting the df columns by the defined order
    for col in mut_utils.all_fields:
        if col in df_columns:
            new_columns.append(col)

    # Collecting the df columns of the selection test by the defined order
    for col in all_shazam_columns:
        if col in df_columns:
            new_shazam_columns.append(col)

    for col in df_columns:
        if (col not in new_columns) and (col not in new_shazam_columns):
            new_columns.append(col)
            log_object.warning(f'Appending the unrecognized column {col} to the output database')  # TEMP

    dfg = df[new_columns]
    dfs = df[new_shazam_columns]

    return dfg, dfs


def count_obs_func(args):
    """
    The main function - counts the mutations in all given trees.
    :param args: The input arguments
    :return: None
    """

    log_object = Logger(args.name, args.silent, 'mutations', args=args)

    linked_trees = []

    if args.json:
        import utils.json_utils as json_utils
        linked_trees = json_utils.collect_trees_and_seqs_from_scheme(args, log_object)

    if args.database:
        linked_trees = db_utils.attach_seqs_by_df(args, log_object)

    if args.fasta:
        linked_trees = fasta_utils.attach_trees_using_fasta(args, log_object)

    if args.sample:
        linked_trees = general_utils.sample_tree_list(linked_trees, args.sample, log_object)

    linked_trees = [i for i in linked_trees if i]  # Remove all None values
    log_object.info(f'Finish collecting and linking {len(linked_trees)} trees')

    # Analysing trees in parallel
    args.ncors = int(args.ncors)
    number_of_workers = min(args.ncors, len(linked_trees))
    if len(linked_trees) > 0:
        with mp.Pool(number_of_workers) as pool:
            # mut_dict_list = pool.map(mut_utils.count_observed_in_tree, linked_trees)
            mut_dict_list = pool.starmap(mut_utils.count_observed_in_tree,
                                         zip(linked_trees, repeat(args.notrunk),
                                             repeat(args.nocdr3), repeat(log_object)))

    else:
        general_utils.handle_errors('Could not link alignment to trees', log_object)
        mut_dict_list = []

    # Converting to dataframes
    mut_df = pd.DataFrame(mut_dict_list)
    mut_df, selection_df = order_mut_df_columns(mut_df, log_object, args.name)
    mut_df = tree_utils.clean_db_from_none_vals(mut_df)
    selection_df = tree_utils.clean_db_from_none_vals(selection_df)

    dir_name = general_utils.create_output_dir(args.name, log_object=log_object)

    # Saving to csv tables
    if not mut_df.empty:
        file_name = f'{dir_name}/{args.name}_mutations.csv'
        mut_df.to_csv(file_name, index=False)
        log_object.info(f'Output saved as {file_name}')

        if args.selection:
            file_name = f'{dir_name}/{args.name}_for_shazam.csv'
            selection_df.to_csv(file_name, index=False)
            log_object.info(f'Output saved as {file_name}')

        if args.plot:
            from utils import mut_plot
            mut_plot.plot_mut_data(mut_df, args.name, log_object=log_object)
            log_object.info(f'Plots saved in the directory: {dir_name}/mutation_plots')
    else:
        log_object.warning('No mutations found')

    log_object.done()


def get_arg_parser():
    """
    Defines the input and the output field help message
    :return: The parser
    """

    desc = "Counts observed mutations"
    accessible_cores = mp.cpu_count()  # IgBlast uses 5 parallel processes

    # Define argument parser
    local_parser = ArgumentParser(description=desc)
    local_parser.add_argument('--version', action='version', version='%(prog)s:' + ' %s:%s' % (__version__, __date__))
    local_parser.add_argument('-t', '--tree', action='store', nargs='+', help='The trees in Newick format',
                              required=True)
    # local_parser.add_argument('-s', '--seq', action='store', nargs='+', help='The txts file', required=True)
    local_parser.add_argument('-d', '--db', action='store', help='The tab file', required=True)
    local_parser.add_argument('-n', '--name', action='store', help='A sample name. '
                                                                   'If the sample name exists - '
                                                                   'files will be overwritten', required=True)
    # local_parser.add_argument('-cf', '--clone_field', help='The clone field. default = ' + _clone_field)
    # local_parser.add_argument('-sf', '--seq_field', help='The sequence field. default = ' + _seq_field)
    # local_parser.add_argument('-gf', '--gl_field', help='The germline field. default = ' + _gl_field)
    # local_parser.add_argument('-if', '--id_field', help='The sequence ID field. default = ' + _id_field)
    local_parser.add_argument('--plot', action='store_true', help='Creates plots')
    local_parser.add_argument('--ncors', help='The number of CPU cores. default ' + str(accessible_cores),
                              default=accessible_cores)
    local_parser.add_argument('--silent', action='store_true', help='Avoid output')
    return local_parser


if __name__ == "__main__":
    # Parse command line arguments
    parser = get_arg_parser()
    global_args = parser.parse_args()

    count_obs_func(global_args)
