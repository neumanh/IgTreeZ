#!/usr/bin/python3

"""
Compares two or more mtree outputs
"""

# Info
__author__ = 'Hadas Neuman'
__version__ = '1.1'
__date__ = '15/9/2020'

from argparse import ArgumentParser

import numpy as np
import pandas as pd
from objects import Logger
# Project imports
from utils import plot_utils, general_utils, stat_utils


def find_shortest_df(dfs: list):
    """
    Finds the length of the shortest dataframe
    :param dfs: The dataframes
    :return: The length of the shortest dataframe
    """
    min_df = np.Inf  # Infinity
    for df in dfs:
        curr_df_len = df.shape[0]  # Number of rows
        if curr_df_len < min_df:
            min_df = curr_df_len  # Save the smallest value
    return min_df


def create_scatter_plots(dfs: list, scat_fields: list, out_dir: str, log_object: Logger):
    """
    Creates scatter plots for each pair of variables.
    :param dfs: The dataframes list
    :param scat_fields: The field to analyze
    :param out_dir: The output directory
    :param log_object: The log object
    :return: None
    """

    log_object.info('Plotting scatter plots')

    for idx, fx in enumerate(scat_fields):
        for idy, fy in enumerate(scat_fields):
            if idx < idy:  # Each pair is handled once
                if (fx not in [general_utils.sample_str, 'id', 'tree_id']) and (
                        fy not in [general_utils.sample_str, 'id', 'tree_id']):  # Ignore the sample and id column
                    values = []
                    # print('Comparing', fy, 'to', fx)
                    for df in dfs:  # Collecting the data of the specific columns
                        if (fx in df) and (fy in df):
                            if (pd.api.types.is_numeric_dtype(df[fx]) and pd.api.types.is_numeric_dtype(
                                    df[fy])):  # Numeric values
                                curr_val = []
                                sample_name = df.at[0, general_utils.sample_str]
                                xvals = df[fx].to_numpy()
                                yvals = df[fy].to_numpy()
                                curr_val.append(sample_name)
                                curr_val.append(xvals)
                                curr_val.append(yvals)
                                values.append(curr_val)
                            if values:
                                figname = f'{out_dir}/{fx}_to_{fy}.pdf'  # Defining the output name
                                plot_utils.scatter_plot(values, xlab=fx, ylab=fy, fig_name=figname)


def create_box_plots_and_stats(dfs: list, box_fields: list, out_dir: str, name: str, log_object: Logger):
    """
    Creates the box plots and the statistic analysis.
    :param dfs: The dataframes list
    :param box_fields: The field to analyze
    :param out_dir: The output directory
    :param name: the analysis name
    :param log_object: The log object
    :return: None
    """
    log_object.info('Creating the box plots and the statistical analysis.')

    pvals = {}
    pval_pair_dict = {}
    fdr_pair_dict = {}

    stat_file_name = f'{out_dir}/{name}_statistics.txt'  # Defining the output name
    stat_table_name_pvals = f'{out_dir}/{name}_statistics_2_samples.csv'  # Defining the output name
    stat_table_name_fdr = f'{out_dir}/{name}_statistics_2_samples_fdr.csv'  # Defining the output name

    general_utils.remove_files([stat_file_name])

    for field in box_fields:
        if field not in [general_utils.sample_str, general_utils.id_str, 'id']:  # Ignore the sample and id column
            # print("Comparing by the field", field)
            dic = {}

            for df in dfs:  # Collecting the data of the specific columns
                if (field in list(df)) and (pd.api.types.is_numeric_dtype(df[field])):  # A numeric value
                    key = df.at[0, general_utils.sample_str]
                    field_data = df[field].dropna()
                    dic[key] = field_data.to_numpy()
            if dic:
                figname = f'{out_dir}/{name}_{field}.pdf'  # Defining the output name
                plot_utils.box_plot(dic, fig_name=figname, xlab='Sample', ylab=field)
                table_name = f'{out_dir}/{name}_{field}.csv'  # Defining the output name
                sum_df = pd.DataFrame.from_dict(dic, orient='index')
                sum_df = sum_df.transpose()
                sum_df.to_csv(table_name, index=False)
                pval, pval_dict, fdr_dict = None, None, None

                try:
                    pval, pval_dict, fdr_dict = stat_utils.analyze_df(sum_df, name=field, out_file=stat_file_name,
                                                            log_object=log_object)
                except Exception as e:
                    log_object.warning(f'Could not perform statistical analysis on the field {field}: {e}')

                if pval is not None:
                    pvals[field] = pval
                if pval_dict is not None:
                    pval_pair_dict[field] = pval_dict
                if fdr_dict is not None:
                    fdr_pair_dict[field] = fdr_dict

    pval_df = pd.DataFrame(pval_pair_dict)
    pval_df.to_csv(stat_table_name_pvals)
    fdr_df = pd.DataFrame(fdr_pair_dict)
    fdr_df.to_csv(stat_table_name_fdr)


def comp_mtree(args):
    """
    Compares mtree results in graphical representation
    :param args: The input argumens
    :return: None
    """

    dfs = []
    log_object = Logger(args.name, args.silent, 'compare', args)

    # Collecting the mtree tables
    for mtree_file in args.csv:
        temp_df = pd.read_csv(mtree_file)
        dfs.append(temp_df)

    field_set = set([])
    for df in dfs:
        field_set = field_set.union(set(list(df)))

    if args.field:
        fields = args.field

        for field in fields:
            if field not in field_set:
                general_utils.handle_errors(f'The field {field} was not found in any database. Exiting', exit_stat=True)
    else:
        fields = list(field_set)  # The column names of the first dataframe

    # Filtering the dataframes
    if args.filter_field or args.filter_val:
        if not (args.filter_field and args.filter_val):
            general_utils.handle_errors('The filter_field and filter_val arguments must be used together.',
                                        exit_stat=True)
        elif len(args.filter_val) != 2:
            general_utils.handle_errors('The filter_val must receive minimum and maximum values.', exit_stat=True)

        elif args.filter_field not in field_set:
            general_utils.handle_errors(f'The filtering field {args.filter_field} was not found.', exit_stat=True)

        min_node_number = args.filter_val[0]
        max_node_number = args.filter_val[1]

        temp_dfs = []
        for df in dfs:
            if args.filter_field in list(df):
                filtered_df = df[
                    (df[args.filter_field] >= int(min_node_number)) & (df[args.filter_field] <= int(max_node_number))]
                temp_dfs.append(filtered_df)
            dfs = temp_dfs

    # Sampling the dataframes
    if args.sample:
        sample_size = find_shortest_df(dfs)
        temp_dfs = []
        for df in dfs:
            short_df = df.sample(n=sample_size)
            temp_dfs.append(short_df)
        dfs = temp_dfs

    # Re-indexing the dataframe
    temp_dfs = []
    for df in dfs:
        indexed_df = df.reset_index(drop=True)
        temp_dfs.append(indexed_df)
    dfs = temp_dfs

    # TODO: Check input sanity

    out_dir = general_utils.create_output_dir(args.name, log_object=log_object)

    # TODO:  Check if the field is in the database

    # Run only if the field is numeric

    numeric_fields = fields
    box_field = numeric_fields  # TEMP

    # Boxplots
    create_box_plots_and_stats(dfs, box_field, out_dir, args.name, log_object)

    # Scatter plots
    scat_fields = numeric_fields
    create_scatter_plots(dfs, scat_fields, out_dir, log_object)

    log_object.done()


def get_arg_parser():
    """
    Defines the input and the output field help message
    :return: The parser
    """

    desc = "Counts observed mutations"

    # Define argument parser
    local_parser = ArgumentParser(description=desc)
    local_parser.add_argument('--version', action='version', version='%(prog)s:' + ' %s:%s' % (__version__, __date__))
    local_parser.add_argument('-c', '--csv', action='store', nargs='+', help='The mtree outputs in csv format',
                              required=True)
    local_parser.add_argument('-n', '--name', action='store', help='The analysis name', required=True)
    local_parser.add_argument('-f', '--field', action='store', nargs='+', help='The comparison field. Default - all')

    local_parser.add_argument('-ff', '--filter_field', action='store',
                              help='Filters by the given filed. To use with the -filter_val argument')
    local_parser.add_argument('-fv', '--filter_val', action='store', nargs='+',
                              help='Filters by the given filed values (min max). '
                                   'To use with the -filter_field argument')

    local_parser.add_argument('--sample', action='store_true', help='Samples by the minimum mtree result size')

    local_parser.add_argument('--silent', action='store_true', help='Avoid output')
    return local_parser


if __name__ == "__main__":
    # Parse command line arguments
    parser = get_arg_parser()
    global_args = parser.parse_args()

    comp_mtree(global_args)
