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
from utils import plot_utils, general_utils, stat_utils

_test_str = 'test'
_sample_str = 'sample'


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


def create_box_plots_and_stats(dfs: list, box_fields: list, out_dir: str, name: str, log_object: Logger,
                               font: str = None, subsample: bool = False):
    """
    Creates the box plots and the statistic analysis.
    :param dfs: The dataframes list
    :param box_fields: The field to analyze
    :param out_dir: The output directory
    :param name: the analysis name
    :param log_object: The log object
    :param font: Font size
    :param subsample: Whether to use the subsample column
    :return: None
    """
    log_object.info('Creating the box plots and the statistical analysis.')

    pvals = {}
    pval_pair_dict = {}
    fdr_pair_dict = {}

    if subsample:
        subgroup = True
        for df in dfs:
            if not ('subsample' in list(df)):
                subgroup = False
    else:
        subgroup = False

    stat_file_name = f'{out_dir}/{name}_statistics.txt'  # Defining the output name
    stat_table_name_pvals = f'{out_dir}/{name}_statistics_2_samples.csv'  # Defining the output name
    stat_table_name_fdr = f'{out_dir}/{name}_statistics_2_samples_fdr.csv'  # Defining the output name

    general_utils.remove_files([stat_file_name])
    grouped_box_df = pd.DataFrame()

    for field in box_fields:
        if field not in [general_utils.sample_str, general_utils.id_str, 'id']:  # Ignore the sample and id column
            # print("Comparing by the field", field)
            dic = {}

            for df in dfs:  # Collecting the data of the specific columns
                key = df.at[0, general_utils.sample_str]  # The sample name
                if subgroup:
                    df = df.groupby('subsample').mean()
                    filename = f'{out_dir}/{key}_grouped.csv'  # Defining the output name
                    df = df.sort_values(by=['subsample']).reset_index()  # Removing the row index
                    df.to_csv(filename, index=False)
                    df[general_utils.sample_str] = key
                if (field in list(df)) and (pd.api.types.is_numeric_dtype(df[field])):  # A numeric value
                    field_data = df[field].dropna()
                    dic[key] = field_data.to_numpy()
            if dic:
                figname = f'{out_dir}/{name}_{field}.pdf'  # Defining the output name
                plot_utils.box_plot(dic, fig_name=figname, xlab='Sample', ylab=field, font=font)
                table_name = f'{out_dir}/{name}_{field}.csv'  # Defining the output name
                sum_df = pd.DataFrame.from_dict(dic, orient='index')
                sum_df = sum_df.transpose()
                sum_df.to_csv(table_name, index=False)
                grouped_box_df = sum_dfs_to_grouped_box_plot(grouped_box_df, sum_df, field)
                pval, pval_dict, fdr_dict = None, None, None

                try:
                    pval, pval_dict, fdr_dict = stat_utils.analyze_df(sum_df, name=field,
                                                                      out_file=stat_file_name, log_object=log_object)
                except Exception as e:
                    log_object.error(f'Could not perform statistical analysis on the field {field}: {e}')

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

    # Plotting the grouped boxplot
    figname = f'{out_dir}/grouped_box_plot.pdf'  # Defining the output name
    plot_utils.grouped_box_plot(grouped_box_df, x=_test_str, y='value', hue=_sample_str, fig_name=figname)


def test_input(args):
    """
    Tests if the input is valid
    :param args: The input argument.
    :return: The argument error type. If no error was found - return None
    """
    arg_error = None

    # Check if the input files exist
    missing_file = general_utils.test_file_existence(args.csv)
    if missing_file:
        arg_error = f'The file {missing_file} does not exist'
    return arg_error


def collect_dfs(files: list, fields: list):
    """
    Collects the dataframes and the relevant fields
    :param files: The input files.
    :param fields: The fields to analyzed.
    :return: A list of dataframes
    """
    dfs = []
    # Collecting the csv tables
    for file in files:
        temp_df = pd.read_csv(file)
        dfs.append(temp_df)

    field_set = set([])
    for df in dfs:
        field_set = field_set.union(set(list(df)))

    if fields:
        for field in fields:
            if field not in field_set:
                general_utils.handle_errors(f'The field {field} was not found in any database. Exiting', exit_stat=True)
    else:
        fields = list(field_set)  # The column names of all dataframes

    return dfs, fields, field_set


def filter_dfs(dfs: list, filter_field: str, filter_val: list, field_set: set):
    """
    Filter tha dataframes in a filter field is given
    :param dfs: The dataframes.
    :param filter_field: The fields to filter by.
    :param filter_val: The values of the field.
    :param field_set: All the dataframe fields.
    :return: A list of dataframes
    """
    # Filtering the dataframes
    if filter_field or filter_val:
        if not (filter_field and filter_val):
            general_utils.handle_errors('The filter_field and filter_val arguments must be used together.',
                                        exit_stat=True)
        elif len(filter_val) != 2:
            general_utils.handle_errors('The filter_val must receive minimum and maximum values.', exit_stat=True)

        elif filter_field not in field_set:
            general_utils.handle_errors(f'The filtering field {filter_field} was not found.', exit_stat=True)

        min_node_number = filter_val[0]
        max_node_number = filter_val[1]

        temp_dfs = []
        for df in dfs:
            if filter_field in list(df):
                filtered_df = df[
                    (df[filter_field] >= int(min_node_number)) & (df[filter_field] <= int(max_node_number))]
                temp_dfs.append(filtered_df)
            dfs = temp_dfs
        return dfs


def sample_dfs(dfs: list):
    """
    Samples the dataframes by the shortest
    :param dfs: The dataframes.
    :return: A list of dataframes
    """
    sample_size = find_shortest_df(dfs)
    temp_dfs = []
    for df in dfs:
        short_df = df.sample(n=sample_size)
        temp_dfs.append(short_df)
    dfs = temp_dfs
    return dfs


def comp_reps(args):
    """
    Compares mtree results in graphical representation
    :param args: The input argumens
    :return: None
    """
    log_object = Logger(args.name, args.silent, 'compare', args)

    arg_error = test_input(args)
    if arg_error:
        general_utils.handle_errors(arg_error, exit_stat=True)

    # Collecting the dataframes
    dfs, fields, field_set = collect_dfs(args.csv, args.field)

    # Filtering the dataframes if asked
    if args.filter_field or args.filter_val:
        dfs = filter_dfs(dfs, args.filter_field, args.filter_val, field_set)

    # Sampling the dataframes if asked
    if args.sample:
        dfs = sample_dfs(dfs)

    # Re-indexing the dataframe
    temp_dfs = []
    for df in dfs:
        indexed_df = df.reset_index(drop=True)
        temp_dfs.append(indexed_df)
    dfs = temp_dfs

    out_dir = general_utils.create_output_dir(args.name, log_object=log_object)

    # Run only if the field is numeric
    numeric_fields = fields
    box_field = numeric_fields  # TEMP

    # Boxplots
    create_box_plots_and_stats(dfs, box_field, out_dir, args.name, log_object, args.font, args.subsample)

    # Scatter plots
    if args.scatter:
        scat_fields = numeric_fields
        create_scatter_plots(dfs, scat_fields, out_dir, log_object)

    log_object.done()


def sum_dfs_to_grouped_box_plot(gdf: pd.DataFrame, df: pd.DataFrame, test_name: str):
    """
    Groups the dataframes into one big dataframes for the grouped box plot
    :param gdf: The grouped dataframe
    :param df: The test dataframe
    :param test_name: The test name
    :return: None
    """

    melted_df = df.melt()  # Melting the dataframe, so that every row represent all the columns
    melted_df[_test_str] = test_name  # Adding the test information
    melted_df.rename(columns={"variable": _sample_str}, inplace=True)  # Renaming the columns
    gdf = pd.concat([gdf, melted_df], ignore_index=True)
    return gdf


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

    local_parser.add_argument('--font', action='store', help='Font size')

    local_parser.add_argument('--sample', action='store_true', help='Samples by the minimum mtree result size')
    local_parser.add_argument('--scatter', action='store_true', help='Creates scatter plots as well')
    local_parser.add_argument('--subsample', action='store_true', help='Compares by sub-samples if the databases '
                                                                       'include a \'sabsample\' column.')

    local_parser.add_argument('--silent', action='store_true', help='Avoid output')
    return local_parser


if __name__ == "__main__":
    # Parse command line arguments
    parser = get_arg_parser()
    global_args = parser.parse_args()

    comp_reps(global_args)
