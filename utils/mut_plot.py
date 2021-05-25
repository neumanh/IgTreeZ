#!/usr/bin/python3
"""
Plots Ig statistics
"""

# Info
__author__ = 'Hadas Neuman'
__version__ = '1.61'
__date__ = '19/10/20'

# Imports
import pandas as pd
from objects import Logger
from utils import general_utils, plot_utils


def plot_mut_data(mut_df: pd.DataFrame, name, log_object: Logger = None):
    """
    Plots mutation data
    :param mut_df: A dataframe of the count_mut results
    :param name: The analysis name
    :param log_object: The Logger object
    :return: None
    """
    if not log_object:
        log_object = general_utils.create_general_log_object()

    dir_name = general_utils.create_output_dir(name, 'mutation_plots', log_object)
    basic_file_name = f'{dir_name}/{name}'
    plot_mut_array(mut_df, basic_file_name)


def df_columns_to_dict(df: pd.DataFrame, cols: list, names: list = None):
    """
    Collects the relevant columns from a given dataframe and returns them as a dictionary.
    :param df: The dataframe.
    :param cols: The dataframe columns to collect.
    :param names: The corresponding names of the groups.
    :return: A dictionary of the columns.
    """
    mut_dic = {}
    if not names:
        names = cols

    if len(names) != len(cols):
        mut_dic = None
    else:
        for idx, col in enumerate(cols):
            if col in df:
                name = names[idx]
                mut_dic[name] = df[col].dropna().to_numpy()

    return mut_dic


def plot_mut_array(mut_df: pd.DataFrame, name: str):
    """
    Plots mutation statistics.
    :param mut_df: The mutation dataframe
    :param name: The analysis name
    :return: None
    """

    # Hydrophobic/ hydrophilic source
    fig_title = '_hydro_source_'
    cols = ['source_hydrophobic', 'source_hydrophilic']
    names = ['Hydrophobic source', 'Hydrophilic source']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)

    # Hydrophobic/ hydrophilic source
    fig_title = '_hydro_dest_'
    cols = ['dest_hydrophobic', 'dest_hydrophilic']
    names = ['Hydrophobic destination', 'Hydrophilic destination']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)

    # Hydrophobic/ hydrophilic change
    fig_title = '_hydro_change_'
    cols = ['hydro_keep', 'hydro_change']
    names = ['Same hydropathy', 'Switched hydropathy']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)

    # Charge change
    fig_title = '_charge_change_'
    cols = ['charge_keep', 'charge_change']
    names = ['Same charge', 'Switched charge']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)

    # Source charge
    fig_title = '_charge_source_'
    cols = ['source_positive', 'source_neutral', 'source_negative']
    names = ['Positive source', 'Neutral source', 'Negative source']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)

    # Destination charge
    fig_title = '_charge_dest_'
    cols = ['dest_positive', 'dest_neutral', 'dest_negative']
    names = ['Positive destination', 'Neutral destination', 'Negative destination']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)

    # CDR/FWR
    fig_title = '_region_'
    cols = ['CDR', 'FWR']
    names = ['CDR region', 'FWR region']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)

    # CDR1/CDR2/CDR3
    fig_title = '_CDR_region_'
    cols = ['cdr1', 'cdr2', 'cdr3']
    names = ['CDR1 region', 'CDR2 region', 'CDR3 region']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)

    # FWR1/FWR2/FWR3
    fig_title = '_FWR_region_'
    cols = ['fwr1', 'fwr2', 'fwr3']
    names = ['FWR1 region', 'FWR2 region', 'FWR3 region']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)

    # CDR1/CDR2/CDR3/FWR1/FWR2/FWR3
    fig_title = '_FWR_CDR_region_'
    cols = ['fwr1', 'fwr2', 'fwr3', 'cdr1', 'cdr2', 'cdr3']
    names = ['FWR1 region', 'FWR2 region', 'FWR3 region', 'CDR1 region', 'CDR2 region', 'CDR3 region']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)

    # Transition / Tranversion
    fig_title = '_trans_'
    cols = ['transition', 'transversion']
    names = ['Transition', 'Tranversion']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)

    # Nucleotide source
    fig_title = '_nucleotide_source_'
    cols = ['source_a', 'source_c', 'source_g', 'source_t']
    names = ['Source A', 'Source C', 'Source G', 'Source T']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)

    # R / S
    fig_title = '_rs_'
    cols = ['r', 's']
    names = ['Replacement Mutation', 'Silent Mutation']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)

    # # CDR / FWR / R / S
    # fig_title = '_region_and_type_'
    # cols = ['obs_r_cdr', 'obs_s_cdr', 'obs_r_fwr', 'obs_s_fwr']
    # names = ['CDR-R', 'CDR-S', 'FWR-R', 'FWR-S']
    # plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)


def plot_three_plots(dic: dict, fig_title: str):
    """
    Sends the given dictionary for bar plot, box plot and pie plot.
    :param dic: The data dictionary
    :param fig_title: The analysis name
    :return: None
    """
    plot_utils.box_plot(dic, fig_name=fig_title + 'box_plot.pdf', xlab='Mutation type', ylab='Mutations per tree')
    plot_utils.bar_plot(dic, fig_name=fig_title + 'bar_plot.pdf', xlab='Mutation type', ylab='Overall mutations',
                        sum_flag=True)
    plot_utils.basic_pie_plot(dic, fig_name=fig_title + 'pie_plot.pdf', sum_flag=True)
