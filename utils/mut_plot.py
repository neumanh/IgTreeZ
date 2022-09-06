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
    cols = ['replacement', 'silent']
    names = ['Replacement Mutation', 'Silent Mutation']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)

    plot_conservation(mut_df, name)
    plot_hydro(mut_df, name)
    plot_charge(mut_df, name)
    plot_chemical(mut_df, name)
    plot_polarity(mut_df, name)
    plot_region(mut_df, name)
    plot_volume(mut_df, name)
    plot_hydro_da(mut_df, name)


def plot_volume(mut_df: pd.DataFrame, name: str):
    """
    Plots mutation statistics.
    :param mut_df: The mutation dataframe
    :param name: The analysis name
    :return: None
    """
    # Source charge
    fig_title = '_volume_source_'
    cols = ['vs_volume_source', 'small_volume_source', 'medium_volume_source', 'large_volume_source',
            'vl_volume_source']
    names = ['Very small source', 'Small source', 'Medium source', 'Large source', 'Very large source']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)

    # Target charge
    fig_title = '_volume_target_'
    cols = ['vs_volume_target', 'small_volume_target', 'medium_volume_target', 'large_volume_target',
            'vl_volume_target']
    names = ['Very small target', 'Small target', 'Medium target', 'Large target', 'Very large target']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)

    # Charge change
    fig_title = '_volume_change_'
    cols = ['volume_decrease', 'volume_increase', 'volume_keep']
    names = ['Decreased volume', 'Increased volume', 'Same volume']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)


def plot_chemical(mut_df: pd.DataFrame, name: str):
    """
    Plots mutation statistics.
    :param mut_df: The mutation dataframe
    :param name: The analysis name
    :return: None
    """
    # Source charge
    fig_title = '_chemical_source_'
    cols = ['aliphatic_source', 'aromatic_source', 'sulfur_source', 'hydroxyl_source', 'basic_source', 'acidic_source',
            'amide_source']
    names = ['Aliphatic source', 'Aromatic source', 'Sulfur_source', 'Hydroxyl source', 'Basic source', 'Acidic source',
             'Amide source']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)

    # Target charge
    fig_title = '_chemical_target_'
    cols = ['aliphatic_target', 'aromatic_target', 'sulfur_target', 'hydroxyl_target', 'basic_target', 'acidic_target',
            'amide_target']
    names = ['Aliphatic target', 'Aromatic target', 'target', 'Hydroxyl target', 'Basic target', 'Acidic target',
             'Amide target']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)

    # Charge change
    fig_title = '_chemical_change_'
    cols = ['chemical_keep', 'chemical_change']
    names = ['Same chemical', 'Switched chemical']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)


def plot_hydro_da(mut_df: pd.DataFrame, name: str):
    """
    Plots mutation statistics.
    :param mut_df: The mutation dataframe
    :param name: The analysis name
    :return: None
    """
    # Source hydro-donor-acceptor
    fig_title = '_hydro_donor_acceptor_source_'
    cols = ['hydro_donor_source', 'hydro_acceptor_source', 'hydro_donor_acceptor_source', 'hydro_da_none_source']
    names = ['Hydrogen donor source', 'Hydrogen acceptor source', 'Hydrogen donor and acceptor source',
             'No hydrogen donor or acceptor source']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)

    # Target hydro-donor-acceptor
    fig_title = '_hydro_donor_acceptor_target_'
    cols = ['hydro_donor_target', 'hydro_acceptor_target', 'hydro_donor_acceptor_target', 'hydro_da_none_target']
    names = ['Hydrogen donor target', 'Hydrogen acceptor target', 'Hydrogen donor and acceptor target',
             'No hydrogen donor or acceptor target']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)

    # Change
    fig_title = '_hydro_da_change_'
    cols = ['hydro_da_keep', 'hydro_da_change']
    names = ['Same hydrogen-donor-acceptor tendency', 'Switched hydrogen-donor-acceptor tendency']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)


def plot_hydro(mut_df: pd.DataFrame, name: str):
    """
    Plots mutation statistics.
    :param mut_df: The mutation dataframe
    :param name: The analysis name
    :return: None
    """
    # Hydrophobic/ hydrophilic source
    fig_title = '_hydro_source_'
    cols = ['hydrophobic_source', 'neutral_hydro_source', 'hydrophilic_source']
    names = ['Hydrophobic source', 'Neutral hydropathy source', 'Hydrophilic source']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)

    # Hydrophobic/ hydrophilic source
    fig_title = '_hydro_target_'
    cols = ['hydrophobic_target', 'neutral_hydro_target', 'hydrophilic_target']
    names = ['Hydrophobic target', 'Neutral hydropathy target', 'Hydrophilic target']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)

    # Hydrophobic/ hydrophilic change
    fig_title = '_hydro_change_'
    cols = ['hydro_keep', 'hydro_change']
    names = ['Same hydropathy', 'Switched hydropathy']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)


def plot_charge(mut_df: pd.DataFrame, name: str):
    """
    Plots mutation statistics.
    :param mut_df: The mutation dataframe
    :param name: The analysis name
    :return: None
    """
    # Charge change
    fig_title = '_charge_change_'
    cols = ['charge_keep', 'charge_change']
    names = ['Same charge', 'Switched charge']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)

    # Source charge
    fig_title = '_charge_source_'
    cols = ['positive_source', 'neutral_charge_source', 'negative_source']
    names = ['Positive source', 'Neutral source', 'Negative source']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)

    # Target charge
    fig_title = '_charge_target_'
    cols = ['positive_target', 'neutral_charge_target', 'negative_target']
    names = ['Positive target', 'Neutral target', 'Negative target']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)


def plot_region(mut_df: pd.DataFrame, name: str):
    """
    Plots mutation statistics.
    :param mut_df: The mutation dataframe
    :param name: The analysis name
    :return: None
    """
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
    fig_title = '_all_regions_'
    cols = ['fwr1', 'fwr2', 'fwr3', 'cdr1', 'cdr2', 'cdr3']
    names = ['FWR1 region', 'FWR2 region', 'FWR3 region', 'CDR1 region', 'CDR2 region', 'CDR3 region']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)


def plot_conservation(mut_df: pd.DataFrame, name: str):
    """
    Plots mutation statistics.
    :param mut_df: The mutation dataframe
    :param name: The analysis name
    :return: None
    """
    # Source
    # TODO


def plot_polarity(mut_df: pd.DataFrame, name: str):
    """
    Plots mutation statistics.
    :param mut_df: The mutation dataframe
    :param name: The analysis name
    :return: None
    """
    # Source
    fig_title = '_polarity_source_'
    cols = ['polar_source', 'non_polar_source']
    names = ['Polar source', 'Non-polar source']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)

    # Target
    fig_title = '_polarity_target_'
    cols = ['polar_target', 'non_polar_target']
    names = ['Polar target', 'Non-polar target']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)

    # Polarity change
    fig_title = '_polarity_change_'
    cols = ['polarity_keep', 'polarity_change']
    names = ['Same polarity', 'Switched polarity']
    plot_three_plots(df_columns_to_dict(mut_df, cols, names), name + fig_title)


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
