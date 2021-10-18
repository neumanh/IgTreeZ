#!/usr/bin/python3
"""
Counts the mutations in Newick trees.
"""

# Info
__author__ = 'Hadas Neuman'
__version__ = '1.0'
__date__ = '13/10/20'

# Imports
from argparse import ArgumentParser

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams

matplotlib.use('Agg')

rcParams.update({'figure.autolayout': True})  # Adjusting the x labels o the figure


def box_plot(dicti: dict, fig_name: str = 'figure.pdf', xlab: str = 'Transitions', ylab: str = 'Distance (mutations)'):
    """
    Plots the box plot
    :param fig_name: The output figure name
    :param dicti: The transition / population dictionary
    :param xlab: The x lable of the created graph
    :param ylab: The y lable of the created graph
    :return: None
    """

    # Collect the data
    all_data, labels = dict_to_lists(dicti)
    # print(all_data)
    # print(labels)

    if len(labels) > 1:
        fig, ax = plt.subplots()

        # rectangular box plot with customized median lines
        medianprops = dict(linewidth=1, color='black')
        bplot = ax.boxplot(all_data,
                           vert=True,  # vertical box alignment
                           patch_artist=True,  # fill with color
                           labels=labels,  # will be used to label x-ticks
                           medianprops=medianprops)

        # fill with colors
        colors = get_color_arr(labels)
        for patch, color in zip(bplot['boxes'], colors):
            patch.set_facecolor(color)

        # Adjusting the plot for large number of populations
        if len(labels) > 20:
            plt.setp(ax.get_xticklabels(), rotation=90)  # Rotate in 30 degrees
            plt.gcf().subplots_adjust(bottom=0.3)  # Add the bottom space
        elif len(labels) > 15:
            plt.setp(ax.get_xticklabels(), rotation=45,
                     horizontalalignment='right')  # Rotate in 90 degrees (vertically)
            plt.gcf().subplots_adjust(bottom=0.3)  # Add the bottom space
        elif len(labels) > 4:
            plt.setp(ax.get_xticklabels(), rotation=30,
                     horizontalalignment='right')  # Rotate in 90 degrees (vertically)
            plt.gcf().subplots_adjust(bottom=0.15)  # Add the bottom space

        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)

        # Adjusting the limits
        plt.gcf().subplots_adjust(left=0.15)

        plt.savefig(fig_name)
        plt.close()

    else:
        print("Could not create a boxplot for less than 2 populations")


def bar_plot(dicti: dict, fig_name: str = 'bar_plot.pdf', xlab='Populations', ylab='Count', sum_flag=False):
    """
    Plots the bar plot of transition number
    :param fig_name: The output figure name
    :param dicti: The transition / population dictionary
    :param xlab: The x lable of the created graph
    :param ylab: The y lable of the created graph
    :param sum_flag: Indicated weather to sum the data rather than counting it
    :return: None
    """

    # Collecting the data
    all_data, labels = dict_to_lists(dicti)
    count_data = []
    for data in all_data:  # Collecting the count data
        if sum_flag:
            count_data.append(np.sum(data))
        else:  # Count the data
            count_data.append(len(data))

    fig, ax = plt.subplots()
    ax.bar(labels, count_data)

    # Adjusting the plot for large number of populations
    if len(labels) > 20:
        plt.setp(ax.get_xticklabels(), rotation=90)  # Rotate in 30 degrees
        plt.gcf().subplots_adjust(bottom=0.5)  # Add the bottom space
    elif len(labels) > 15:
        plt.setp(ax.get_xticklabels(), rotation=45, horizontalalignment='right')  # Rotate in 90 degrees (vertically)
        plt.gcf().subplots_adjust(bottom=0.3)  # Add the bottom space
    elif len(labels) > 4:
        plt.setp(ax.get_xticklabels(), rotation=30, horizontalalignment='right')  # Rotate in 90 degrees (vertically)
        plt.gcf().subplots_adjust(bottom=0.15)  # Add the bottom space

    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)

    # Adjusting the limits
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.gcf().subplots_adjust(left=0.15)

    plt.savefig(fig_name)

    plt.close()


def get_color_arr(pops: list):
    """
    Gets the color array
    :param pops:
    :return:
    """
    colors = []
    delta = 1 / len(pops)
    col = 0.0 + (delta / 2)
    color_palt = plt.cm.Spectral

    for _ in pops:
        colors.append(color_palt(col))
        col = col + delta
    return colors


def norm_bar_plot(dicti: dict, fig_name: str = 'norm_bar_plot.pdf'):
    """
    Plots the bar plot of the normalized transition number
    :return: None
    """
    # Collect the data
    count_data, labels = dict_to_lists(dicti)

    fig, ax = plt.subplots()
    ax.bar(labels, count_data)

    xlab = 'Normalized Transitions'

    # Adjusting the plot for large number of populations
    if len(labels) > 20:
        plt.setp(ax.get_xticklabels(), rotation=90)  # Rotate in 30 degrees
    elif len(labels) > 15:
        plt.setp(ax.get_xticklabels(), rotation=45, horizontalalignment='right')  # Rotate in 90 degrees (vertically)
    elif len(labels) > 4:
        plt.setp(ax.get_xticklabels(), rotation=30, horizontalalignment='right')  # Rotate in 90 degrees (vertically)

    ax.set_xlabel(xlab)
    ax.set_ylabel('Count')

    # Adjusting the limits
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.gcf().subplots_adjust(left=0.15)

    plt.savefig(fig_name)
    plt.close()


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

    for pop in pops:
        color_dict[pop] = color_palt(col)
        col = col + delta
    return color_dict


def sunburst_pie_plot(trans_dict: dict, pop_dict: dict, fig_name: str = 'pie_plot.pdf'):
    """
    Plot the sun-burst plot
    :param pop_dict: The population dictionary
    :param fig_name: The output figure name
    :param trans_dict: The transition dictionary
    :return: None
    """
    in_labels = []
    in_values = np.array([])
    out_labels = []
    out_values = np.array([])

    colors_dict = create_color_dict(pop_dict)
    in_colors = []
    out_colors = []

    # print(trans_dict)

    # Creating the two arrays
    for pop1 in pop_dict:
        temp_trans_list = np.array([])
        in_labels.append(pop1)
        in_colors.append(colors_dict[pop1])
        for pop2 in pop_dict:
            if pop1 != pop2:
                out_labels.append(pop2)  # Saving the label
                out_colors.append(colors_dict[pop2])  # Saving the color
                trans_count = 0  # The default
                if (pop1, pop2) in trans_dict:
                    if isinstance(trans_dict[(pop1, pop2)], np.ndarray):
                        trans_count = len(trans_dict[(pop1, pop2)])
                    elif isinstance(trans_dict[(pop1, pop2)], float):
                        trans_count = trans_dict[(pop1, pop2)]  # For the normalized list
                # Appending the transition counter to the temporary transition list
                temp_trans_list = np.append(temp_trans_list, trans_count)
        out_values = np.append(out_values, temp_trans_list.flatten())
        in_values = np.append(in_values, (sum(temp_trans_list)))
    #
    # print('in_values:', in_values)
    # print('out_values:', out_values)

    # Plotting the pies

    fig, ax = plt.subplots()
    size = 0.3

    # The outer circle
    ax.pie(out_values, radius=1, colors=out_colors,
           wedgeprops=dict(width=size, edgecolor='w'))

    # The inner circle
    mypie, _ = ax.pie(in_values, radius=1 - size, colors=in_colors,
                      wedgeprops=dict(width=size, edgecolor='w'))

    # # Changeing the plot location to make room to the legend
    # chartBox = ax.get_position()
    # ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.7, chartBox.height])

    # Add legend
    ax.legend(mypie, in_labels, loc='upper right', bbox_to_anchor=(1.2, 0.8))

    plt.savefig(fig_name)
    plt.close()


def basic_pie_plot(dicti: dict, fig_name: str = 'basic_pie.pdf', sum_flag=False):
    """
    Plots the bar plot of transition number
    :param dicti: The transition / population dictionary
    :param fig_name: The output figure name
    :param sum_flag: A flag to indicate whether to sum-up the dictionary or to present it length.
    :return: None
    """

    # Collecting the data
    all_data, labels = dict_to_lists(dicti)
    count_data = []
    for data in all_data:  # Collecting the count data
        if sum_flag:  # Summarize the data
            count_data.append(np.sum(data))
        else:  # Count the data
            count_data.append(len(data))

    fig, ax = plt.subplots()
    ax.pie(count_data, labels=labels, autopct='%1.1f%%',
           shadow=True, startangle=90, colors=get_color_arr(labels))
    ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

    plt.tight_layout()
    plt.savefig(fig_name)

    plt.close()


def dict_to_lists(dic: dict):
    """
    Converts the dictionary into 2 lists - the keys (labels) the the values
    :param dic: The input dictionary
    :return: The values list and the keys list
    """
    labels = []
    values = []
    for key in dic:
        if isinstance(key, tuple):
            label = f'{key[0]}->{key[1]}'
        else:
            label = key
        labels.append(label)
        values.append(dic[key])

    return values, labels


def simple_scatter_plot(xvals, yvals, xlab, ylab, fig_name: str = 'scatter_plot.pdf'):
    """
    Plots a scatter plot with a trend line.
    :param xvals: The X values
    :param yvals: The Y values
    :param xlab: The X lables
    :param ylab: The Y label
    :param fig_name: The output figure name
    :return: None
    """
    plt.scatter(xvals, yvals)
    # plt.tight_layout()

    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.savefig(fig_name)

    plt.close()


def scatter_plot(vals: list, xlab, ylab, fig_name: str = 'scatter_plot.pdf'):
    """
    Plots a scatter plot with a trend line.
    :param vals: list of sample-names, x-values, and y-values
    :param xlab: The X lables
    :param ylab: The Y label
    :param fig_name: The output figure name
    :return: None
    """

    fig, ax = plt.subplots()

    labels = []
    for val in vals:
        label = val[0]
        x = val[1]
        y = val[2]

        ax.scatter(x, y, label=label, alpha=0.5)

    # plt.tight_layout()
    ax.legend()
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    plt.savefig(fig_name)

    plt.close()


def plot_box(args):
    """
    Plots a box plot of the input CSV file.
    :param args: The input arguments.
    :return: None
    """
    df = pd.read_csv(args.csv)
    dic = {}
    for col in df:
        dic[col] = df[col].dropna().to_numpy()

    fig_name = f'{args.csv}.pdf'

    if args.sort:
        new_dic = {}
        for s_col in args.sort:
            if s_col in dic:
                new_dic[s_col] = dic[s_col]
            else:
                print(f'Warning: could not find the column {s_col} in the database')
        dic = new_dic

    box_plot(dic, fig_name, xlab='Population')

    print('Saved to', fig_name)


def get_arg_parser():
    """
    Defines the input and the output field help message
    :return: The parser
    """

    desc = "Plots box from a CSV file"

    # Define argument parser
    local_parser = ArgumentParser(description=desc)
    local_parser.add_argument('-c', '--csv', help='The CSV file to plot', required=True)
    local_parser.add_argument('-s', '--sort', help='Sort the columns in the input order', nargs='+')
    local_parser.add_argument('--version', action='version', version='%(prog)s:' + ' %s:%s' % (__version__, __date__))

    return local_parser


if __name__ == "__main__":
    # Parse command line arguments
    parser = get_arg_parser()
    global_args = parser.parse_args()

    plot_box(global_args)
