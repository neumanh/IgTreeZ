#!/usr/bin/python3
"""
Plots Ig statistics
"""

# Info
__author__ = 'Hadas Neuman'
__mail__ = 'hadas.doron@gmail.com'
__version__ = '1.7.0'
__date__ = '3/8/21'

import multiprocessing as mp
from argparse import ArgumentParser


def check_arguments(args):
    """
    Test tha arguments' validity
    :param args: The input arguments
    :return: validity (true or false)
    """
    # Check mutual exclusion
    valid = False
    if (args.command != 'compare') and args.tree and args.json:
        print('Error: TREE and JSON: only one of the arguments can be used')
    elif (args.command == 'mutations') and args.fasta and args.database:
        print('Error: FASTA and DATABASE: only one of the arguments can be used')
    elif (args.command == 'mutations') and args.json and args.database:
        print('Error: JSON and DATABASE: only one of the arguments can be used')
    elif (args.command != 'compare') and (not args.tree) and (not args.json):
        print('Error: Trees must given as Newick files or as part of an AIRR scheme json file')
    elif (args.command == 'mutations') and not args.json and (not args.database) and (not args.fasta):
        print('Error: Sequences must given as Fasta files or as ChangeO/AIRR database, or as AIRR scheme json file')
    elif (args.command == 'mutations') and (not args.database) and (args.clone_field or args.gl_field or args.seq_field
                                                                    or args.id_field):
        print('Error: The FIELD arguments can only be used with the DATABASE argument')
    else:
        valid = True

    return valid


def get_parent():
    """
    Gets the parent = the common arguments for all subcommands
    :return: The parent arguments
    """
    parser = ArgumentParser(add_help=False)

    accessible_cores = mp.cpu_count()

    parser.add_argument('-t', '--tree', action='store', nargs='+', help='The trees in Newick format')
    parser.add_argument('-j', '--json', action='store',
                        help='A jason file in the AIRR \'Clone and Lineage Tree Schema\' format,'
                             'as described in '
                             'https://docs.airr-community.org/en/stable/datarep/clone.html')
    parser.add_argument('-n', '--name', action='store', help='The analysis name', required=True)
    parser.add_argument('--sample', help='Samples the input tree list by the given number')
    parser.add_argument('--ncors', help='The number of CPU cores. default ' + str(accessible_cores),
                        default=accessible_cores)
    parser.add_argument('--silent', action='store_true', help='Avoid output')
    parser.add_argument('--version', action='version', version='%(prog)s:' + ' %s %s' % (__version__, __date__))

    return parser


def get_arg_parser():
    """
    Defines the input and the output field help message
    :return: The parser
    """

    desc = "Plot analysis results"

    # Define argument parser
    parents = get_parent()

    parser = ArgumentParser(description=desc)
    subparsers = parser.add_subparsers(title='subcommands', dest='command', help='Analysis type')
    subparsers.required = True

    # The mutation count
    parser_mut = subparsers.add_parser('mutations', parents=[parents],
                                       help='Plot the IgTreez mutation outputs.')
    parser_mut.add_argument('-f', '--fasta', action='store', nargs='+',
                            help='Fasta files for each clone\'s sequences and germinal sequence')
    parser_mut.add_argument('-d', '--database', action='store',
                            help='A database in the ChangeO or AIRR format contians the sequences and germline '
                                 'sequences. Obligatory fields:\n'
                                 '\t SEQUENCE_IMGT / sequence_alignment \n'
                                 '\t GERMLINE_IMGT_D_MASK / germline_alignment \n'
                                 '\t CLONE / clone_id \n'
                                 '\t SEQUENCE_ID / sequence_id \n'
                                 'If the CDR3_IMGT / cdr3_imgt column also exists and the \'--nocdr3\' parameter '
                                 'is avoided - defines the CDR3 region as well.')
    parser_mut.add_argument('-dbf', '--dbformat', help='The database format (AIRR/ChangeO). default - changeo',
                            choices=['airr', 'changeo'], default='changeo')
    parser_mut.add_argument('--notrunk', action='store_true',
                            help='Do not analyze nodes that descend directly from the root. '
                                 'To avoid incorrect mutation profiling resulting from an incorrect annotation.')
    parser_mut.add_argument('--nocdr3', action='store_true',
                            help='Do not analyze the CDR3 region. Relevant only with the -d parameter.')
    parser_mut.add_argument('--selection', action='store_true',
                            help='Create an output table for the selection analysis using the ShazaM package')
    parser_mut.add_argument('--plot', action='store_true', help='Creates plots')
    parser_mut.add_argument('--illumina', action='store_true', help='If the sequence names contain colons (:) '
                                                                    'or semicolons (;), they will be replaced '
                                                                    'in under-lines')
    parser_mut.add_argument('-cf', '--clone_field', help='The clone field. default - as in ChangeO / AIRR format')
    parser_mut.add_argument('-sf', '--seq_field', help=f'The sequence field. default - as in ChangeO / AIRR format')
    parser_mut.add_argument('-gf', '--gl_field', help=f'The germline field. default - as in ChangeO / AIRR format')
    parser_mut.add_argument('-if', '--id_field', help=f'The sequence ID field. default - as in ChangeO / AIRR format')
    parser_mut.add_argument('-gl', '--gl_name', help=f'The germline name. default - GL', default='GL')
    parser_mut.set_defaults(function=call_mutations)

    # The mtree analysis
    parse_mtree = subparsers.add_parser('mtree', parents=[parents],
                                        help='Describes the trees\' shape properties.')
    parse_mtree.set_defaults(function=call_mtree)

    # The poptree analysis
    parse_poptree = subparsers.add_parser('poptree', parents=[parents],
                                          help='Characterises the populations in '
                                               'trees and the transitions between them')
    parse_poptree.add_argument('-p', '--pops', action='store', nargs='+',
                               help='The population to search for in the trees',
                               required=True)
    parse_poptree.add_argument('--plot', action='store_true', help='Plots the populations statistics')
    parse_poptree.add_argument('-z', '--zero', action='store_true', help='Analyse zero distance transitions as well')

    parse_poptree.set_defaults(function=call_poptree)
    # The filter analysis
    parse_filter = subparsers.add_parser('filter', parents=[parents],
                                         help='Draws the trees using the dot command')
    parse_filter.add_argument('-AND', action='store', nargs='+', help='Collect trees with all the given populations.')
    parse_filter.add_argument('-OR', action='store', nargs='+', help='Collect trees with one of the given populations.')
    parse_filter.add_argument('-NOT', action='store', nargs='+',
                              help='Collect trees without the given populations.')
    parse_filter.add_argument('-nodes', action='store', nargs='+',
                              help='Collect trees with number of nodes in the size limits')
    parse_filter.add_argument('-leaves', action='store', nargs='+',
                              help='Collect trees with number of leaves in the size limits')

    parse_filter.add_argument('--copy', action='store_true',
                              help='If given, copies the collected trees in addition to listing them.')
    parse_filter.set_defaults(function=call_filter)

    # The draw analysis
    parse_draw = subparsers.add_parser('draw', parents=[parents],
                                       help='Draws the trees using the dot command')
    parse_draw.add_argument('-p', '--pops', action='store', nargs='+',
                            help='The population to search for in the trees')
    parse_draw.add_argument('--dot', help='The dot program path. default: dot', default='dot')
    parse_draw.add_argument('--format', help='The output format. default: png', default='png')
    parse_draw.add_argument('--keep', action='store_true', help='Keep the intermediate dot-format files')
    parse_draw.add_argument('-s', '--fontsize', action='store', help='The font size in the output figure. '
                                                                     'Default - 40', default='40')
    parse_draw.add_argument('-w', '--linewidth', action='store', help='The line width in the output figure node names. '
                                                                      'Default - 30', default='30')

    parse_draw.add_argument('-c', '--colors', action='store', nargs='+',
                            help='The colors for the populations in the same order of the populations. '
                                 'Must be in the same number of populations. '
                                 'Should be in a HEX format (for example: #5e4fa2) or a color name. '
                                 'The color names appears in https://www.graphviz.org/doc/info/colors.html '
                                 'Default - depends on the number of populations')
    parse_draw.set_defaults(function=call_draw)

    return parser


def call_draw(args):
    """
    Calls the draw function
    :param args: The input arguments
    :return: None
    """
    import draw_trees
    draw_trees.draw_trees(args)


def call_poptree(args):
    """
    Calls the draw function
    :param args: The input arguments
    :return: None
    """
    import poptree
    poptree.poptree(args)


def call_mtree(args):
    """
    Calls the draw function
    :param args: The input arguments
    :return: None
    """
    import mtree
    mtree.mtree(args)


def call_mutations(args):
    """
    Calls the draw function
    :param args: The input arguments
    :return: None
    """
    import count_muts
    count_muts.count_obs_func(args)


def call_filter(args):
    """
    Calls the draw function
    :param args: The input arguments
    :return: None
    """
    import filter
    filter.filter_rep(args)


if __name__ == "__main__":
    # Parse command line arguments
    gloabl_parser = get_arg_parser()
    arguments = gloabl_parser.parse_args()

    if check_arguments(arguments):
        arguments.function(arguments)
