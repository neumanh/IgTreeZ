#!/usr/bin/python3
"""
Functions for linking tree nodes and sequence alignments
"""

# Info

__author__ = 'Hadas Neuman'
__version__ = '1.0'
__date__ = '13/10/20'

import multiprocessing as mp

import pandas as pd
from ete3 import Tree
from objects import Logger
from utils import general_utils, tree_utils

# Global variables
_seq_field = ''
_gl_field = ''
_clone_field = ''
_id_field = ''
_sample_field = ''
_cdr3_len_field = ''

# Germlines optional names
gls = ["GL", "G.L", 'G.L.', "Germline", "Germline_UCA", 'EX;G.L', 'EX_G.L']
current_log_object = None


def initiate_db_fileds_by_changeo():
    """
    Initiates the field names by Change-O convention
    :return: None
    """
    global _clone_field, _seq_field, _gl_field, _id_field, _cdr3_len_field
    _seq_field = 'SEQUENCE_IMGT'
    _gl_field = 'GERMLINE_IMGT_D_MASK'
    _clone_field = 'CLONE'
    _id_field = 'SEQUENCE_ID'
    _cdr3_len_field = 'CDR3_IMGT'


def initiate_db_fileds_by_airr():
    """
    Initiates the field names by the AIRR convention
    :return: None
    """
    global _clone_field, _seq_field, _gl_field, _id_field, _cdr3_len_field
    _seq_field = 'sequence_alignment'
    _gl_field = 'germline_alignment'
    _clone_field = 'clone_id'
    _id_field = 'sequence_id'
    _cdr3_len_field = 'cdr3'


def create_input_db(args, log_object: Logger):
    """
    Creates the sequences db
    :param args: The input arguments
    :param log_object: The Logger object
    :return: The dataframe
    """
    global _clone_field, _seq_field, _gl_field, _id_field

    # Initiating the field variables
    if args.dbformat == 'airr':
        initiate_db_fileds_by_airr()
    else:  # ChangeoO format
        initiate_db_fileds_by_changeo()

    # Update the field names by input arguments
    if args.clone_field:
        _clone_field = args.clone_field
    if args.seq_field:
        _seq_field = args.seq_field
    if args.gl_field:
        _gl_field = args.gl_field
    if args.id_field:
        _id_field = args.id_field

    extra_fields = [_cdr3_len_field]
    min_fields = [_id_field, _clone_field, _seq_field, _gl_field]
    max_fields = min_fields + extra_fields

    # db_file = args.db
    db_file = args.database
    df = pd.DataFrame()  # An empty dataframe

    if not general_utils.test_file_existence([db_file]):
        try:
            df = pd.read_csv(db_file, sep='\t', usecols=max_fields)
        except Exception:
            log_object.warning(f'Could not find some of the columns {extra_fields} in the database {db_file}.'
                               f' taking the other columns.')  # TEMP
            try:
                df = pd.read_csv(db_file, sep='\t', usecols=min_fields)
            except Exception as e:
                general_utils.handle_errors(f'Could not read the database file {db_file}: {e}. Exiting.', log_object,
                                            exit_stat=True)
    else:
        general_utils.handle_errors(f'The file {db_file} does not exist. Exiting.', log_object,
                                    exit_stat=True)

    return df


def attach_gl_to_tree(t: Tree, gl_seq: str):
    """
    Attaches the germline sequences to nodes
    :param t: a Tree
    :param gl_seq: The germline sequence
    :return: The updated tree
    """
    for node in t.traverse("preorder"):
        if not hasattr(node, 'sequence'):  # The node is not linked to any sequence
            if hasattr(node, 'name') and (node.name in gls):
                if gl_seq:
                    node.sequence = gl_seq
                    # print('Attaching gl to the node', node.name) # DEBUG

    return t


def link_alignment_to_tree_by_df(t: Tree, df: pd.DataFrame, log_object=current_log_object):
    """
    Links sequences to tree
    :param df: The partial dataframe, containing only the relevant clone
    :param t: A ETE tree object
    :param log_object: The log object
    :return: The updated tree, with the sequence attached
    """
    if not log_object:
        log_object = general_utils.create_general_log_object()

    if t:
        t = general_utils.replace_new_line_in_ig_tree(t)
        all_seq_assigned_flag = True

        # replacing the : with _ in the column name
        df[_id_field] = df[_id_field].str.replace(':', '-')

        # print(t.get_ascii(show_internal=True))

        # Attaching the germline sequences to the tree
        gl_seq, cdr3_len, first_seq = find_gl_seq_in_df(df, general_utils.get_clone_number(t.id),
                                                        get_first_none_root_node_name(t))

        # gl_seq = find_gl_seq_in_df(df, general_utils.get_clone_number(t.id))
        seq_num = 0
        if gl_seq:
            gl_seq = gl_seq.replace('.', '-')  # Replacing the dot with dashes
            # t = attach_gl_to_tree(t, gl_seq)
            t.sequence = gl_seq.upper()

            root = t

            # Attaching all the other sequences to the tree
            for node in t.traverse("preorder"):  # Going from root to leaves
                if (node != root) and (hasattr(node, 'name')):  # This is not the root
                    seq = find_seq_in_df(node.name, df)
                    if seq:
                        node.sequence = seq.replace('.', '-').upper()  # Replacing IMGT's '.' with standart '-'
                        seq_num += 1

                else:
                    all_seq_assigned_flag = False

            # Attaching the hypothetical nodes
            alright_flag = True
            if hasattr(t, 'sequence') and (not all_seq_assigned_flag):

                for node in t.traverse(
                        'preorder'):  # Going from root to leaves to fill all the unknown sequences by father
                    if hasattr(node, 'sequence'):
                        node.sequence = node.sequence.replace('.', '-').upper()
                    # If the node is identical to its father
                    elif (node.get_distance(node.up) == 0) and (hasattr(node.up, 'sequence')):
                        node.sequence = node.up.sequence.replace('.', '-').upper()

                for node in t.traverse("postorder"):  # Going from leaves to root to generate missing sequences
                    if alright_flag:
                        if not hasattr(node, 'sequence'):  # this in a hypothetical node w/o sequence
                            generated_seq = tree_utils.generate_missing_sequence(node.children, node.up, log_object)
                            if generated_seq:
                                node.sequence = generated_seq
                            else:
                                general_utils.handle_errors(f'Could not generate a sequence to the node {node.name}.'
                                                            f'\nThe tree will not be analysed.',
                                                            log_object=log_object, exit_stat=False)
                                alright_flag = False
        else:
            general_utils.handle_errors(f'Could not find a germline sequence for tree {t.id}',
                                        log_object, exit_stat=False)
            alright_flag = False
        t.seq_num = seq_num

        if cdr3_len:
            #cdr3_end = get_adjusted_cdr3_end(general_utils.imgt_regions['cdr3'], cdr3_len, first_seq)
            cdr3_end = general_utils.imgt_regions['cdr3'] + cdr3_len
            t.cdr3_end = cdr3_end

        if not alright_flag:
            t = None
    else:  # Could not find germline sequence
        t = None

    return t


def get_first_none_root_node_name(t):
    """
    Gets the first non-root node name
    :param t: An ete3 tree
    :return: The first node name
    """
    name = None
    root = t
    for node in t.traverse():
        if (node != root) and (not name):
            if hasattr(node, 'name'):
                name = node.name
    return name


def find_gl_seq_in_df(df: pd.DataFrame, clone_number: int = None, first_seq_name: str = None):
    """
    Attaches the germline sequences to nodes
    :param df: The dataframe with all sequences
    :param clone_number: The clone number
    :param first_seq_name: The sequence name of the first node
    :return: The germline sequence
    """

    gl_seq, cdr3_length, first_seq = None, None, None
    gl_df = pd.DataFrame()  # An empty dataframe

    clones = df[_clone_field].unique()  # Number of different clones
    # print(clones)

    if not df.empty:
        if len(clones) == 1:  # The database contains only one clone
            gl_df = df
        else:  # The database contains several clones
            if first_seq_name:  # Find the GL by the sequence name
                gl_df = df[df[_id_field] == first_seq_name].dropna()
            if gl_df.empty:
                if clone_number:
                    # Look for the GL_FIELD values where the clone number is an the input, and take the first value
                    gl_df = df[df[_clone_field] == clone_number].dropna()  # TODO: NEED TO BE TESTED!!!

        if (not gl_df.empty) and (gl_df.shape[0] > 0):
            # Get the GL of the first line (using values because the index is not permanent)
            gl_seq = gl_df[_gl_field].values[0]
            if _cdr3_len_field in gl_df:
                cdr3_length = len(gl_df[_cdr3_len_field].values[0])
                first_seq = gl_df[_seq_field].values[0]

    return gl_seq, cdr3_length, first_seq


def find_seq_in_df(seq_name, df):
    """
    Find the sequence by it's name in a df
    :param seq_name: The sequence name
    :param df: The partial dataframe, for a specific clone
    :return: The sequence
    """

    seq = None
    # temp_df = df[SEQ_FIELD].where(df[ID_FIELD] == seq_name).dropna() # 25 % Time
    temp_df = df[df[_id_field] == seq_name].dropna()  # TODO: NEED TO BE TESTED!!!

    if not temp_df.empty:
        seq = temp_df[_seq_field].iloc[0]  # TODO: TEST
    else:
        for row in zip(df[_id_field], df[_seq_field]):
            if row[0] in seq_name:  # row[0] - df[ID_FIELD]
                seq = row[1]  # row[1] - df[SEQ_FIELD]

    return seq


def generate_trees_and_clone_dfs(tree_list: list, df: pd.DataFrame, log_object: Logger = None):
    """
    Creates the tree lists and dataframe lists
    :param tree_list: The input tree list
    :param df: The input dataframe
    :param log_object: The Logger object
    :return: Two lists - ete tree list, and pandas df list
    """
    if not log_object:
        log_object = general_utils.create_general_log_object()

    for t in tree_list:
        clone_num = general_utils.get_clone_number(t.id)
        clone_df = create_clone_df(df, clone_num)

        yield t, clone_df, log_object


def create_clone_df(full_df: pd.DataFrame, clone_number: int):
    """
    Creates temporary dataframe of a specific clone
    :param full_df: The all sequences df
    :param clone_number: the clone number
    :return: The partial dataframe
    """
    temp_df = full_df[full_df[_clone_field] == clone_number].dropna()
    temp_df = temp_df.reset_index(drop=True)  # Resetting the index
    if temp_df.empty:
        temp_df = full_df

    return temp_df


def attach_seqs_by_df(args, log_object: Logger = None):
    """
    Collect trees and attaches sequences to tree nodes.
    :param args: The input arguments
    :param log_object: The Logger object
    :return: A list of linked trees
    """
    global current_log_object

    if not log_object:
        current_log_object = Logger(args.name, args.silent)

    tree_list = general_utils.collect_trees(args, log_object)

    # Creates the relevant database
    df = create_input_db(args, log_object)

    if not df.empty:  # The returned dataframe is not empty
        # Link alignment to nodes in parallel
        args.ncors = int(args.ncors)
        number_of_workers = min(args.ncors, len(tree_list))
        with mp.Pool(number_of_workers) as pool:
            # Using the generate_trees_and_clone_dfs generator
            trees = pool.starmap(link_alignment_to_tree_by_df,
                                 generate_trees_and_clone_dfs(tree_list, df, log_object))
    else:
        general_utils.handle_errors(f'Could not create a dataframe from the file {args.database}', log_object,
                                    exit_stat=True)
        trees = []

    return trees


def update_cdr3_length(t: Tree, cdr3_length):
    pass


def get_adjusted_cdr3_end(cdr3_start_pos: int, cdr3_length: int, seq: str):
    """
    Adjust the cdr3 length by the gaped sequence.
    :param cdr3_start_pos: The CDR3 coordinate
    :param cdr3_length: The CDR3 length (without gapes)
    :param seq: The aligned sequence
    :return: The adjusted CDR3 end position
    """
    short_seq = seq[cdr3_start_pos:]
    idx, non_gap_idx = 0, 0
    while non_gap_idx < cdr3_length:
        if short_seq[idx] != '-':
            non_gap_idx += 1
        idx += 1
    cdr3_end = cdr3_start_pos + idx
    return cdr3_end
