#!/usr/bin/python3
"""
Tests statistic functions
"""
import pandas as pd
from objects import Logger
from scipy import stats

alpha = 0.05


def shapiro(data_list: list):
    """
    Runs the Shapiro test on all data series
    :param data_list: A list of series
    :return: The smallest shapiro test results
    """
    sh_results = []
    for ser in data_list:
        sh = stats.shapiro(ser)[1]  # If pval > 0.05 => noramlly dist.
        sh_results.append(sh)

    min_sh = min(sh_results)

    return min_sh


def levene(data_list: list):
    """
    Calculates the levene test of homogeneity of variances
    :param data_list: A list of series
    :return: The tests result
    """
    lev = stats.levene(*data_list)

    return lev[1]


def is_paramteric(data_list: list):
    """
    Tests if the data is parametric by testing noramllity and homogeneity of variances
    :param data_list: A list of series
    :return: True if the data is normal & with homogeny variances, or False otherwise
    """
    out_str = ''
    parametric = False
    if shapiro(data_list) > alpha:  # If pval > 0.05 => normally dist.
        if levene(data_list) > alpha:  # If pval > 0.05 => homogenous variances
            parametric = True
            out_str += 'The data is normally distributed and has a homogeneity of variances'
        else:
            out_str += 'The assumption of homogeneity of variances is not met.\n'
    else:
        out_str += 'The normality assumption is not met.\n'

    return parametric, out_str


def benj_hoch(pval_dict: dict, out_file: str):
    """
    Calculates the Benjamin-Hochberg test
    :param pval_dict: A doctionary of label:p-value
    :param out_file: The output file
    :return: The adjusted P values
    """
    q = alpha  # The false discovery rate
    pval_tuples = sorted(pval_dict.items(), key=lambda kv: kv[1])
    lables, pvals = (zip(*pval_tuples))
    pval_num = len(pvals)
    bh_dict = {}
    rank = 1
    for _ in pvals:
        bh = q * (rank / pval_num)
        bh_dict[lables[rank - 1]] = bh
        rank += 1
    str_out = f'Benjamin-Hochberg false detection rate results: \n{bh_dict}'
    with open(out_file, 'a') as f:
        f.write(str_out)


def ttest(arg_list):
    """
    Performs T-test for each pair of data
    :param arg_list: The lists to compare
    :return: A list of all T-test
    """
    out_str = ''
    pval_dict = {}
    for idx1, arg1 in enumerate(arg_list):
        for idx2, arg2 in enumerate(arg_list):
            if idx1 < idx2:
                tt = stats.ttest_ind(arg1, arg2, nan_policy='omit')
                out_str += f'comparing {arg1.name} with {arg2.name} using T-test:\t {tt}\n'
                dict_label = f'{arg1.name}_to_{arg2.name}'
                pval_dict[dict_label] = tt[1]

    return out_str, pval_dict


def mwtest(arg_list):
    """
    Performs T-test for each pair of data
    :param arg_list: The lists to compare
    :return: A list of all T-test
    """
    # arg_list = [args]
    # print(len(arg_list))
    out_str = ''
    pval_dict = {}
    for idx1, arg1 in enumerate(arg_list):
        for idx2, arg2 in enumerate(arg_list):
            if idx1 < idx2:
                mw = stats.mannwhitneyu(arg1, arg2)
                out_str += f'comparing {arg1.name} with {arg2.name} using Mann Whitney test:\t {mw}\n'
                dict_label = f'{arg1.name}_to_{arg2.name}'
                pval_dict[dict_label] = mw[1]

    return out_str, pval_dict


def all_identicle(data_list: list):
    """
    Checks if all values are identicle
    :param data_list: A list of series
    :return: True if all values are identicle
    """
    all_identicle_flag = True
    uniq_vals = data_list[0].unique()
    for data in data_list:
        if all_identicle_flag:
            curr_uniq_vals = data.unique()
            if (len(curr_uniq_vals) > 1) or (curr_uniq_vals[0] != uniq_vals[0]):
                all_identicle_flag = False
    return all_identicle_flag


def analyze_df(df: pd.DataFrame, name: str, out_file: str = None, log_object: Logger = None):
    """
    Run statistic analysis on Dataframe
    :param df: A Dataframe
    :param name: The comparison field
    :param out_file: The output file
    :param log_object: The Logger object, for documentation
    :return: The ANOVA/Kruskal Wallis p-value
    """
    non_na_num = len(df.dropna())
    if non_na_num > 5000:
        non_na_num = 5000
        log_object.warning('The significance analysis was performed on 5000 randomly selected trees.')

    do_no_analyze_fields = ['nodes', 'leaves']
    df = df.dropna()
    data_list = [df[col].sample(n=non_na_num) for col in df]  # A list of series
    pval = None
    pval_dict = None

    # name = data_list[0].name

    if not out_file:
        out_file = f'{name}_stats.txt'

    if name not in do_no_analyze_fields:
        out_str = f'Analyzing {name}\n\n'

        # Describe the data
        desc = str(df.describe())
        out_str += f'{desc}\n'

        if not all_identicle(data_list):
            parametric, param_str = is_paramteric(data_list)
            out_str += param_str

            if parametric:
                # ANOVA
                anova = stats.f_oneway(*data_list)
                pval = anova[1]

                out_str += f'ANOVA results: {anova}\n'

                # Post Hoc
                curr_str, pval_dict = ttest(data_list)
                out_str += curr_str
                # out_str += ttest(data_list)

            else:
                # Converting the series to lists
                data_list_list = []
                # print(data_list)
                for d in data_list:
                    data_list_list.append(d.to_list())

                # Kruskal-Wallis H-test
                kw = stats.mstats.kruskal(*data_list_list)
                pval = kw[1]

                out_str += f'Kruskal Wallis results: {kw}\n'

                # Post Hoc
                curr_str, pval_dict = mwtest(data_list)
                out_str += curr_str
                # out_str += mwtest(data_list)
        else:
            out_str += 'All values are identicle. Not analyzing siginificance.\n'

        out_str += '-------------------------------------------------------------\n\n'

        with open(out_file, 'a') as f:
            f.write(out_str)

    return pval, pval_dict
