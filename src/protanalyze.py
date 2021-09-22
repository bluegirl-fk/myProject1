# this is to analyze proteins of Human, Brain and NDD phenotypes creating arrays of intersection and unition of these.
import numpy as np
import pandas as pd
import config as cfg
import brain as brn
import varAnalysis as var


def phens_acc_dict_maker(phen_lst, phens_df):  # TODO: later try to use these two methods inside eachother
    phen_dict = dict.fromkeys(phen_lst)  # this dict has phen_name as key and lst of corresponding ACCs as values
    for each_phen in phen_lst:
        tmp_lst = list(set(phens_df.loc[phens_df['Phenotype'] == each_phen, 'acc'].tolist()))
        phen_dict[each_phen] = tmp_lst
    return phen_dict


def human_brain_acc_adder(phens_dict, vartype):
    mobidb, brain, _ = var.all_vars_or_vars_inidr(vartype)  # this creates dfs from varanalysis class
    # (all vars or only vars that are in IDR)
    human_lst = list(set(mobidb['acc'].tolist()))
    brain_lst = list(set(brain['acc'].tolist()))
    phens_dict['Human'] = human_lst
    phens_dict['Brain'] = brain_lst
    return phens_dict


def phens_intersect_df_maker(phen_acc_dic):
    # gets dictionary with phen as key and their uniprot IDs as values,
    intersection_df = pd.DataFrame(index=phen_acc_dic.keys(), columns=phen_acc_dic.keys())
    for key1 in phen_acc_dic.keys():
        intersection_count_lst = []
        for key2 in phen_acc_dic.keys():
            count_tmp = len(list(set(phen_acc_dic[key1]) & set(phen_acc_dic[key2])))  # this lst contains the ACCs
            # can achieve the mutual ACCs by storing the list in the upper line in a dict of dics or df cell
            intersection_count_lst.append(count_tmp)
        intersection_df[key1] = intersection_count_lst
    return intersection_df


def phens_union_df_maker(phen_acc_dic):
    union_df = pd.DataFrame(index=phen_acc_dic.keys(), columns=phen_acc_dic.keys())
    for key1 in phen_acc_dic.keys():
        union_count_lst = []
        for key2 in phen_acc_dic.keys():
            count_tmp = len(list(set(phen_acc_dic[key1]) | set(phen_acc_dic[key2])))
            union_count_lst.append(count_tmp)
        union_df[key1] = union_count_lst
    return union_df


if __name__ == '__main__':
    phens_lst = ['Human', 'Brain', 'ASD', 'EE', 'ID', 'DD', 'SCZ', 'NDDs', 'Control']
    _, _, ndd_subdf = var.all_vars_or_vars_inidr('all')
    ndd_subdf = ndd_subdf.drop_duplicates()  #
    # Dictionary with phen as key and their corresponding  list of ACCs as value
    phens_acc_dict_all = human_brain_acc_adder(phens_acc_dict_maker(phens_lst, ndd_subdf), 'all')
    phens_acc_dict_idr = human_brain_acc_adder(phens_acc_dict_maker(phens_lst, ndd_subdf), 'idr')
    # todo: change this class's code to have all and idr vars at the same time, purpose = to be divided by each other

    ## Intersection and Union
    phens_inters_all_df = phens_intersect_df_maker(phens_acc_dict_all)
    phens_inters_idr_df = phens_intersect_df_maker(phens_acc_dict_idr)
    # phens_inters_df.style.bar(subset=["ASD"], color='#FFA07A')
    # https://github.com/dexplo/dataframe_image very useful for saving pandas styler figures or jupyterlab to pdf
    phens_union_all_df = phens_union_df_maker(phens_acc_dict_all)
    phens_union_idr_df = phens_union_df_maker(phens_acc_dict_idr)

    ## percentage of vars that are occuring in IDR regions
    # divided and * 100