# this is to perform previous analysis which was done in mobibool and protanalyse but for 'variants in IDR' dataset
# Sep 17th
from itertools import count

import numpy as np
import pandas as pd
from pandas import DataFrame
import config as cfg
import brain as brn
import seaborn as sns


def vars_multiple_df_generator(input_df):
    # gets df as input (either merged mobidb with all variants or the variants in IDR), then creates ndd and brain subdf
    # based on proteins in the var dataset, basically deletes ACCs from ndd and brain that don't exist in the human var
    # list, then these dfs will be used in mobibool or protanalysis
    var_pr_lst = input_df['acc'].unique().tolist()
    ## ndds in variants
    ndd_in_var = list(set(ndd_pr_lst).intersection(var_pr_lst))
    ndd_in_var_subdf = ndd_subdf[ndd_subdf.acc.isin(ndd_in_var)]
    # del ndd_in_var_subdf['Unnamed: 0']
    ## Brain
    brain_var_lst = list(set(brain_prot_lst).intersection(var_pr_lst))
    brain_var_subdf = DataFrame(brain_var_lst, columns=['acc'])

    return input_df, brain_var_subdf, ndd_in_var_subdf


def all_vars_or_vars_inidr(all_or_idr):
    # input "all" all vars are needed and input 'idr' when only variants inside idr are needed
    vars_in_idr_df = dis_maj_filtered.loc[dis_maj_filtered['isin_idr'] == 1]
    all_vars = pd.read_csv(cfg.data['vars'] + '/all-vars-mrg-mobidb-with-isin_idr-column.csv')
    all_vars = all_vars.loc[all_vars['feature'] == 'prediction-disorder-th_50']
    if all_or_idr == 'all':
        return vars_multiple_df_generator(all_vars)
    elif all_or_idr == 'idr':
        return vars_multiple_df_generator(vars_in_idr_df)


def df_feature_filterer(feature):  # this generates a filtered df with the desired mobidb feature filtered then it will
    # be the input of var_idr_percentage_creator() method, input is the mobidb feature we want variants for
    vars_all_features = pd.read_csv(cfg.data['vars'] + '/all-vars-mrg-mobidb-with-isin_idr-column.csv')
    vars_my_feature = vars_all_features.loc[vars_all_features['feature'] == feature]
    return vars_my_feature


def var_countcol_creator(df):  # a percentage column for variants in IDR / all Vars
    # input it the out put of df_feature_filtered method
    var_count_dic = (df.pivot_table(columns=['acc'], aggfunc='size')).to_dict()
    # dict values to list so that I could append other data to dict values
    for key, value in var_count_dic.items():
        var_count_dic[key] = [value]
    for k in var_count_dic:  # adds the count of each acc var being in idr
        each_acc_vars_in_idr_count = df.loc[df['acc'] == k, 'isin_idr'].sum()
        var_count_dic[k].append(each_acc_vars_in_idr_count)
    # dict to df
    var_count_df = pd.DataFrame.from_dict(var_count_dic, orient='index', columns=['total_vars', 'in_idr_vars'])
    var_count_df['out_idr_vars'] = (var_count_df['total_vars'] - var_count_df['in_idr_vars'])
    var_count_df = var_count_df.reset_index()
    var_count_df = var_count_df.rename(columns={'index': 'acc'})
    mrg_var_and_count_df = pd.merge(df, var_count_df, on='acc')
    mrg_var_and_count_df = mrg_var_and_count_df.drop(columns=['Unnamed: 0', 'Unnamed: 0_x', 'Unnamed: 0_y'])
    return mrg_var_and_count_df


def var_cnt_residue_normaliezer(df):  # gets df with count columns for vars in/out idr and normalize them based on
    # number of disordered residues or non disordered, respectively
    df['in_idr_vars_perc'] = (df['in_idr_vars'] / df['content_count'])
    df['out_idr_vars_perc'] = (df['out_idr_vars'] / (df['length'] - df['content_count']))
    ## now we will divide each by sum of the calculated data for in+out IDRs to produce complementary % values for cols
    sum_of_normalized_vars = df['in_idr_vars_perc']+df['out_idr_vars_perc']
    df['in_idr_vars_perc'] = df['in_idr_vars_perc']/sum_of_normalized_vars
    df['out_idr_vars_perc'] = df['out_idr_vars_perc']/sum_of_normalized_vars
    return df


## NDD and brain original import
ndd_subdf = pd.read_csv(cfg.data['phens-fdr'] + '/acc-phen-5percentFDR.csv')
ndd_subdf = ndd_subdf.drop_duplicates()  # (4531, 3)
ndd_pr_lst = ndd_subdf['acc'].unique().tolist()  # 1308 proteins
brain_prot_lst = brn.brain_pr_lst_generator()  # n: 8320
## Mobidb disorder
# disorder_majority = var_cnt_residue_normaliezer(var_countcol_creator(df_feature_filterer('prediction-disorder-th_50')))
# disorder_majority.to_csv(cfg.data['vars'] + '/disorder-majority-inout-idr-vars-count-normalized.csv')
disorder_majority = pd.read_csv(cfg.data['vars'] + '/disorder-majority-inout-idr-vars-count-normalized.csv')

# 4631 unique ACCs
dis_maj_filtered = disorder_majority.loc[disorder_majority['in_idr_vars'] >= disorder_majority['out_idr_vars']]


