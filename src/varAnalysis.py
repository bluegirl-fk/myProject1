# this is to perform previous analysis which was done in mobibool and protanalyse but for 'variants in IDR' dataset
# Sep 17th
from itertools import count

import numpy as np
import pandas as pd
from pandas import DataFrame
import config as cfg
import brain as brn


def vars_df_generator(input_df):
    # gets df as input (either merged mobidb with all variants or the variants in IDR), then creates ndd and brain subdf
    # based on proteins in the var dataset, basically deletes ACCs from ndd and brain that don't exist in the human var
    # list, then these dfs will be used in mobibool or protanalysis
    var_pr_lst = input_df['acc'].unique().tolist()
    ## ndds in variants
    ndd_in_var = list(set(ndd_pr_lst).intersection(var_pr_lst))
    ndd_in_var_subdf = ndd_subdf[ndd_subdf.acc.isin(ndd_in_var)]
    del ndd_in_var_subdf['Unnamed: 0']
    ## Brain
    brain_var_lst = list(set(brain_prot_lst).intersection(var_pr_lst))
    brain_var_subdf = DataFrame(brain_var_lst, columns=['acc'])

    return input_df, brain_var_subdf, ndd_in_var_subdf


def all_vars_or_vars_inidr(input):
    # input "all" all vars are needed and input 'idr' when only variants inside idr are needed
    vars_in_idr_df = pd.read_csv(cfg.data['vars'] + '/all-vars-all-features.csv')
    all_vars = pd.read_csv(cfg.data['vars'] + '/all-vars-mrg-mobidb-with-isin_idr-column.csv')
    del all_vars['Unnamed: 0_x']
    del all_vars['Unnamed: 0_y']
    ## filtering only the disorder features mobidb lite and majority
    feaures_lst = ['prediction-disorder-mobidb_lite', 'prediction-disorder-th_50']
    vars_in_idr_df = vars_in_idr_df[vars_in_idr_df.feature.isin(feaures_lst)]
    all_vars = all_vars[all_vars.feature.isin(feaures_lst)]
    if input == 'all':
        return vars_df_generator(all_vars)
    elif input == 'idr':
        return vars_df_generator(vars_in_idr_df)


def df_feature_filterer(feature):  # this generates a filtered df with the desired mobidb feature filtered then it will
    # be the input of var_idr_percentage_creator() method, input is the mobidb feature we want variants for
    vars_all_features = pd.read_csv(cfg.data['vars'] + '/all-vars-mrg-mobidb-with-isin_idr-column.csv')
    vars_my_feature = vars_all_features.loc[vars_all_features['feature'] == feature]
    return vars_my_feature


def var_idr_percentage_creator(df):  # a percentage column for variants in IDR / all Vars
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
    var_count_df['percentage'] = (var_count_df['in_idr_vars'] * 100) / var_count_df['total_vars']
    var_count_df = var_count_df.reset_index()
    var_count_df = var_count_df.rename(columns={'index': 'acc'})
    mrg_var_and_percentage_df = pd.merge(df, var_count_df, on='acc')
    mrg_var_and_percentage_df = mrg_var_and_percentage_df.drop(columns=['Unnamed: 0', 'Unnamed: 0_x', 'Unnamed: 0_y'])
    return mrg_var_and_percentage_df


## NDD and brain original import
ndd_subdf = pd.read_csv(cfg.data['phens-fdr'] + '/acc-phen-5percentFDR.csv')
ndd_subdf = ndd_subdf.drop_duplicates()  # (4531, 3)
ndd_pr_lst = ndd_subdf['acc'].unique().tolist()  # 1308 proteins
brain_prot_lst = brn.brain_pr_lst_generator()  # n: 8320

mobidb_lite_var_count_df = var_idr_percentage_creator(df_feature_filterer('prediction-disorder-mobidb_lite'))
mobidb_lite_var_count_df.to_csv(cfg.data['vars'] + '/mobidb-lite-idrvars-percentage.csv')
