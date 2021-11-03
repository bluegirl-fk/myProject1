# this is to perform previous analysis which was done in mobibool and protanalyse but for 'variants in IDR' dataset
# Sep 17th
from itertools import count

import numpy as np
import pandas as pd
from pandas import DataFrame
import config as cfg
import brain as brn
import seaborn as sns
import matplotlib.pyplot as plt


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
    all_vars = pd.read_csv(cfg.data['vars'] + '/all-vars-mrg-mobidb-with-isin_idr-column.csv')
    all_vars = all_vars.loc[all_vars['feature'] == 'prediction-disorder-th_50']
    vars_in_idr_df = all_vars.loc[all_vars['isin_idr'] == 1]
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
    sum_of_normalized_vars = df['in_idr_vars_perc'] + df['out_idr_vars_perc']
    df['in_idr_vars_perc'] = df['in_idr_vars_perc'] / sum_of_normalized_vars
    df['out_idr_vars_perc'] = df['out_idr_vars_perc'] / sum_of_normalized_vars
    return df


def var_residue_stats_table_generator(input_vardf):
    # the values of number of vars in idr and out idr, number of residues for vars inidr and outidr generated,
    # converted to dict and then data frame, it was kept in mind that even though several vars could exist per protein,
    # residues of disorder or order regions are counted once
    vars_in_count = len(input_vardf.loc[input_vardf['isin_idr'] == 1])
    vars_out_count = len(input_vardf.loc[input_vardf['isin_idr'] == 0])

    inidr_vars_residue = input_vardf.loc[input_vardf['isin_idr'] == 1, ['acc', 'content_count']]
    inidr_vars_residue = inidr_vars_residue.drop_duplicates(subset='acc')
    inidr_vars_residue = int(inidr_vars_residue['content_count'].sum())

    outidr_vars_residue = input_vardf.loc[input_vardf['isin_idr'] == 0, ['acc', 'length', 'content_count']]
    outidr_vars_residue = outidr_vars_residue.drop_duplicates(subset='acc')
    outidr_vars_residue['order_res'] = (outidr_vars_residue['length'] - outidr_vars_residue['content_count'])
    outidr_vars_residue = int(outidr_vars_residue['order_res'].sum())
    # length_df = input_vardf[['acc', 'length']]
    # length = length_df['length'].sum()
    vars_count_dict = {'vars_count': [vars_in_count, vars_out_count, (vars_in_count + vars_out_count)],
                       'residues': [inidr_vars_residue, outidr_vars_residue,
                                    (inidr_vars_residue + outidr_vars_residue)]}
    vars_count_df = pd.DataFrame.from_dict(vars_count_dict, orient='index', columns=['disordered', 'ordered', 'sum'])
    return vars_count_df


def draw_barplot(x, y, data, xticklabel, yscale):  # input is DF, not list
    plt.figure(figsize=(12, 6))  # bigger figsize to have xticklabels shown
    sns.set_style("ticks")
    g = sns.barplot(x=x, y=y, data=data)
    # sns.despine(trim=True, offset=2)
    g.set_xticklabels(xticklabel, rotation=0, va="center", position=(0, -0.02))
    sns.color_palette("pastel")
    plt.yscale(yscale)
    plt.tight_layout()
    # plt.savefig(save_rout)
    plt.show()
    return


def residue_heatmapper(df, hmap_title, filename):
    # input df is disorder_majority or filtered_dis_maj if preferred
    res_df = df[['acc', 'var_id', 'orig_aa', 'var_aa']]
    aa_lst = res_df['orig_aa'].unique().tolist()
    res_df = res_df.loc[res_df.var_aa.isin(aa_lst)]
    residue_ser = res_df.groupby(['orig_aa', 'var_aa']).size().to_frame(name='size').reset_index()
    residue_pivot_df = pd.pivot_table(residue_ser, values='size', index='orig_aa', columns='var_aa')
    sns.set()
    fig = plt.figure(figsize=(20, 15))
    ax = plt.axes()
    ax.set_title('lalala')
    sns.heatmap(residue_pivot_df, linewidths=.5, cmap='viridis', vmin=0, vmax=1000, annot=True, fmt='g', ax=ax)
    ax.set_title(hmap_title, fontsize=30)
    plt.xlabel('Variant residues', fontsize=20)
    plt.ylabel('Original residues', fontsize=20)
    plt.savefig(cfg.plots['var-hms'] + '/' + filename + '.png', dpi=120)
    plt.show()
    # plt.savefig(cfg.plots['var-hms'] + '/' + filename + '.png', dpi=120)
    return residue_pivot_df


## NDD and brain original import
ndd_subdf = pd.read_csv(cfg.data['phens-fdr'] + '/acc-phen-5percentFDR.csv')
ndd_subdf = ndd_subdf.drop_duplicates()  # (4531, 3)
ndd_pr_lst = ndd_subdf['acc'].unique().tolist()  # 1308 proteins
brain_prot_lst = brn.brain_pr_lst_generator()  # n: 8320
## Mobidb disorder
# disorder_majority = var_cnt_residue_normaliezer(var_countcol_creator(df_feature_filterer('prediction-disorder-th_50')))
# # deleting proteins with content count of less than 20 residues
# disorder_majority = disorder_majority.loc[disorder_majority['content_count'] >= 20]
# disorder_majority.to_csv(cfg.data['vars'] + '/disorder-majority-inout-idr-vars-count-normalized.csv')
disorder_majority = pd.read_csv(cfg.data['vars'] + '/disorder-majority-inout-idr-vars-count-normalized.csv')

## Table of Vars inside and outside plus residues
var_residue_sum_table = var_residue_stats_table_generator(disorder_majority)
# def total_var_res_percentage(df): # to generate percentage for the whole dataset

# 4631 unique ACCs
dis_maj_filtered = disorder_majority.loc[disorder_majority['in_idr_vars'] >= disorder_majority['out_idr_vars']]
print(dis_maj_filtered['orig_aa'].value_counts()[:5].index.tolist())
print(dis_maj_filtered['var_aa'].value_counts()[:5].index.tolist())
orig_aa_count = dis_maj_filtered.groupby('orig_aa').count()
orig_aa_count = orig_aa_count.reset_index()
var_aa_count = dis_maj_filtered.groupby('var_aa').count()
var_aa_count = var_aa_count.reset_index()
aa_lst = disorder_majority['orig_aa'].unique().tolist()
var_aa_count = var_aa_count.loc[var_aa_count.var_aa.isin(aa_lst)]

# Distribution of altered residues, what about the new residues?
# turn this into a method and change y columns name to count or something more relevant
# draw_barplot(x='orig_aa', y='acc', data=orig_aa_count, xticklabel=orig_aa_count['orig_aa'].tolist(), yscale='linear')
# draw_barplot(x='var_aa', y='acc', data=var_aa_count, xticklabel=var_aa_count['var_aa'].tolist(), yscale='linear')

residue_heatmapper(disorder_majority, 'Residue Variations', 'Res-var-disorder_majority')
# residue_heatmapper(dis_maj_filtered)
