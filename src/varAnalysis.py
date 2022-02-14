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
import collections
from itertools import product
import humsavar as hvar


def vars_multiple_df_generator(input_df):
    # gets df as input (either merged mobidb with all variants or the variants in IDR), then creates ndd and brain subdf
    # based on proteins in the var dataset, basically deletes ACCs from ndd and brain that don't exist in the human var
    # list, then these dfs will be used in mobibool or protanalysis
    var_pr_lst = input_df['acc'].unique().tolist()
    ## ndds in variants
    ndd_subdf = pd.read_csv(cfg.data['phens-fdr'] + '/acc-phen-5percentFDR.csv')
    ndd_pr_lst = ndd_subdf['acc'].unique().tolist()  # 1308 proteins
    ndd_in_var = list(set(ndd_pr_lst).intersection(var_pr_lst))
    ndd_in_var_subdf = ndd_subdf[ndd_subdf.acc.isin(ndd_in_var)]
    # del ndd_in_var_subdf['Unnamed: 0']
    ## Brain
    brain_prot_lst = brn.brain_pr_lst_generator()  # n: 8320
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
    ## this where method conditions is to prevent nan values for fully disordered proteins (length-content_count == 0)
    # in other proteins cases it just divides number of vars on number of residues
    df['in_idr_vars_perc'] = np.where((df['length'] == df['content_count']), 1,
                                      (df['in_idr_vars'] / df['content_count']))
    df['out_idr_vars_perc'] = np.where((df['length'] == df['content_count']), 0,
                                       (df['out_idr_vars'] / (df['length'] - df['content_count'])))

    # df['in_idr_vars_perc'] = (df['in_idr_vars'] / df['content_count'])
    # df['out_idr_vars_perc'] = (df['out_idr_vars'] / (df['length'] - df['content_count']))
    ## now we will divide each by sum of the calculated data for in+out IDRs to produce complementary % values for cols
    sum_of_normalized_vars = df['in_idr_vars_perc'] + df['out_idr_vars_perc']
    df['in_idr_vars_perc'] = df['in_idr_vars_perc'] / sum_of_normalized_vars
    df['out_idr_vars_perc'] = df['out_idr_vars_perc'] / sum_of_normalized_vars
    return df


def var_residue_stats_table_generator(input_vardf, inp_var):
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
    vars_count_dict = {inp_var+' variants count': [vars_in_count, vars_out_count, (vars_in_count + vars_out_count)],
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
    # g.bar_label(g.containers[0])
    plt.yscale(yscale)
    plt.tight_layout()
    # plt.savefig(save_rout)
    plt.show()
    return


def residue_heatmapper(df_lst, hmap_title_lst, filename):
    # input df is disorder_majority or filtered_dis_maj if preferred
    new_pivot_df_lst, aa_symbols_lst = heatmap_pivotdf_maker(df_lst)
    sns.set(font_scale=2.4)
    fig, axes = plt.subplots(len(new_pivot_df_lst), 1, figsize=(30, 20 * len(new_pivot_df_lst)))
    for i, (ax, d, t) in enumerate(zip(axes.reshape(-1), new_pivot_df_lst, hmap_title_lst)):
        sb = sns.heatmap(d, cmap="viridis", annot=True, fmt='g', linewidth=1.1, linecolor='black', ax=ax, square=True,
                         cbar_kws={'label': 'Transition percentage'}
                         )
        ax.set_title(t, fontsize=40)
        ax.set_xlabel('Variant residues', fontsize=30)
        ax.set_ylabel('Original residues', fontsize=30)
        ax.tick_params(axis='x', colors='red')
        ax.set_xticklabels(aa_symbols_lst)

        # if i < (len(new_pivot_df_lst) - 1):
        #     sb.set(xticklabels=[])
        #     sb.set(xlabel=None)
    plt.tight_layout()
    plt.savefig(cfg.plots['var-hms'] + '/' + filename + '.png', dpi=120)
    plt.show()

    # ax.set_title(hmap_title, fontsize=40)
    # plt.xlabel('Variant residues', fontsize=25)
    # # ax.xaxis.label.set_color('red')
    # plt.ylabel('Original residues', fontsize=25)
    # # ax.yaxis.label.set_color('blue')
    return


def heatmap_pivotdf_maker(df_lst):
    # todo check if the transitions are counted correctly, because the numbers here are different
    new_pivot_df_lst = []
    for df in df_lst:
        res_df = df[['acc', 'var_id', 'orig_aa', 'var_aa']]
        aa_categ_order_x = ['C', 'M', 'V', 'L', 'I', 'W', 'Y', 'F', 'H', 'T', 'N', 'Q', 'K', 'R', 'D', 'E', 'A', 'S',
                            'G',
                            'P']
        res_df = res_df.loc[res_df.var_aa.isin(aa_categ_order_x)]
        residue_ser = res_df.groupby(['orig_aa', 'var_aa']).size().to_frame(name='size').reset_index()
        residue_pivot_df = pd.pivot_table(residue_ser, values='size', index='orig_aa', columns='var_aa')
        aa_categ_order_y = aa_categ_order_x.__reversed__()
        residue_pivot_df = residue_pivot_df.reindex(columns=aa_categ_order_x)
        residue_pivot_df = residue_pivot_df.reindex(aa_categ_order_y)
        # max_value = residue_pivot_df.max() # this finds max of each column and the result percentage is comparative
        # only for each column
        max_value = residue_ser['size'].max()
        residue_pivot_df = (residue_pivot_df.div(max_value)).round(3) * 100
        residue_pivot_df = residue_pivot_df.apply(pd.to_numeric)
        new_pivot_df_lst.append(residue_pivot_df)
    difference_df = new_pivot_df_lst[0].subtract(new_pivot_df_lst[1])
    new_pivot_df_lst.append(difference_df)
    return new_pivot_df_lst, aa_categ_order_x


if __name__ == '__main__':
    ## NDD and brain original import
    ndd_subdf = pd.read_csv(cfg.data['phens-fdr'] + '/acc-phen-5percentFDR.csv')
    ndd_subdf = ndd_subdf.drop_duplicates()  # (4531, 3)
    ndd_pr_lst = ndd_subdf['acc'].unique().tolist()  # 1308 proteins
    brain_prot_lst = brn.brain_pr_lst_generator()  # n: 8320
    ## Mobidb disorder
    # mobidb_lite = var_cnt_residue_normaliezer(var_countcol_creator(df_feature_filterer('prediction-disorder-mobidb_lite')))
    # mobidb_lite.to_csv(cfg.data['vars'] + '/mobidb_lite-inout-idr-vars-count-normalized.csv')
    mobidb_lite = pd.read_csv(cfg.data['vars'] + '/mobidb_lite-inout-idr-vars-count-normalized.csv')
    disorder_majority = pd.read_csv(cfg.data['vars'] + '/disorder-majority-inout-idr-vars-count-normalized.csv')
    ndd_mobilite = mobidb_lite.loc[mobidb_lite.acc.isin(ndd_pr_lst)]
    ndd_dismaj = disorder_majority.loc[disorder_majority.acc.isin(ndd_pr_lst)]

    ## Table of Vars inside and outside plus residues
    var_residue_sum_table_all = var_residue_stats_table_generator(disorder_majority, 'all')
    var_residue_sum_table_ndd = var_residue_stats_table_generator(ndd_dismaj, 'NDD')
    # def total_var_res_percentage(df): # to generate percentage for the whole dataset

    # 4631 unique ACCs
    dismaj_vars_in = disorder_majority.loc[disorder_majority['isin_idr'] == 1]  # (8120, 18)
    ndd_dismaj_vars_in = dismaj_vars_in.loc[dismaj_vars_in.acc.isin(ndd_pr_lst)]

    # mobilite_vars_out = mobidb_lite.loc[mobidb_lite['isin_idr'] == 0]  # (37672, 18)
    # ndd_mobilite_vars_out = mobilite_vars_out.loc[mobilite_vars_out.acc.isin(ndd_pr_lst)]

    ## Variants description (in disorder)
    ndd_vars_in_descriptions = ndd_dismaj_vars_in['description'].unique().tolist()
    desc_counter = collections.Counter(ndd_vars_in_descriptions)
    # from: https://stackoverflow.com/questions/2161752/how-to-count-the-frequency-of-the-elements-in-an-unordered-list
    # parse uniprot-xml keyword disease part
    print(desc_counter)  # later you can somehow categorize them based on type of pathology and stuff

    print(dismaj_vars_in['orig_aa'].value_counts()[:5].index.tolist())
    # most_occuring_aa_transition_inIDR = mobilite_vars_in.groupby(['orig_aa', 'var_aa']).size().idxmax(5)

    ## heatmaps
    residue_heatmapper([ndd_dismaj_vars_in, dismaj_vars_in], ['Residue Variations - in IDRs (NDDs)',
                                                                  'Residue Variations- in IDRs (Homo sapiens)',
                                                                  'Difference (NDD - Homo sapiens)'],
                       'heatmap-inidr-HS&NDD-dismaj-fb1')

