# this file is supposed to get mobidb+mutposition, is_in_startend and ndd bool, index*, and give a df in the outcome
# with added columns: brain bool, ASD bool, EE bool, ID bool, SCZ bool
# *index (mutual idx of mut_acc_merged_df and mutinfo subdf which is useful to merge mobidb & g4dn_mutinfo_acc,
# could be useless now :/ not sure
## this description will be changed! (August 4th)

import pandas as pd
import numpy as np
from pandas import DataFrame
import config as cfg
import sys
import brain as bd
import matplotlib.pyplot as plt
from matplotlib import ticker as mticker
import seaborn as sns


def box_plotter(data, save_route):
    plt.figure(figsize=(60, 60))  # bigger figsize to have xticklabels shown
    g = sns.catplot(data=data, kind="box")
    sns.set_style("ticks")
    g.set_xticklabels(rotation=45, va="center", position=(0, -0.02))
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(save_route)
    plt.close()


def violin_plotter(data, save_route):
    plt.figure(figsize=(12, 12))
    sns.set_theme(style="whitegrid")
    g = sns.violinplot(data=data, split=True, bw=.1)
    sns.set_style("ticks")
    # g.set_xticklabels(labels=df['Phenotype'].unique().tolist(),rotation=45, va="center", position=(0, -0.02))
    # plt.yscale('log')
    g.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
    # g.yaxis.set_ticks([np.log10(x) for p in range(-4, 5) for x in np.linspace(10 ** p, 10 ** (p + 1), 10)],minor=True)
    plt.tight_layout()
    plt.savefig(save_route)
    plt.close()


if __name__ == '__main__':
    mobidb = pd.read_csv(cfg.data[''] + '/mobidb_result.tsv', sep='\t')
    # adding a human column for all the proteins in mobidb
    mobidb['human'] = 'Human'
    ## selected features
    features_lst = ['prediction-disorder-mobidb_lite', 'prediction-low_complexity-merge',
                    'prediction-lip-anchor', 'homology-domain-merge']
    ## brain proteins
    brain_prot_lst = bd.brain_pr_lst_generator()  # n: 8428
    brain_pr_df = DataFrame(brain_prot_lst, columns=['acc'])
    brain_pr_df['brain'] = 'Brain'
    mobidb = mobidb.merge(brain_pr_df, how='left', on='acc')

    ## NDD proteins, could also specify index_col= ..., and pass a list for multiple idxs
    phen_df = pd.read_csv(cfg.data['gene4'] + '/positive_cand_g4mobi_concat.csv', usecols=["acc_x", "Phenotype"])
    phen_df = phen_df.drop_duplicates()
    ## merge to get Phenotypes column
    mobidb = mobidb.merge(phen_df, how='left', left_on='acc', right_on='acc_x')
    mobidb = mobidb.drop(columns=['start..end', 'acc_x'])

    ## content_fraction
    # multi-level index Phenotypes: from: https://www.youtube.com/watch?v=tcRGa2soc-c
    mobi_feature_df = mobidb[mobidb.feature.isin(features_lst)]  # (194383, 7)
    cf_phens_series = mobi_feature_df.groupby(['acc', 'feature', 'Phenotype']).content_fraction.mean()
    cf_phen_3d_df = cf_phens_series.unstack().sort_index()
    # multi-idx brain to be merged with multi-id phen df
    mobi_br_subdf = mobidb[['acc', 'feature', 'content_fraction', 'content_count', 'length', 'brain']]
    mobi_br_subdf = mobi_br_subdf[mobi_br_subdf.feature.isin(features_lst)]
    cf_brain_3d_df = mobi_br_subdf.groupby(['acc', 'feature', 'brain']).content_fraction.mean().unstack().sort_index()
    # multi-idx add human
    mobi_human_subdf = mobidb[['acc', 'feature', 'content_fraction', 'content_count', 'length', 'human']]
    mobi_human_subdf = mobi_human_subdf[mobi_human_subdf.feature.isin(features_lst)]
    cf_human_3d_df = mobi_human_subdf.groupby(['acc', 'feature', 'human']).content_fraction.mean().unstack().sort_index()
    ## merge brain and phens+ merge that one with human
    mobi_lim_cf_mrg_df = pd.merge(cf_brain_3d_df, cf_phen_3d_df, left_index=True, right_index=True, how='left')
    mobi_cf_mrg_df = pd.merge(cf_human_3d_df, mobi_lim_cf_mrg_df, left_index=True, right_index=True, how='left')

    phens_lst = ['Brain', 'ASD', 'EE', 'ID', 'DD', 'SCZ', 'Mix', 'NDDs', 'Control']
    # for feature in features_lst:
    #     box_plotter(data=mobi_cf_mrg_df.loc[(slice(None), feature), phens_lst],
    #                 save_route=(cfg.plots['box-cf'] + '/' + feature + '-cf-all-log1' + '.png'))
    # for feature in features_lst:
    #     violin_plotter(data=mobi_cf_mrg_df.loc[(slice(None), feature), phens_lst],
    #                    save_route=(cfg.plots['vio-cf'] + '/' + feature + '-cf-all-log1' + '.png'))

    ## Content count
    cc_phens_3d_df = mobi_feature_df.groupby(['acc', 'feature', 'Phenotype']).content_count.mean().unstack().sort_index()
    cc_brain_3d_df = mobi_br_subdf.groupby(['acc', 'feature', 'brain']).content_count.mean().unstack().sort_index()
    cc_human_3d_df = mobi_human_subdf.groupby(['acc', 'feature', 'human']).content_count.mean().unstack().sort_index()
    mobi_lim_cc_mrg_df = pd.merge(cc_brain_3d_df, cc_phens_3d_df, left_index=True, right_index=True, how='left')
    mobi_cc_mrg_df = pd.merge(cc_human_3d_df, mobi_lim_cc_mrg_df, left_index=True, right_index=True, how='left')

    # for feature in features_lst:
    #     box_plotter(data=mobi_cc_mrg_df.loc[(slice(None), feature), phens_lst],
    #                 save_route=(cfg.plots['box-cc'] + '/' + feature + '-cc-all-log1' + '.png'))
    for feature in features_lst:
        violin_plotter(data=mobi_cc_mrg_df.loc[(slice(None), feature), phens_lst],
                       save_route=(cfg.plots['vio-cc'] + '/' + feature + '-cc-all-log2' + '.png'))

    ## Length
    len_phens_3d_df = mobi_feature_df.groupby(
        ['acc', 'feature', 'Phenotype']).length.mean().unstack().sort_index()
    len_brain_3d_df = mobi_br_subdf.groupby(['acc', 'feature', 'brain']).length.mean().unstack().sort_index()
    len_human_3d_df = mobi_human_subdf.groupby(['acc', 'feature', 'human']).length.mean().unstack().sort_index()
    mobi_len_lim_mrg_df = pd.merge(len_brain_3d_df, len_phens_3d_df, left_index=True, right_index=True, how='left')
    mobi_len_mrg_df = pd.merge(len_human_3d_df, mobi_len_lim_mrg_df, left_index=True, right_index=True, how='left')
    # for feature in features_lst:
    #     box_plotter(data=mobi_len_mrg_df.loc[(slice(None), feature), phens_lst],
    #                 save_route=(cfg.plots['box-len'] + '/' + feature + '-len-all-log1' + '.png'))
    for feature in features_lst:
        violin_plotter(data=mobi_len_mrg_df.loc[(slice(None), feature), phens_lst],
                       save_route=(cfg.plots['vio-len'] + '/' + feature + '-len-all-log2' + '.png'))




