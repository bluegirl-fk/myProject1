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


def mobi_phens_col_maker(df1_mobi, df2, df3):
    # mobidb
    df1_mobi = df1_mobi.drop(columns='start..end')
    subdf1 = df1_mobi[['acc']]
    subdf1['phenotype'] = 'Human'
    # brain
    df2['phenotype'] = 'Brain'
    # ndd
    df3 = df3.drop_duplicates().rename(columns={'acc_x': 'acc', 'Phenotype': 'phenotype'})
    all_phens = subdf1.append([df2, df3]).drop_duplicates()
    new_mobi_all_phens = df1_mobi.merge(all_phens, how='left', on='acc')
    return new_mobi_all_phens


def box_plotter(data, save_route):
    plt.figure(figsize=(60, 60))  # bigger figsize to have xticklabels shown
    g = sns.catplot(data=data, kind="box")
    sns.set_style("ticks")
    g.set_xticklabels(rotation=45, va="center", position=(0, -0.02))
    # plt.yscale('log')
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
    # g.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
    # g.yaxis.set_ticks([np.log10(x) for p in range(-4, 5) for x in np.linspace(10 ** p, 10 ** (p + 1), 10)],minor=True)
    plt.tight_layout()
    plt.savefig(save_route)
    plt.close()


if __name__ == '__main__':
    ## selected features
    features_lst = ['prediction-disorder-mobidb_lite', 'prediction-low_complexity-merge',
                    'prediction-lip-anchor', 'homology-domain-merge']

    mobidb = pd.read_csv(cfg.data[''] + '/mobidb_result.tsv', sep='\t')
    ## brain proteins
    brain_prot_lst = bd.brain_pr_lst_generator()  # n: 8428
    brain_subdf = DataFrame(brain_prot_lst, columns=['acc'])
    ## NDD proteins, could also specify index_col= ..., and pass a list for multiple idxs
    ndd_subdf = pd.read_csv(cfg.data['gene4'] + '/positive_cand_g4mobi_concat.csv', usecols=["acc_x", "Phenotype"])
    # new mobidb with one column for phens of NDD, human, and brain
    mobidb = mobi_phens_col_maker(mobidb, brain_subdf, ndd_subdf)
    mobi_feature_df = mobidb[mobidb.feature.isin(features_lst)]
    # multi-level index Phenotypes: from: https://www.youtube.com/watch?v=tcRGa2soc-c
    # disorder content
    mobi_disorder_df = mobi_feature_df.groupby(['acc', 'feature', 'phenotype']).content_fraction.mean().unstack().sort_index()
    # content count
    mobi_cont_count_df = mobi_feature_df.groupby(['acc', 'feature', 'Phenotype']).content_count.mean().unstack().sort_index()
    # Length, original mobidb, not the one with selected features
    mobi_length_df = mobidb.groupby(['acc', 'feature', 'Phenotype']).length.mean().unstack().sort_index()
    mobi_length_df = mobi_length_df[mobi_length_df < 6000]
    ## Selected phenotypes for plots
    phens_lst = ['Brain', 'ASD', 'EE', 'ID', 'DD', 'SCZ', 'Mix', 'NDDs', 'Control']

    # # for feature in features_lst:
    # #     box_plotter(data=mobi_cf_mrg_df.loc[(slice(None), feature), phens_lst],
    # #                 save_route=(cfg.plots['box-cf'] + '/' + feature + '-cf-all-log1' + '.png'))
    # # for feature in features_lst:
    # #     violin_plotter(data=mobi_cf_mrg_df.loc[(slice(None), feature), phens_lst],
    # #                    save_route=(cfg.plots['vio-cf'] + '/' + feature + '-cf-all-log1' + '.png'))
    #

    # # for feature in features_lst:
    # #     box_plotter(data=mobi_cc_mrg_df.loc[(slice(None), feature), phens_lst],
    # #                 save_route=(cfg.plots['box-cc'] + '/' + feature + '-cc-all-log1' + '.png'))
    # for feature in features_lst:
    #     violin_plotter(data=mobi_cc_mrg_df.loc[(slice(None), feature), phens_lst],
    #                    save_route=(cfg.plots['vio-cc'] + '/' + feature + '-cc-all-log2' + '.png'))
    #
    # ## Length
    # for feature in features_lst:
    #     box_plotter(data=mobi_len_mrg_df.loc[(slice(None), feature), phens_lst],
    #                 save_route=(cfg.plots['box-len'] + '/' + feature + '-len6000' + '.png'))
    # for feature in features_lst:
    #     violin_plotter(data=mobi_len_mrg_df.loc[(slice(None), feature), phens_lst],
    #                    save_route=(cfg.plots['vio-len'] + '/' + feature + '-len6000' + '.png'))
