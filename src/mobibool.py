# this file is supposed to get mobidb+mutposition, is_in_startend and ndd bool, index*, and give a df in the outcome
# with added columns: brain bool, ASD bool, EE bool, ID bool, SCZ bool
# *index (mutual idx of mut_acc_merged_df and mutinfo subdf which is useful to merge mobidb & g4dn_mutinfo_acc,
# could be useless now :/ not sure
## this description will be changed! (August 4th)

import pandas as pd
import numpy as np
from pandas import DataFrame
import config as cfg
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


def multidx_df_maker(input_dfs_lst, idx_lst):
    # multi-level index Phenotypes: from: https://www.youtube.com/watch?v=tcRGa2soc-c
    # supposed to do this:
    # mobi_disorder_df = mobi_feature_df.groupby(
    # ['acc', 'feature', 'phenotype']).content_fraction.mean().unstack().sort_index()
    cf_multidx_df = input_dfs_lst[0].groupby(idx_lst).content_fraction.mean().unstack().sort_index()
    cc_multidx_df = input_dfs_lst[0].groupby(idx_lst).content_count.mean().unstack().sort_index()
    len_multidx_df = input_dfs_lst[1].groupby(idx_lst).length.mean().unstack().sort_index()
    return cf_multidx_df, cc_multidx_df, len_multidx_df


def box_plotter(data, title, ylabel, txt_name, save_route):
    plt.figure(figsize=(60, 60))  # bigger figsize to have xticklabels shown
    g = sns.catplot(data=data, kind="box").set(title=title, xlabel='Phenotypes', ylabel=ylabel)
    sns.set_style("ticks")
    # dict[dict_key] = data.describe().T
    g.set_xticklabels(rotation=45, va="center", position=(0, -0.02))
    # plt.yscale('log')
    plt.tight_layout()
    plt.savefig(save_route)
    plt.close('all')
    return


def violin_plotter(data, title, save_route, ylabel):
    plt.figure(figsize=(12, 12))
    sns.set_theme(style="whitegrid")
    g = sns.violinplot(data=data, split=True, bw=.1).set(title=title, xlabel='Phenotypes', ylabel=ylabel)
    sns.set_style("ticks")
    # g.set_xticklabels(labels=df['Phenotype'].unique().tolist(),rotation=45, va="center", position=(0, -0.02))
    # plt.yscale('log')
    plt.tight_layout()
    plt.savefig(save_route)
    plt.close('all')
    info_df = data.describe().T
    return info_df


if __name__ == '__main__':
    ## selected features  (10)
    features_lst = ['prediction-disorder-mobidb_lite', 'prediction-low_complexity-merge',
                    'prediction-lip-anchor', 'homology-domain-merge',
                    'prediction-signal_peptide-uniprot', 'prediction-transmembrane-uniprot',
                    'curated-conformational_diversity-merge', 'derived-binding_mode_disorder_to_disorder-mobi',
                    'derived-binding_mode_disorder_to_order-mobi', 'derived-mobile_context_dependent-th_90']
    titles_lst = ['Disorder content - Mobidb-lite', 'Low complexity', 'Linear Interacting Peptides - Anchor',
                  'Protein domains', 'Signal peptides - Uniprot', 'Transmembrane helices - Uniprot',
                  'Conformational diversity', 'Binding mode - disorder to disorder',
                  'Binding mode - disorder to order', 'Protein mobility']

    ## selected phenotypes
    phens_lst = ['Human', 'Brain', 'ASD', 'EE', 'ID', 'DD', 'SCZ', 'Mix', 'NDDs', 'Control']
    ## import dfs # (mobidb)
    mobidb = pd.read_csv(cfg.data[''] + '/mobidb_result.tsv', sep='\t')
    # brain
    brain_prot_lst = bd.brain_pr_lst_generator()  # n: 8428
    brain_subdf = DataFrame(brain_prot_lst, columns=['acc'])
    # NDD , could also specify index_col= ..., and pass a list for multiple idxs
    ndd_subdf = pd.read_csv(cfg.data['gene4'] + '/positive_cand_g4mobi_concat.csv', usecols=["acc_x", "Phenotype"])

    # new mobidb with one column for phens of NDD, human, and brain
    mobidb = mobi_phens_col_maker(mobidb, brain_subdf, ndd_subdf)
    mobi_feature_df = mobidb[mobidb.feature.isin(features_lst)]
    # multi-idx-dfs
    mobi_disorder_df, mobi_cont_count_df, mobi_length_df = multidx_df_maker(
        [mobi_feature_df, mobidb], ['acc', 'feature', 'phenotype'])
    # filtering data
    # TODO maybe different filteration per each feature, cause the outcome is different for each, and some still have
    # outliers and are rather difficult to analyse
    mobi_disorder_df = mobi_disorder_df[mobi_disorder_df < (0.9 * mobi_disorder_df.max())] * 100
    # mobi_disorder_df = mobi_disorder_df * 100
    mobi_cont_count_df = mobi_cont_count_df[mobi_cont_count_df <= 1000]
    mobi_length_df = mobi_length_df[mobi_length_df < 6000]

    ## plot (boxplot)
    # disorder content
    for (feature, title) in zip(features_lst, titles_lst):
        box_plotter(data=mobi_disorder_df.loc[(slice(None), feature), phens_lst],
                    save_route=(cfg.plots['box-cf'] + '/' + feature + '-cf90' + '.png'),
                    title=title, ylabel='Content (%)', txt_name=feature)
    # content count
    for (feature, title) in zip(features_lst, titles_lst):
        box_plotter(data=mobi_cont_count_df.loc[(slice(None), feature), phens_lst],
                    save_route=(cfg.plots['box-cc'] + '/' + feature + '-cc1000' + '.png'),
                    title=title, ylabel='Content (residues)', txt_name=feature)
    # length
    box_plotter(data=mobi_length_df.loc[(slice(None)), phens_lst],
                save_route=(cfg.plots['box-len'] + '/length<6000' + '.png'),
                title='Protein sequence length', ylabel='Residues', txt_name='length')
    ## plot (violinplot)
    # disorder content
    for (feature, title) in zip(features_lst, titles_lst):
        violin_plotter(data=mobi_disorder_df.loc[(slice(None), feature), phens_lst],
                       save_route=(cfg.plots['vio-cf'] + '/' + feature + '-cf-90' + '.png'),
                       title=title, ylabel='Content (%)')
    # content count
    for (feature, title) in zip(features_lst, titles_lst):
        violin_plotter(data=mobi_cont_count_df.loc[(slice(None), feature), phens_lst],
                       save_route=(cfg.plots['vio-cc'] + '/' + feature + '-cc-1000' + '.png'),
                       title=title, ylabel='Content (residues)')
    ## Length
    violin_plotter(data=mobi_length_df.loc[(slice(None)), phens_lst],
                   save_route=(cfg.plots['vio-len'] + '/length-below6000' + '.png'),
                   title='Protein sequence length', ylabel='Residues')
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    for feature in features_lst:
        mobi_disorder_df.loc[(slice(None), feature), phens_lst].describe().T\
            .to_csv(cfg.data['phens'] + '/' + feature + '.csv')
        # ds.to_csv(os.path.join('..', 'data', 'test.csv'))
        # .describe().T.to_csv(cfg.data['phens']+'/test.csv')
