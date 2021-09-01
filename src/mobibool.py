# this file is supposed to get mobidb+mutposition, is_in_startend and ndd bool, index*, and give a df in the outcome
# with added columns: brain bool, ASD bool, EE bool, ID bool, SCZ bool
# *index (mutual idx of mut_acc_merged_df and mutinfo subdf which is useful to merge mobidb & g4dn_mutinfo_acc,
# could be useless now :/ not sure
## this description will be changed! (August 4th)
import sys

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
    df3 = df3.drop_duplicates().rename(columns={'Phenotype': 'phenotype'})
    all_phens = subdf1.append([df2, df3]).drop_duplicates()
    new_mobi_all_phens = df1_mobi.merge(all_phens, how='left', on='acc')
    return new_mobi_all_phens


def multidx_df_maker(input_dfs_lst, idx_lst):
    # multi-level index Phenotypes: from: https://www.youtube.com/watch?v=tcRGa2soc-c
    # supposed to do this:
    # mobi_disorder_df = mobi_feature_df.groupby(
    # ['acc', 'feature', 'phenotype']).content_fraction.mean().unstack().sort_index()
    cf_multidx_df = input_dfs_lst[0].groupby(idx_lst).content_fraction.mean().unstack().sort_index()
    cf_multidx_df = cf_multidx_df * 100  # converting cf to percentage
    cc_multidx_df = input_dfs_lst[0].groupby(idx_lst).content_count.mean().unstack().sort_index()
    len_multidx_df = input_dfs_lst[1].groupby(idx_lst).length.mean().unstack().sort_index()
    return cf_multidx_df, cc_multidx_df, len_multidx_df


def box_plotter(data, title, ylabel, save_route):
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
    plt.figure(figsize=(12, 6))
    sns.set_theme(style="whitegrid")
    g = sns.violinplot(data=data, split=True, bw=.1).set(title=title, xlabel='Phenotypes', ylabel=ylabel)
    sns.set_style("ticks")
    # g.set_xticklabels(labels=df['Phenotype'].unique().tolist(),rotation=45, va="center", position=(0, -0.02))
    # plt.yscale('log')
    # g.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
    # g.yaxis.set_ticks([np.log10(x) for p in range(-4, 5) for x in np.linspace(10 ** p, 10 ** (p + 1), 10)],minor=True)
    plt.tight_layout()
    plt.savefig(save_route)
    plt.close('all')
    return


def several_plotter(plot_type, inputdf):  # TODO better (shorter) way to write this
    # disorder content
    if plot_type == 'box-cf':
        inputdf = inputdf[inputdf < (0.9 * inputdf.max())]
        for key in feature_dict:
            box_plotter(data=inputdf.loc[(slice(None), key), phens_lst],
                        save_route=(cfg.plots['box-cf'] + '/' + key + '-cf90-x' + '.png'),
                        title=feature_dict[key][0], ylabel='Content (%)')
    # content count
    elif plot_type == 'box-cc':
        for key in feature_dict:
            box_plotter(data=mobi_cont_count_df.loc[(slice(None), key), phens_lst],
                        save_route=(cfg.plots['box-cc'] + '/' + key + '-cc1000-x' + '.png'),
                        title=feature_dict[key][0], ylabel='Content (residues)')
    elif plot_type == 'box-len':
        # length
        inputdf = inputdf[inputdf < 6000]
        box_plotter(data=inputdf.loc[(slice(None)), phens_lst],
                    save_route=(cfg.plots['box-len'] + '/length<6000-x' + '.png'),
                    title='Protein sequence length', ylabel='Residues')
    elif plot_type == 'viol-cf':
        # disorder content
        for key in feature_dict:
            violin_plotter(data=inputdf.loc[(slice(None), key), phens_lst],
                           save_route=(cfg.plots['vio-cf'] + '/' + key + '-cf-90-x' + '.png'),
                           title=feature_dict[key][0], ylabel='Content (%)')
    elif plot_type == 'viol-cc':
        # content count
        for key in feature_dict:
            violin_plotter(data=mobi_cont_count_df.loc[(slice(None), key), phens_lst],
                           save_route=(cfg.plots['vio-cc'] + '/' + key + '-cc-1000-x' + '.png'),
                           title=feature_dict[key][0], ylabel='Content (residues)')
    elif plot_type == 'viol-len':
        ## Length
        inputdf = inputdf[inputdf < 6000]
        violin_plotter(data=inputdf.loc[(slice(None)), phens_lst],
                       save_route=(cfg.plots['vio-len'] + '/length-below6000-x' + '.png'),
                       title='Protein sequence length', ylabel='Residues')


# def cc_feature_filterer(cc_org_df, lim, selected_features):  # what if I take this out
#     cc_new_df = cc_org_df[cc_org_df <= lim]
#     cc_new_df = cc_new_df.loc[(slice(None), selected_features), phens_lst]
#     return cc_new_df


# # def plot_distinguisher(input_plot_func, input_filterer_func):
# #     def inner_plotter(data, title, ylabel, save_route):
# #         def filtere(data, lim, feats):
# #
# #     return inner_plotter()


if __name__ == '__main__':
    ## selected features  (10)
    features_lst = ['prediction-disorder-mobidb_lite', 'prediction-low_complexity-merge',
                    'prediction-lip-anchor', 'homology-domain-merge',
                    'prediction-signal_peptide-uniprot', 'prediction-transmembrane-uniprot',
                    'curated-conformational_diversity-merge', 'derived-binding_mode_disorder_to_disorder-mobi',
                    'derived-binding_mode_disorder_to_order-mobi', 'derived-mobile_context_dependent-th_90',
                    'prediction-disorder-th_50']
    titles_lst = ['Disorder content - Mobidb-lite', 'Low complexity', 'Linear Interacting Peptides - Anchor',
                  'Protein domains', 'Signal peptides - Uniprot', 'Transmembrane helices - Uniprot',
                  'Conformational diversity', 'Binding mode - disorder to disorder',
                  'Binding mode - disorder to order', 'Protein mobility', 'Disorder - majority']
    titles_lst = [[i] for i in titles_lst]
    # this is a dict with feature names as key and their plot titles as value (then dfs for cc will be added as values)
    feature_dict = dict(zip(features_lst, titles_lst))
    cc_lim_lst = [1000, 800, 1000, 1200, 70, 600, 700, 250, 100, 500, 1200]
    cc_lim_feature_dict = dict(zip(features_lst, cc_lim_lst))
    cf_lim_lst = [100, 60, 100, 100, 100, 100, 100, 100, 100, 100, 100]
    # now add content_count limit of each feature to your dict as second value (idx=1)
    # and content fraction limit as idx=2
    for feature in features_lst:
        for cf_lim in cf_lim_lst:
            feature_dict[feature].append(cc_lim_feature_dict[feature])
            feature_dict[feature].append(cf_lim)

    # in the end use decorator thing with the @
    phens_lst = ['Human', 'Brain', 'ASD', 'EE', 'ID', 'DD', 'SCZ', 'NDDs', 'Mix', 'Control']
    ## import dfs # (mobidb)
    mobidb = pd.read_csv(cfg.data[''] + '/mobidb_result.tsv', sep='\t')
    # NDD , could also specify index_col= ..., and pass a list for multiple idxs
    ndd_subdf = pd.read_csv(cfg.data['phens-fdr'] + '/acc-phen-5percentFDR.csv')
    ndd_subdf = ndd_subdf.drop_duplicates()  # (4531, 3)
    ndd_pr_lst = ndd_subdf['acc'].unique().tolist()  # 1308 proteins  => these are all, not the selected phens
    ## brain
    brain_prot_lst = bd.brain_pr_lst_generator()  # n: 8297
    brain_subdf = DataFrame(brain_prot_lst, columns=['acc'])
    # mutual_brain_ndd_prs_lst = [i for i in brain_prot_lst if i in ndd_pr_lst]  # 455

    # new mobidb with one column for phens of NDD, human, and brain
    mobidb = mobi_phens_col_maker(mobidb, brain_subdf, ndd_subdf)
    mobi_feature_df = mobidb[mobidb.feature.isin(features_lst)]
    # multi-idx-dfs
    mobi_disorder_df, mobi_cont_count_df, mobi_length_df = multidx_df_maker(
        [mobi_feature_df, mobidb], ['acc', 'feature', 'phenotype'])

    # filtering data
    # def cc_feature_filterer(cc_org_df, lim):
    #     cc_new_df = cc_org_df[cc_org_df <= lim]
    #     return cc_new_df

    # # box plots for disorder_content, content_count and length
    # several_plotter('box-cf', mobi_disorder_df)
    # several_plotter('box-cc', mobi_cont_count_df)  # this still needs to be filtered
    # several_plotter('box-len', mobi_length_df)
    # # violin plots for disorder_content, content_count and length
    # several_plotter('viol-cf', mobi_disorder_df)
    # several_plotter('viol-cc', mobi_cont_count_df)  # needs filteration
    # several_plotter('viol-len', mobi_length_df)
    #
    # ## writing data statistics to CSV
    # pd.set_option('display.max_columns', None)
    # pd.set_option('display.max_rows', None)
    # for each_f in features_lst:
    #     mobi_disorder_df.loc[(slice(None), each_f), phens_lst].describe().T. \
    #         to_csv(cfg.data['phens'] + '/' + each_f + '-cf.csv')
    #     mobi_cont_count_df.loc[(slice(None), each_f), phens_lst].describe().T. \
    #         to_csv(cfg.data['phens'] + '/' + each_f + '-cc.csv')
    # mobi_length_df.loc[slice(None), phens_lst].describe().T.to_csv(cfg.data['phens'] + '/length-stats.csv')
