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
import varAnalysis as var


def mobi_phens_col_maker(df1_mobi, df2, df3):
    # mobidb
    # df1_mobi = df1_mobi.drop(columns='start..end')
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
    len_multidx_df = input_dfs_lst[1].groupby(['acc', 'phenotype']).length.mean().unstack().sort_index()  # this is to
    # avoid duplications cuz length is the same and does not differ with mobidb features,
    # so the features should not be taken into account
    return cf_multidx_df, cc_multidx_df, len_multidx_df


def box_plotter(data, title, ylabel, ylim, save_route):
    plt.figure(figsize=(60, 60))  # bigger figsize to have xticklabels shown
    g = sns.catplot(data=data, kind="box").set(title=title, xlabel='Phenotypes', ylabel=ylabel, ylim=(-5, ylim))
    sns.set_style("ticks")
    g.set_xticklabels(rotation=45, va="center", position=(0, -0.02))
    plt.tight_layout()
    plt.savefig(save_route)
    plt.close('all')
    return


def violin_plotter(data, title, save_route, ylabel, ylim):
    plt.figure(figsize=(12, 6))
    sns.set_theme(style="whitegrid")
    g = sns.violinplot(data=data, split=True, bw=.1).set(title=title, xlabel='Phenotypes', ylabel=ylabel)
    sns.set_style("ticks")
    # g.set_xticklabels(rotation=45, va="center", position=(0, -0.02))
    plt.ylim(-5, ylim)
    plt.tight_layout()
    plt.savefig(save_route)
    plt.close('all')
    return


def draw_barplot(x, y, data, xticklabel, yscale, save_rout):  # input is DF, not list
    plt.figure(figsize=(12, 6))  # bigger figsize to have xticklabels shown
    sns.set_style("ticks")
    g = sns.barplot(x=x, y=y, data=data)
    # sns.despine(trim=True, offset=2)
    g.set_xticklabels(xticklabel, rotation=0, va="center", position=(0, -0.02))
    sns.color_palette("pastel")
    plt.yscale(yscale)
    plt.tight_layout()
    plt.savefig(save_rout)
    plt.show()
    return


def sig_pep_percent_df_maker():
    # this gets a df with the selected peptides count and then divides them by all proteins of that phenotype or
    # brn, human, ... and creates percentage, all Prs count is taken from length_df cuz it does not have redundancy
    sig_peptide_subdf = pd.read_csv(
        cfg.data['all-var-desc-cc'] + '/prediction-signal_peptide-uniprot-cc.csv',
        usecols=['phenotype', 'count'])
    all_phen_pr_count_df = pd.read_csv(cfg.data['all-var-desc-len'] + '/length-stats.csv', usecols=['phenotype', 'count'])
    sig_pep_mrged_df = pd.merge(sig_peptide_subdf, all_phen_pr_count_df, on='phenotype')
    sig_pep_mrged_df = sig_pep_mrged_df.rename(
        columns={'phenotype': 'Phenotypes', 'count_x': 'Signal peptide count', 'count_y': 'count_all'})
    sig_pep_mrged_df['Signal peptide percentage'] = (sig_pep_mrged_df['Signal peptide count'] * 100) / \
                                                    sig_pep_mrged_df['count_all']
    return sig_pep_mrged_df


def transmem_pr_percent_df_maker():
    all_phen_pr_count_df = pd.read_csv(cfg.data['all-var-desc-len'] + '/length-stats.csv', usecols=['phenotype', 'count'])
    transmemb_subdf = pd.read_csv(
        cfg.data['all-var-desc-cc'] + '/prediction-transmembrane-uniprot-cc.csv',
        usecols=['phenotype', 'count'])
    transmem_mrg_df = pd.merge(transmemb_subdf, all_phen_pr_count_df, on='phenotype')
    transmem_mrg_df = transmem_mrg_df.rename(
        columns={'phenotype': 'Phenotypes', 'count_x': 'Transmembrane protein count', 'count_y': 'count_all'})
    transmem_mrg_df['Transmembrane protein percentage'] = (transmem_mrg_df['Transmembrane protein count'] * 100) / \
                                                          transmem_mrg_df['count_all']
    return transmem_mrg_df


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
    cc_lim_lst = [1010, 810, 1210, 1810, 75, 610, 900, 360, 105, 510, 1810]
    cc_lim_feature_dict = dict(zip(features_lst, cc_lim_lst))
    cf_lim_lst = [105, 65, 105, 105, 105, 105, 105, 105, 105, 105, 105]
    # cf_lim_lst = [110, 65, 110, 110, 110, 110, 110, 110, 110, 110, 110] for violins I used 110 instead of 105
    cf_lim_feature_dict = dict(zip(features_lst, cf_lim_lst))
    # now add content_count limit of each feature to your dict as second value (idx=1)
    # and content fraction limit as idx=2
    for feature in features_lst:
        feature_dict[feature].append(cc_lim_feature_dict[feature])
        feature_dict[feature].append(cf_lim_feature_dict[feature])

    # in the end use decorator thing with the @
    phens_lst = ['Human', 'Brain', 'ASD', 'EE', 'ID', 'DD', 'SCZ', 'NDDs', 'Control']
    ## import dfs # (mobidb)
    # mobidb = pd.read_csv(cfg.data['phens'] + '/mobidb-results+ndd-tsv-damiano-shared.tsv')
    # # NDD , could also specify index_col= ..., and pass a list for multiple idxs
    # ndd_subdf = pd.read_csv(cfg.data['phens-fdr'] + '/acc-phen-5percentFDR.csv')
    # ndd_subdf = ndd_subdf.drop_duplicates()  # (4531, 3)
    # ndd_pr_lst = ndd_subdf['acc'].unique().tolist()  # 1308 proteins  => these are all, not the selected phens
    # ## brain
    # brain_prot_lst = bd.brain_pr_lst_generator()  # n: 8320
    # brain_subdf = DataFrame(brain_prot_lst, columns=['acc'])
    # # mutual_brain_ndd_prs_lst = [i for i in brain_prot_lst if i in ndd_pr_lst]  # 455

    ## Mobidb, Brain and Ndd dfs with all variations (in/out of IDR)
    mobidb, brain_subdf, ndd_subdf = var.all_vars_or_vars_inidr('all')

    # new mobidb with one column for phens of NDD, human, and brain
    mobidb = mobi_phens_col_maker(mobidb, brain_subdf, ndd_subdf)
    mobidb.to_csv(cfg.data['vars'] + '/all-vars-mobidb-plus-phenotype-column.csv')

    mobi_feature_df = mobidb[mobidb.feature.isin(features_lst)]
    # multi-idx-dfs
    mobi_disorder_df, mobi_cont_count_df, mobi_length_df = multidx_df_maker(
        [mobi_feature_df, mobidb], ['acc', 'feature', 'phenotype'])
    ## Boxplots
    # content count
    for key in feature_dict.keys():
        box_plotter(data=mobi_cont_count_df.loc[(slice(None), key), phens_lst],
                    save_route=(cfg.plots['avb-cc'] + '/' + key + '-cc' + '.png'),
                    title=feature_dict[key][0], ylabel='Residues count', ylim=feature_dict[key][1])
    # content fraction
    for key in feature_dict.keys():
        box_plotter(data=mobi_disorder_df.loc[(slice(None), key), phens_lst],
                    save_route=(cfg.plots['avb-cf'] + '/' + key + '-cf' + '.png'),
                    title=feature_dict[key][0], ylabel='Content (%)', ylim=feature_dict[key][2])
    # Length
    box_plotter(data=mobi_length_df.loc[(slice(None)), phens_lst],
                save_route=(cfg.plots['avb-len'] + '/' + 'len4200' + '.png'),
                title='Protein sequence length', ylabel='Residues count', ylim=4200)
    ## Violin plots
    # content count
    for key in feature_dict.keys():
        violin_plotter(data=mobi_cont_count_df.loc[(slice(None), key), phens_lst],
                       save_route=(cfg.plots['avv-cc'] + '/' + key + '-cc' + '.png'),
                       title=feature_dict[key][0], ylabel='Residues count', ylim=feature_dict[key][1])
    # content fraction
    for key in feature_dict.keys():
        violin_plotter(data=mobi_disorder_df.loc[(slice(None), key), phens_lst],
                       save_route=(cfg.plots['avv-cf'] + '/' + key + '-cf' + '.png'),
                       title=feature_dict[key][0], ylabel='Content (%)', ylim=feature_dict[key][2])
    # Length
    violin_plotter(data=mobi_length_df.loc[(slice(None)), phens_lst],
                   save_route=(cfg.plots['avv-len'] + '/' + 'len4200' + '.png'),
                   title='Protein sequence length', ylabel='Residues count', ylim=4200)

    ## writing data statistics to CSV
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    for each_f in features_lst:
        mobi_disorder_df.loc[(slice(None), each_f), phens_lst].describe().T. \
            to_csv(cfg.data['all-var-desc-cf'] + '/' + each_f + '-cf.csv')
        mobi_cont_count_df.loc[(slice(None), each_f), phens_lst].describe().T. \
            to_csv(cfg.data['all-var-desc-cc'] + '/' + each_f + '-cc.csv')
    mobi_length_df.loc[slice(None), phens_lst].describe().T.to_csv(cfg.data['all-var-desc-len'] + '/length-stats.csv')

    sig_pep_percent_df = sig_pep_percent_df_maker()
    transmem_pr_percent_df = transmem_pr_percent_df_maker()

    ## Barplots
    # Signal peptide
    draw_barplot(x='Phenotypes', y='Signal peptide percentage', data=sig_pep_percent_df, xticklabel=phens_lst,
                 yscale='linear', save_rout=cfg.plots['avar-bar'] + '/sig-peptide-percent.png')
    draw_barplot(x='Phenotypes', y='Signal peptide count', data=sig_pep_percent_df, xticklabel=phens_lst,
                 yscale='log', save_rout=cfg.plots['avar-bar'] + '/sig-peptide-count.png')
    # Transmembrane protein
    draw_barplot(x='Phenotypes', y='Transmembrane protein percentage', data=transmem_pr_percent_df, xticklabel=phens_lst,
                 yscale='linear', save_rout=cfg.plots['avar-bar'] + '/transmemb-prots-percent.png')
    draw_barplot(x='Phenotypes', y='Transmembrane protein count', data=transmem_pr_percent_df, xticklabel=phens_lst,
                 yscale='log', save_rout=cfg.plots['avar-bar'] + '/transmemb-prots-count.png')
