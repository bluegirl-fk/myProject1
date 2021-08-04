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
import seaborn as sns


if __name__ == '__main__':
    mobidb = pd.read_csv(cfg.data[''] + '/mobidb_result.tsv', sep='\t')
    ## brain proteins
    brain_prot_lst = bd.brain_pr_lst_generator()  # n: 8428
    brain_pr_df = DataFrame(brain_prot_lst, columns=['acc'])
    brain_pr_df['brain'] = 1
    mobidb = mobidb.merge(brain_pr_df, how='left', on='acc')
    ## NDD proteins, could also specify index_col= ..., and pass a list for multiple idxs
    phen_df = pd.read_csv(cfg.data['gene4'] + '/positive_cand_g4mobi_concat.csv', usecols=["acc_x", "Phenotype"])
    phen_df = phen_df.drop_duplicates()
    ## merge to get Phenotypes column
    mobidb = mobidb.merge(phen_df, how='left', left_on='acc', right_on='acc_x')
    mobidb = mobidb.drop(columns=['start..end', 'acc_x'])
    # trying multi-level index here: from: https://www.youtube.com/watch?v=tcRGa2soc-c
    # TODO: in the end you should add the brain to the Phenotype column
    features_lst = ['prediction-disorder-mobidb_lite', 'prediction-low_complexity-merge',
                    'prediction-lip-anchor', 'homology-domain-merge']
    mobi_feature_df = mobidb[mobidb.feature.isin(features_lst)]  # (194383, 7)

    mobi_3d_series = mobi_feature_df.groupby(['acc', 'feature', 'Phenotype']).content_fraction.mean()
    mobi_3d_df = mobi_3d_series.unstack()
    ## how to access the values
    # ser.loc['A6NHU9']
    # ser.loc['A6NHU9', 'prediction-lip-anchor']
    # ser.loc[:, 'prediction-lip-anchor']
    ## works the same with the dataframe
    # ser_df.loc['A6NHU9']
    # ser_df.loc[('A6NHU9', 'prediction-lip-anchor'), :]
    # ser_df.loc[('A6NHU9', 'prediction-lip-anchor'), 'ASD']
    # ser_df.loc[('A6NHU9', ['prediction-lip-anchor', 'prediction-low_complexity-merge']), :]
    # ser_df.loc[('A6NHU9', ['prediction-lip-anchor', 'prediction-low_complexity-merge']), ['ASD', 'EE', 'ID', 'SCZ']]
    # here slice(None) is used to get all ACCs but for specific items of the other idx (feature), so we don't use : here
    # ser_df.loc[(slice(None), ['prediction-lip-anchor', 'prediction-low_complexity-merge']), :]
    ## for merging, could be useful for brain! try it later
    # mrged_df = pd.merge(brain_subdf_sameidx, phens_subdf_same idxkinda, left_index=True, right_index=True, how='left')

    # pivot_df = df.pivot_table(values='content_fraction', index='acc', columns='Phenotype') # not the same cuz not
    # taking into account the features

    # plt.figure(figsize=(60, 60))  # bigger figsize to have xticklabels shown
    # plot = sns.catplot(data=ser_df.loc[(slice(None), 'prediction-disorder-mobidb_lite'), :], kind="box")
    # sns.set_style("ticks")
    # plot.set_xticklabels(rotation=45, va="center", position=(0, -0.02))
    # plt.tight_layout()
    # plt.show()

    plt.figure(figsize=(12, 12))
    sns.set_theme(style="whitegrid")
    # Load the example tips dataset
    # Draw a nested violinplot and split the violins for easier comparison
    plot = sns.violinplot(data=mobi_3d_df.loc[(slice(None), 'prediction-disorder-mobidb_lite'), :])
    sns.set_style("ticks")
    # plot.set_xticklabels(labels=df['Phenotype'].unique().tolist(),rotation=45, va="center", position=(0, -0.02))
    plt.tight_layout()
    plt.show()


