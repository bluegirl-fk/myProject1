# this file is supposed to get mobidb+mutposition, is_in_startend and ndd bool, index*, and give a df in the outcome
# with added columns: brain bool, ASD bool, EE bool, ID bool, SCZ bool
# *index (mutual idx of mut_acc_merged_df and mutinfo subdf which is useful to merge mobidb & g4dn_mutinfo_acc,
# could be useless now :/ not sure

import pandas as pd
import numpy as np
from pandas import DataFrame
import config as cfg
import sys
import brain as bd
import matplotlib.pyplot as plt
import seaborn as sns

# def phen_col_merger(phenotype, mobi_df):
#     phen_lst = positive_cand_g4mobi_df.loc[positive_cand_g4mobi_df['Phenotype'] == phenotype, 'acc_x'].unique().tolist()
#     phen_cand_df = DataFrame(phen_lst, columns=['acc'])
#     phen_cand_df[phenotype] = 1
#     mobi_bool_df = mobi_df.merge(phen_cand_df, how='left', on='acc')
#     mobi_bool_df[phenotype] = mobi_bool_df[phenotype].fillna(0)
#     return mobi_bool_df, phen_lst


if __name__ == '__main__':
    df = pd.read_csv(cfg.data[''] + '/mobidb_result.tsv', sep='\t')
    ## brain proteins
    brain_prot_lst = bd.brain_pr_lst_generator()  # n: 8428
    brain_pr_df = DataFrame(brain_prot_lst, columns=['acc'])
    brain_pr_df['brain'] = 1
    df = df.merge(brain_pr_df, how='left', on='acc')
    # df = df.set_index("acc")

    # ## NDD proteins
    df_g4 = pd.read_csv(cfg.data['gene4'] + '/positive_cand_g4mobi_concat.csv', usecols=["acc_x", "Phenotype"])
    df_g4 = df_g4.drop_duplicates()
    # df_g4 = df_g4.set_index("acc_x")

    # # TODO think about a flatter (maybe more redundant) dataframe. Use the Phenotype as index in addition to acc
    df = df.merge(df_g4, how='left', left_on='acc', right_on='acc_x')

    # trying multi-level index here: from: https://www.youtube.com/watch?v=tcRGa2soc-c
    # TODO: (in the end you should add the brain to the Phenotype column) merge it to the phens existing column
    df = df.drop(columns=['start..end', 'acc_x'])
    features_lst = ['prediction-disorder-mobidb_lite', 'prediction-low_complexity-merge',
                    'prediction-lip-anchor', 'homology-domain-merge']
    midx = df[df.feature.isin(features_lst)]  # (194383, 7)

    ser = midx.groupby(['acc', 'feature', 'Phenotype']).content_fraction.mean()
    ser_df = ser.unstack()
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





    # pivot_df = df.pivot_table(values='content_fraction', index='acc', columns='Phenotype') # not the same cuz not
    # taking into account the features


    # ## prev file
    # prev_df = pd.read_csv(cfg.data['csv'] + '/mobi_phens_bool.csv')
    # prev_df.brain[prev_df.brain == True] = prev_df.content_fraction
    # prev_df.ndd[prev_df.ndd == 1] = prev_df.content_fraction
    # prev_df.ASD[prev_df.ASD == 1] = prev_df.content_fraction
    # prev_df.EE[prev_df.EE == 1] = prev_df.content_fraction
    # prev_df.ID[prev_df.ID == 1] = prev_df.content_fraction
    # prev_df.SCZ[prev_df.SCZ == 1] = prev_df.content_fraction
    # cols = ["ndd", "ASD", "EE", "ID", "SCZ"]
    # prev_df[cols] = prev_df[cols].replace({0: np.nan, 0: np.nan})
    #

    #
    # mobi_lite_df = prev_df.loc[prev_df['feature'] == 'prediction-disorder-mobidb_lite']
    # mobi_lite_df = mobi_lite_df[['acc', 'brain', 'ndd', 'ASD', 'EE', 'ID', 'SCZ']]
    # # mobi_lite_df = mobi_lite_df.rename(columns={'content_fraction': 'human'})
    #
    # plt.figure(figsize=(20, 20))  # bigger figsize to have xticklabels shown
    # plot = sns.catplot(data=mobi_lite_df, kind="box")
    # # plt.tight_layout()
    # plt.show()
    #
    # # plt.figure(figsize=(20, 20))
    # # sns.set_theme(style="whitegrid")
    # # # Load the example tips dataset
    # # # Draw a nested violinplot and split the violins for easier comparison
    # # sns.violinplot(data=mobi_lite_df)
    # # plt.show()
    # # #
    # # df_g4['val'] = True
    # # df_g4 = df_g4.set_index(["acc_x", "Phenotype"]).unstack()
    # # df_g4 = df_g4.loc[:, 'val']
    # # df_g4['ndd'] = True
    # # df = df.merge(df_g4, how='left', left_on='acc', right_on='acc_x')
    # #
    # # # disorder content
    # # columns = ['acc', 'content_fraction', 'ndd', 'brain', 'ASD']
    # # df_dc = df[df['feature'] == 'prediction-disorder-mobidb_lite'][columns]
    # # df_dc = df_dc.set_index(['acc', 'content_fraction'])
    # # df_dc[df_dc.notna()] = df_dc.index.get_level_values('content_fraction')
