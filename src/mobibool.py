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
    brain_pr_df['brain'] = True
    df = df.merge(brain_pr_df, how='left', on='acc')
    # df = df.set_index("acc")

    # ## NDD proteins
    df_g4 = pd.read_csv(cfg.data['gene4'] + '/positive_cand_g4mobi_concat.csv', usecols=["acc_x", "Phenotype"])
    df_g4 = df_g4.drop_duplicates()
    # df_g4 = df_g4.set_index("acc_x")


    #TODO think about a flatter (maybe more redundant) dataframe. Use the Phenotype as index in addition to acc
    df = df.merge(df_g4, how='left', left_on='acc', right_on='acc_x')

    # #
    # df_g4['val'] = True
    # df_g4 = df_g4.set_index(["acc_x", "Phenotype"]).unstack()
    # df_g4 = df_g4.loc[:, 'val']
    # df_g4['ndd'] = True
    # df = df.merge(df_g4, how='left', left_on='acc', right_on='acc_x')
    #
    # # disorder content
    # columns = ['acc', 'content_fraction', 'ndd', 'brain', 'ASD']
    # df_dc = df[df['feature'] == 'prediction-disorder-mobidb_lite'][columns]
    # df_dc = df_dc.set_index(['acc', 'content_fraction'])
    # df_dc[df_dc.notna()] = df_dc.index.get_level_values('content_fraction')





