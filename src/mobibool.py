# this file is supposed to get mobidb+mutposition, is_in_startend and ndd bool, index*, and give a df in the outcome
# with added columns: brain bool, ASD bool, EE bool, ID bool, SCZ bool
# *index (mutual idx of mut_acc_merged_df and mutinfo subdf which is useful to merge mobidb & g4dn_mutinfo_acc,
# could be useless now :/ not sure

import pandas as pd
from pandas import DataFrame
import config as cfg
import sys
import brain as bd

if __name__ == '__main__':
    mobidb_results_df = pd.read_csv(cfg.data[''] + '/mobidb_result.tsv', sep='\t')
    ## brain proteins
    brain_prot_lst = bd.brain_pr_lst_generator()  # n: 8428
    brain_pr_df = DataFrame(brain_prot_lst, columns=['acc'])
    brain_pr_df = brain_pr_df.dropna()  # 8337
    brain_pr_df['brain'] = 1
    mobi_brain_bool_df = mobidb_results_df.merge(brain_pr_df, how='left', on='acc')
    mobi_brain_bool_df['brain'] = mobi_brain_bool_df['brain'].fillna(0)  # (1212280, 7)

    ## NDD proteins
    positive_cand_g4mobi_df = pd.read_csv(cfg.data['gene4'] + '/positive_cand_g4mobi_concat.csv')
    ndd_pos_cand_acc_lst = positive_cand_g4mobi_df['acc_x'].unique().tolist()  # 177
    ndd_pos_cand_acc_df = DataFrame(ndd_pos_cand_acc_lst, columns=['acc'])
    ndd_pos_cand_acc_df = ndd_pos_cand_acc_df.drop_duplicates()
    ndd_pos_cand_acc_df['ndd'] = 1
    mobi_ndd_bool_df = mobi_brain_bool_df.merge(ndd_pos_cand_acc_df, how='left', on='acc')
    mobi_ndd_bool_df['ndd'] = mobi_ndd_bool_df['ndd'].fillna(0)  # (1212289, 8)

    ## Phenotypes
    # 163
    asd_pos_cand_acc_lst = positive_cand_g4mobi_df.loc[positive_cand_g4mobi_df['Phenotype']=='ASD', 'acc_x'].unique().tolist()
    asd_pos_cand_acc_df = DataFrame(asd_pos_cand_acc_lst, columns=['acc'])
    asd_pos_cand_acc_df['asd'] = 1
    mobi_asd_bool_df = mobi_ndd_bool_df.merge(asd_pos_cand_acc_df, how='left', on='acc')
    mobi_asd_bool_df['asd'] = mobi_asd_bool_df['asd'].fillna(0)  # (1212289, 9)
    print("hi")

    sys.exit()
