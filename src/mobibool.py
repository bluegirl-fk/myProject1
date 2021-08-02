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
    # this file is generated in categorizer.py with mobi_mut_inidr_checker()
    mobi_ndd_idr_df = pd.read_csv(cfg.data['gene4'] + '/mut-pos-mobi2000.csv', low_memory=False)
    mobi_ndd_idr_df = mobi_ndd_idr_df.drop(['Unnamed: 0', 'acc.1'], axis=1)  # (1212280, 10)
    ## brain proteins
    brain_prot_lst = bd.brain_pr_lst_generator()  # n: 8428
    brain_pr_df = DataFrame(brain_prot_lst, columns=['brain_acc'])
    brain_pr_df = brain_pr_df.dropna()  # 8337
    brain_pr_df['brain'] = 1
    mobi_brain_bool_df = pd.concat([mobi_ndd_idr_df, brain_pr_df], axis=1)
    mobi_brain_bool_df['brain'] = mobi_brain_bool_df['brain'].fillna(0)

    sys.exit()
