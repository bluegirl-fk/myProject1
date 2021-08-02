# this file is supposed to get mobidb+mutposition, is_in_startend and ndd bool, index*, and give a df in the outcome
# with added columns: brain bool, ASD bool, EE bool, ID bool, SCZ bool
# *index (mutual idx of mut_acc_merged_df and mutinfo subdf which is useful to merge mobidb & g4dn_mutinfo_acc,
# could be useless now :/ not sure

import pandas as pd
import config as cfg
import sys
import brain as bd


if __name__ == '__main__':
    # this file is generated in categorizer.py with mobi_mut_inidr_checker()
    mobi_ndd_idr_df = pd.read_csv(cfg.data['gene4'] + '/mut-pos-mobi1.csv')
    mobi_ndd_idr_df = mobi_ndd_idr_df.drop(['Unnamed: 0', 'acc.1'], axis=1)  # (1212280, 10)
    ## brain proteins
    brain_prot_lst = bd.brain_pr_lst_generator()  # n: 8428

    sys.exit()
