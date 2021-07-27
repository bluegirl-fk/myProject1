# this file is supposed to import the mutations inside and outside IDRs separately and perform different operations:
# a df with mutations filtered based on mobidb-lite feature
# disorder residue length + complement
# July 27th
import pandas as pd
import config as cfg
import sys


def main():
    print("whatever")
    return


if __name__ == '__main__':

    mobi_mut_in_idr_df = pd.read_csv(cfg.data['gene4'] + '/mobidb-mut-pos-true.csv')  # (994383, 11)
    mobi_mut_out_idr_df = pd.read_csv(cfg.data['gene4'] + '/mobidb-mut-pos-false.csv')  # (3238392, 11)

    # (25961, 11)
    mut_in_mobi_lite_df = mobi_mut_in_idr_df.loc[mobi_mut_in_idr_df['feature'] == 'prediction-disorder-mobidb_lite']
    # (102495, 11)
    mut_out_mobi_lite_df = mobi_mut_out_idr_df.loc[mobi_mut_out_idr_df['feature'] == 'prediction-disorder-mobidb_lite']

    ## merge with g4dn
    mut_acc_mrg_df = pd.read_csv(cfg.data['gene4'] + '/mut-acc-mrg-df.csv', low_memory=False)
    in_mobi_g4dn_df = pd.merge(mut_in_mobi_lite_df, mut_acc_mrg_df, on='acc')
    # out_mobi_g4dn_df = pd.merge(mut_out_mobi_lite_df, mut_acc_mrg_df, on='acc')

# TODO: deal with this error: Process finished with exit code 137 (interrupted by signal 9: SIGKILL)
    sys.exit(main())
