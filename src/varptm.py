import pandas as pd
import config as cfg


def mobidb_var_in_ptm_checker(df):  # input is merged 1-variants+in/out disorder columns and 2-ptms df from mobidb
    # gets df as input, checkes if variation position is in ptm position or not, adds column of booleans
    array_is_in = []  # will be filled with boolean of 0,1
    for index, row in df.iterrows():
        if row.ptm_pos is not None:
            set_ptm_pos = set((str(row.ptm_pos)).split(','))
            if str(row.position) in set_ptm_pos:
                array_is_in.append('1')
            else:
                array_is_in.append('0')
        else:
            array_is_in.append('0')
    df['var_in_ptm'] = array_is_in
    return df


def uniprot_var_in_ptm_checker(df):  # input is merged 1-variants+in/out disorder columns and 2-ptms df from uniprot
    array_isin = []
    for index, row in df.iterrows():
        if '&' not in row.ptm_pos:
            if str(row.position) == str(row.ptm_pos):
                array_isin.append('1')
            else:
                array_isin.append('0')
        elif '&' in row.ptm_pos:
            set_ptm_pos = set(str(row.ptm_pos).split('&'))
            if str(row.position) in set_ptm_pos:
                array_isin.append('1')
            else:
                array_isin.append('0')
    df['var_in_ptm'] = array_isin
    return df


def var_in_ptm_checker(input_df, resource):
    if resource == 'mobidb':
        return mobidb_var_in_ptm_checker(input_df)
    elif resource == 'uniprot':
        return uniprot_var_in_ptm_checker(input_df)


if __name__ == '__main__':
    ptms_df = pd.read_csv(cfg.data['ptm-u'] + '/uniprot-ptms-all.csv', usecols=['acc', 'ptm_pos', 'ptm_type'])
    disorder_maj = pd.read_csv(cfg.data['vars'] + '/disorder-majority-inout-idr-vars-count-normalized.csv', usecols=
    ['acc', 'var_id', 'orig_aa', 'var_aa', 'position', 'isin_idr', 'total_vars', 'in_idr_vars', 'out_idr_vars'])
    dismaj_ptm_df = pd.merge(disorder_maj, ptms_df, on='acc')

    # ptm_checked_dismaj_df = var_in_ptm_checker(dismaj_ptm_df, 'uniprot')
    # ptm_checked_dismaj_df.to_csv(cfg.data['ptm-u'] + '/uniprot-vars-inptm-checked.csv')
    ptm_checked_dismaj_df = pd.read_csv(cfg.data['ptm-u'] + '/uniprot-vars-inptm-checked.csv')
    del ptm_checked_dismaj_df['Unnamed: 0']
    var_in_ptm_df = ptm_checked_dismaj_df.loc[ptm_checked_dismaj_df['var_in_ptm'] == 1]  # why the 1 is not string!
    var_in_ptm_lst = var_in_ptm_df['acc'].unique().tolist()

    ## NDDs with variation in PTM sites: n= 35
    ndd_subdf = pd.read_csv(cfg.data['phens-fdr'] + '/acc-phen-5percentFDR.csv')
    ndd_subdf = ndd_subdf.drop_duplicates()  # (4531, 3)
    ndd_pr_lst = ndd_subdf['acc'].unique().tolist()  # 1308 proteins
    ndd_var_in_ptm_lst = list(set(ndd_pr_lst).intersection(var_in_ptm_lst))  # n:49
    ndd_var_in_ptm_subdf = ndd_subdf[ndd_subdf.acc.isin(ndd_var_in_ptm_lst)]  # (191,3)
