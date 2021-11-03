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


def var_in_ptm_checker(input_df, source):
    if source == 'mobidb':
        return mobidb_var_in_ptm_checker(input_df)
    elif source == 'uniprot':
        return uniprot_var_in_ptm_checker(input_df)


def dismaj_var_in_ptm_df_generator():
    # this is just to organize the code, disorder_majority and ptms df are being merged and then with
    # var_in_ptm_checker a new df is produce that shows var positions that are in ptm sites (ptm pos)
    ptms_df = pd.read_csv(cfg.data['ptm-u'] + '/uniprot-ptms-all.csv',
                          usecols=['acc', 'ptm_pos', 'description', 'ptm_type'])  # (140538, 4)
    disorder_maj = pd.read_csv(cfg.data['vars'] + '/disorder-majority-inout-idr-vars-count-normalized.csv', usecols=
    ['acc', 'var_id', 'orig_aa', 'var_aa', 'position', 'isin_idr', 'total_vars', 'content_count'])  # (66843, 7)
    dismaj_ptm_df = pd.merge(disorder_maj, ptms_df, on='acc')
    ptm_checked_dismaj_df = var_in_ptm_checker(dismaj_ptm_df, 'uniprot')
    ptm_checked_dismaj_df.to_csv(cfg.data['ptm-u'] + '/uniprot-vars-inptm-checked.csv')
    return ptm_checked_dismaj_df


def ndd_idrvar_in_ptm_lst_df_generator(all_ptm_checked_pr_lst):
    ndd_subdf = pd.read_csv(cfg.data['phens-fdr'] + '/acc-phen-5percentFDR.csv')
    # phens_lst = ['ASD', 'EE', 'ID', 'DD', 'SCZ', 'NDDs', 'Control']
    # ndd_subdf = ndd_subdf.loc[ndd_subdf.Phenotype.isin(phens_lst)]
    ndd_pr_lst = ndd_subdf['acc'].unique().tolist()  # 1308 proteins
    ptm_idr_var_pr_lst = list(set(ndd_pr_lst).intersection(all_ptm_checked_pr_lst))
    ptm_idr_var_ndd_subdf = ndd_subdf[ndd_subdf.acc.isin(ptm_idr_var_pr_lst)]
    return ptm_idr_var_pr_lst, ptm_idr_var_ndd_subdf


if __name__ == '__main__':
    # var_in_ptm_checked_df = dismaj_var_in_ptm_df_generator()
    var_in_ptm_checked_df = pd.read_csv(cfg.data['ptm-u'] + '/uniprot-vars-inptm-checked.csv')
    # var_in_ptm_df = var_in_ptm_checked_df.loc[var_in_ptm_checked_df['var_in_ptm'] == 1]  # (3170,12) vars, 527 Prs
    ptm_idr_var_all_df = var_in_ptm_checked_df.loc[(var_in_ptm_checked_df['var_in_ptm'] == 1) &
                                                   (var_in_ptm_checked_df['isin_idr'] == 1)]  # (366, 13)
    ptm_idr_var_all_pr_lst = ptm_idr_var_all_df['acc'].unique().tolist()  # 177 Prs.
    ## for NDDs
    ptm_idr_var_ndd_pr_lst, _ = ndd_idrvar_in_ptm_lst_df_generator(ptm_idr_var_all_pr_lst)  # 16 contributes to 77 rows
    # in ndd phens col, meaning each pr is in charge of ~ 5 phens among all phens and not just my desired phenotypes
    print('\n'.join(ptm_idr_var_ndd_pr_lst))
    ptm_type_count = ptm_idr_var_all_df.groupby('ptm_type').count()

    def ptm_type_seperator(inputdf):
        ptms_disulfide_bonds_df = inputdf.loc[inputdf['ptm_type'] == 'disulfide bond']
        the_rest_of_ptms_df = inputdf.loc[inputdf['ptm_type'] != 'disulfide bond']
        return ptms_disulfide_bonds_df, the_rest_of_ptms_df
