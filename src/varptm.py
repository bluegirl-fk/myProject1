import pandas as pd
import config as cfg


def var_in_ptm_checker(input_df):
    # gets df as input, checkes if variation position is in ptm position or not, adds column of booleans
    array_is_in = []  # will be filled with boolean of 0,1
    for index, row in input_df.iterrows():
        if row.ptm_pos is not None:
            set_ptm_pos = set((str(row.ptm_pos)).split(','))
            if str(row.position) in set_ptm_pos:
                array_is_in.append('1')
            else:
                array_is_in.append('0')
        else:
            array_is_in.append('0')
    input_df['ptm_mut'] = array_is_in
    return input_df


if __name__ == '__main__':
    ptms_df = pd.read_csv(cfg.data['ptm'] + '/mobidb-parsed-ptm-info.csv', usecols=['acc', 'ptm_pos'])
    disorder_maj = pd.read_csv(cfg.data['vars'] + '/disorder-majority-inout-idr-vars-count-normalized.csv')
    dismaj_ptm_df = pd.merge(disorder_maj, ptms_df, on='acc')
    del dismaj_ptm_df['Unnamed: 0']
    dismaj_ptm_df = dismaj_ptm_df.fillna(0)
    ## Df with PTM sites checked
    dismaj_ptmsites_df = var_in_ptm_checker(dismaj_ptm_df)
    dismaj_ptmsites_df.to_csv(cfg.data['ptm'] + '/variants-in-ptm-sites.csv')
    var_in_ptm_df = dismaj_ptmsites_df.loc[dismaj_ptmsites_df['ptm_mut'] == '1']
    var_in_ptm_lst = var_in_ptm_df['acc'].unique().tolist()

    ## NDDs with variation in PTM sites: n= 35
    ndd_subdf = pd.read_csv(cfg.data['phens-fdr'] + '/acc-phen-5percentFDR.csv')
    ndd_subdf = ndd_subdf.drop_duplicates()  # (4531, 3)
    ndd_pr_lst = ndd_subdf['acc'].unique().tolist()  # 1308 proteins
    ndd_var_in_ptm_lst = list(set(ndd_pr_lst).intersection(var_in_ptm_lst))
    ndd_var_in_ptm_subdf = ndd_subdf[ndd_subdf.acc.isin(ndd_var_in_ptm_lst)]
