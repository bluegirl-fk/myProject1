# generated on November 23rd to have the previous varptm class more organized
import pandas as pd
import config as cfg
import numpy as np


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
            # here we fist make a list of ptm_pos with range of 3 AAs before and after the ptm-point (not for disulfide)
            ptm_pos_tmp_lst = []
            ptm_pos_tmp_lst = np.arange(int(row.ptm_pos) - 3, int(row.ptm_pos) + 4)
            # (delete negative numbers later)
            if row.position in ptm_pos_tmp_lst:
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


def ndd_idrvar_in_ptm_lst_df_generator(all_ptm_checked_pr_lst, all_ptm_var_filtered_df):
    ndd_subdf = pd.read_csv(cfg.data['phens-fdr'] + '/acc-phen-5percentFDR.csv')
    # phens_lst = ['ASD', 'EE', 'ID', 'DD', 'SCZ', 'NDDs', 'Control']
    # ndd_subdf = ndd_subdf.loc[ndd_subdf.Phenotype.isin(phens_lst)]
    ndd_pr_lst = ndd_subdf['acc'].unique().tolist()  # 1308 proteins
    ptm_idr_var_pr_lst = list(set(ndd_pr_lst).intersection(all_ptm_checked_pr_lst))
    ptm_idr_var_ndd_subdf = ndd_subdf[ndd_subdf.acc.isin(ptm_idr_var_pr_lst)]
    ptm_idr_var_ndd_df = all_ptm_var_filtered_df.loc[all_ptm_var_filtered_df.acc.isin(ndd_pr_lst)]
    return ptm_idr_var_pr_lst, ptm_idr_var_ndd_subdf, ptm_idr_var_ndd_df


def ptm_divider(inputdf):
    disulfide_bonds_df = inputdf.loc[inputdf['ptm_type'] == 'disulfide bond']
    the_rest_of_ptms_df = inputdf.loc[inputdf['ptm_type'] != 'disulfide bond']
    disulfide_bonds_pr_lst = disulfide_bonds_df['acc'].unique().tolist()
    the_rest_of_ptms_pr_lst = the_rest_of_ptms_df['acc'].unique().tolist()
    return disulfide_bonds_df, the_rest_of_ptms_df, disulfide_bonds_pr_lst, the_rest_of_ptms_pr_lst


## from here on, the methods are for checking if PTM is in disorder
# (previous one was giving away var in ptm, and isin_idr separately)


def expand_regions(region_ranges_lst):
    transformed_regions = []
    region_ranges_lst = region_ranges_lst.split(',')
    # region_ranges_lst = list(region_ranges_lst)
    for region in region_ranges_lst:
        start = int(region.split('..')[0])
        end = int(region.split('..')[1])
        while start <= end:
            transformed_regions.append(start)
            start += 1
    return set(transformed_regions)


def ptmIDR_bool_array_maker():
    ## checks if ptm position is in startend disorder region of mobidb or not
    # makke merged df to be checked for ptm in disorder
    ptms_df = pd.read_csv(cfg.data['ptm-u'] + '/uniprot-ptms-all.csv',
                          usecols=['acc', 'ptm_pos', 'description', 'ptm_type'])  # (140538, 4)
    disorder_maj = pd.read_csv(cfg.data['vars'] + '/disorder-majority-inout-idr-vars-count-normalized.csv', usecols=
    ['acc', 'var_id', 'orig_aa', 'var_aa', 'position', 'isin_idr', 'total_vars', 'content_count', 'startend'])
    dismaj_ptm_df = pd.merge(disorder_maj, ptms_df, on='acc')
    # starts here
    array_is_in = []
    for index, row in dismaj_ptm_df.iterrows():
        if '&' not in row.ptm_pos:
            print('one')
            lst_disorder_region = list(expand_regions(row.startend))
            if int(row.ptm_pos) in lst_disorder_region:
                array_is_in.append('1')
            else:
                array_is_in.append('0')
        elif '&' in row.ptm_pos:
            print('two')
            lst_ptm_pos = list(set(str(row.ptm_pos).split('&')))
            lst_disorder_region = list(expand_regions(row.startend))
            lst_disorder_region = [str(i) for i in lst_disorder_region]
            if len(list(set(lst_ptm_pos) & set(lst_disorder_region))) >= 1:  # at least have one value in common
                array_is_in.append('1')
            else:
                array_is_in.append('0')

    dismaj_ptm_df['ptm_inidr'] = array_is_in
    return dismaj_ptm_df


if __name__ == '__main__':
    ptms_df = pd.read_csv(cfg.data['ptm-u'] + '/uniprot-ptms-all.csv',
                          usecols=['acc', 'ptm_pos', 'description', 'ptm_type'])  # (140538, 4)
    ptms_all_lst = ptms_df['acc'].unique().tolist()
    ## PTM in disorder
    ptm_in_idr_checked_df = pd.read_csv(cfg.data['ptm'] + '/ptm_in_idr_checked_(uniprot)-with-disulfide.csv')
    ptm_inidr_df = ptm_in_idr_checked_df.loc[ptm_in_idr_checked_df['ptm_inidr'] == 1]  # 362784
    ## PTM site +-3 in var position
    ptm_var_checked_df = pd.read_csv(cfg.data['ptm-u'] + '/uniprot-vars-inptm-checked+-3res.csv')  # (1392057, 13)
    var_in_ptm_df = ptm_var_checked_df.loc[ptm_var_checked_df['var_in_ptm'] == 1]  # (6419, 13) vars, 527 Prs
    ## merged PTM in disorder in var
    ptm_in_idr_in_var_mrg = pd.merge(ptm_inidr_df, var_in_ptm_df[['var_id', 'var_in_ptm']], on='var_id')
    ptm_idr_var_lst = ptm_in_idr_in_var_mrg['acc'].unique().tolist()  # 1288
    ## divided ptms into disulfide bonds and the rest of ptm types
    _, _, ss_bond_pr_lst, other_ptms_pr_lst = ptm_divider(ptm_in_idr_in_var_mrg)  # 165 # 1234

    ### NDDs
    ndd_subdf = pd.read_csv(cfg.data['phens-fdr'] + '/acc-phen-5percentFDR.csv')
    ndd_pr_lst = ndd_subdf['acc'].unique().tolist()  # 1308 proteins
    ## PTM in disorder
    ndd_ptm_idr_checked_df = ptm_in_idr_checked_df.loc[ptm_in_idr_checked_df.acc.isin(ndd_pr_lst)]  # (303357, 14)
    ndd_ptm_inidr_df = ndd_ptm_idr_checked_df.loc[ndd_ptm_idr_checked_df['ptm_inidr'] == 1]  # (67608)
    ## PTM site +-3 in var position
    ndd_ptm_var_checked_df = ptm_var_checked_df.loc[ptm_var_checked_df.acc.isin(ndd_pr_lst)]  # (303357, 13)
    ndd_var_in_ptm_df = ndd_ptm_var_checked_df.loc[ndd_ptm_var_checked_df['var_in_ptm'] == 1]  # (1244, 12)
    ## merged PTM in disorder in var
    ndd_ptm_idr_var_mrg = pd.merge(ndd_ptm_inidr_df, ndd_var_in_ptm_df[['var_id', 'var_in_ptm']], on='var_id')
    ndd_ptm_idr_var_lst = ndd_ptm_idr_var_mrg['acc'].unique().tolist()  # 121
    ## divided ptms into disulfide bonds and the rest of ptm types
    _, _, ndd_ss_bond_pr_lst, ndd_other_ptms_pr_lst = ptm_divider(ndd_ptm_idr_var_mrg)  # 9 # 119

    ## the ptm division could be performed during other states of analysis as well, shown always in a chart


