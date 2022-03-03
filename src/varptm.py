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
    disulfide_bonds_vid_lst = disulfide_bonds_df['var_id'].unique().tolist()
    the_rest_of_ptms_vid_lst = the_rest_of_ptms_df['var_id'].unique().tolist()
    return disulfide_bonds_df, the_rest_of_ptms_df, disulfide_bonds_vid_lst, the_rest_of_ptms_vid_lst


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
    # var_in_ptm_checked_df = dismaj_var_in_ptm_df_generator()
    # var_in_ptm_checked_df.to_csv(cfg.data['ptm-u'] + '/uniprot-vars-inptm-checked+-3res.csv')
    var_in_ptm_checked_df = pd.read_csv(cfg.data['ptm-u'] + '/uniprot-vars-inptm-checked+-3res.csv')  # (1392057, 13)
    var_in_ptm_df = var_in_ptm_checked_df.loc[var_in_ptm_checked_df['var_in_ptm'] == 1]  # (6419, 13) vars, 527 Prs
    inptm_idr_var_all_df = var_in_ptm_checked_df.loc[(var_in_ptm_checked_df['var_in_ptm'] == 1) &
                                                     (var_in_ptm_checked_df['isin_idr'] == 1)]  # (2126, 13)
    # (find also ptms in disorder, without considering variations)
    inptm_idr_var_all_pr_lst = inptm_idr_var_all_df['acc'].unique().tolist()  # 808 Prs.
    _, _, all_disulfide_vid_lst, all_other_ptms_vid_lst = ptm_divider(inptm_idr_var_all_df)  # 35 # 1582
    ptm_type_count = inptm_idr_var_all_df.groupby('ptm_type').count()
    ## for NDDs
    in_ptm_idr_var_ndd_pr_lst, _, inptm_idr_var_ndd_df = ndd_idrvar_in_ptm_lst_df_generator \
        (inptm_idr_var_all_pr_lst, inptm_idr_var_all_df)  # 76 # (365, 13)
    # contributes to 77 rows in phens col, so each pr is in charge of ~5 phens among all phens and not just my phens
    # print('\n'.join(in_ptm_idr_var_ndd_pr_lst))
    _, _, ndd_disulfide_vid_lst, ndd_other_ptms_vid_lst = ptm_divider(inptm_idr_var_ndd_df)  # 2  # 239

    ## NDD variants in ptm, even if not in disordered regions, can regulate IDP activation
    ndd_subdf = pd.read_csv(cfg.data['phens-fdr'] + '/acc-phen-5percentFDR.csv')
    ndd_pr_lst = ndd_subdf['acc'].unique().tolist()  # 1308 proteins
    ndd_vars_in_ptm = var_in_ptm_checked_df.loc[(var_in_ptm_checked_df.acc.isin(ndd_pr_lst))
                                                & (var_in_ptm_checked_df['var_in_ptm'] == 1)]  # (1244, 13)
    ndd_vars_in_ptm_lst = ndd_vars_in_ptm['acc'].unique().tolist()  # 127
    # print('\n'.join(ndd_other_ptms_pr_lst))

    # cm: maybe find ptms and disease (vars), and checked the relation with disorder and LLPS
    disorder_maj = pd.read_csv(cfg.data['vars'] + '/disorder-majority-inout-idr-vars-count-normalized.csv', usecols=
    ['acc', 'var_id', 'orig_aa', 'var_aa', 'position', 'isin_idr', 'total_vars', 'content_count', 'startend'])
    ptms_df = pd.read_csv(cfg.data['ptm-u'] + '/uniprot-ptms-all.csv',
                          usecols=['acc', 'ptm_pos', 'description', 'ptm_type'])  # (140538, 4)
    ptms_all_lst = ptms_df['acc'].unique().tolist()
    # ptm_in_idr_checked_df = ptmIDR_bool_array_maker()
    ptm_in_idr_checked_df = pd.read_csv(cfg.data['ptm'] + '/ptm_in_idr_checked_(uniprot)-with-disulfide.csv')
    ptm_in_disorder_df = ptm_in_idr_checked_df.loc[ptm_in_idr_checked_df['ptm_inidr'] == 1]  # 362784
    ptm_in_idr_in_var_mrg = pd.merge(ptm_in_disorder_df, var_in_ptm_df[['var_id', 'var_in_ptm']], on='var_id')
    lst_a = ptm_in_idr_in_var_mrg['acc'].unique().tolist()  # 1288
    # maybe make a new one

    disulfide_in_disorder = ptm_in_disorder_df.loc[ptm_in_disorder_df['ptm_type'] == 'disulfide bond']
    disulfide_in_disorder_pr_lst = disulfide_in_disorder['acc'].unique().tolist()  # 524 proteins
    # remember that you are using filtered dismaj based on cc > 20 residues
    # todo: retry with normal disorder majority, also organize the code and method()
    # also this part should come at first, so first ptm in disorder, then var in ptm

    ## list of Vars in PTMs
    # a = var_in_ptm_checked_df.loc[var_in_ptm_checked_df['var_in_ptm'] == 1]
    # print('\n'.join((a['var_id']).unique().tolist()))
    # ## list of Vars in IDR
    # b = var_in_ptm_checked_df.loc[var_in_ptm_checked_df['isin_idr'] == 1]
    # print('\n'.join((b['var_id']).unique().tolist()))
    ## Variants in IDR in PTM
    # c = var_in_ptm_checked_df.loc[(var_in_ptm_checked_df['isin_idr'] == 1) & (var_in_ptm_checked_df['var_in_ptm']==1)]
    # print('\n'.join((c['var_id']).unique().tolist()))
    ## ndd variants
    ndd_variants = var_in_ptm_checked_df.loc[var_in_ptm_checked_df.acc.isin(ndd_pr_lst)]
    print('\n'.join((ndd_variants['var_id']).unique().tolist()))
    # print(len((ndd_variants['var_id']).unique().tolist()))
