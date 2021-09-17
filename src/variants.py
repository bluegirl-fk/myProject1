#  merge with mobidb here
import pandas as pd
import config as cfg
import re


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


def mutidr_bool_array_maker(input_df):
    # input df is merged df of mobidb and mutation positions from mutinfo
    ## checks if mutation position is in startend disorder region of mobidb or not
    array_is_in = []  # will be filled with boolean of 0,1
    for index, row in input_df.iterrows():
        set_disorder_region = expand_regions(row.startend)  # temp set of data, convert each startend lst to a set,
        # write in report
        if row.position in set_disorder_region:
            array_is_in.append('1')
        else:
            array_is_in.append('0')
    return array_is_in


def isin_idr_col_adder(dfs_lst):
    # this gets list of dataframes to be check if they have mut positions in IDR and adds a bool array to them
    for df in dfs_lst:
        array = mutidr_bool_array_maker(df)
        ## add bool array to the df
        df['isin_idr'] = array
    return dfs_lst


def in_idr_df_generator(dfs_lst):
    # generates df with proteins that have variant in IDR, can get as input list of df(s)
    df_only_muinIDR_lst = []
    for df in dfs_lst:
        new_df = df.loc[df['isin_idr'] == '1']
        new_df = new_df.drop(columns=['Unnamed: 0_x', 'Unnamed: 0_y', 'isin_idr'])
        df_only_muinIDR_lst.append(new_df)
    return df_only_muinIDR_lst


## mobidb, filtered based on mobidblite and disorder majority
mobidb = pd.read_csv(cfg.data['phens'] + '/mobidb-results+ndd-tsv-damiano-shared.tsv')
mobidb = mobidb.rename(columns={'start..end': 'startend'})
mobidb_lite = mobidb.loc[mobidb['feature'] == 'prediction-disorder-mobidb_lite']
disorder_majority = mobidb.loc[mobidb['feature'] == 'prediction-disorder-th_50']
## Variants df, (12 proteins not in mobidb)
var_prs_df = pd.read_csv(cfg.data['xml-p'] + '/protein-vars.csv')
## Merged DFs
# mobidb_lite
mobidb_lite_mrg = pd.merge(mobidb_lite, var_prs_df, on='acc')
# Disorder_majority
disorder_majority_mrg = pd.merge(disorder_majority, var_prs_df, on='acc')
## adding isin_idr bool array
[mobidb_lite_mrg, disorder_majority_mrg] = isin_idr_col_adder([mobidb_lite_mrg, disorder_majority_mrg])
## DFs with muts only inside IDR
lite_mut_in_df, dismaj_mut_in_df = in_idr_df_generator([mobidb_lite_mrg, disorder_majority_mrg])
## to CSV
lite_mut_in_df.to_csv(cfg.data['vars'] + '/muts-in-IDR-mobidb_lite.csv')
dismaj_mut_in_df.to_csv(cfg.data['vars'] + '/muts-in-IDR-disorder_majority.csv')

# TODO: insert each of this in mobibool instead of mobidb, but you need to also delete ndd and brain that are not in
#  this list maybe use left merge!
