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


## mobidb, filtered based on mobidblite and disorder majority
mobidb = pd.read_csv(cfg.data['phens'] + '/mobidb-results+ndd-tsv-damiano-shared.tsv',
                     usecols=['acc', 'feature', 'start..end'])
mobidb = mobidb.rename(columns={'start..end': 'startend'})
mobidb_lite = mobidb.loc[mobidb['feature'] == 'prediction-disorder-mobidb_lite']
del mobidb_lite['feature']
disorder_majority = mobidb.loc[mobidb['feature'] == 'prediction-disorder-th_50']
del disorder_majority['feature']
## Variants df, (12 proteins not in mobidb)
var_prs_df = pd.read_csv(cfg.data['xml-p'] + '/protein-vars.csv')
## Merged DFs
# mobidb_lite
mobidb_lite_mrg = pd.merge(mobidb_lite, var_prs_df, on='acc')
del mobidb_lite_mrg['Unnamed: 0']
# Disorder_majority
disorder_majority_mrg = pd.merge(disorder_majority, var_prs_df, on='acc')
del disorder_majority_mrg['Unnamed: 0']
## adding isin_idr bool array
[mobidb_lite_mrg, disorder_majority_mrg] = isin_idr_col_adder([mobidb_lite_mrg, disorder_majority_mrg])

lite_mut_in_df = mobidb_lite_mrg.loc[mobidb_lite_mrg['isin_idr'] == '1']
dismaj_mut_in_df = disorder_majority_mrg.loc[disorder_majority_mrg['isin_idr'] == '1']
