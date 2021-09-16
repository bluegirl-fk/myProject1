#  merge with mobidb here
import pandas as pd
import config as cfg



def expand_regions(region_ranges_lst):
    transformed_regions = []
    for reg in region_ranges_lst:
        start = int(reg.split('..')[0])
        end = int(reg.split('..')[1])
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

array = mutidr_bool_array_maker(mobidb_lite_mrg)

## add bool array to the df
mobidb_lite_mrg['is_in_startend'] = array

