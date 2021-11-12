# checking other variants from uniprot Nov 12th
import pandas as pd
import config as cfg


def humsavar_pos_adder():
    hsvar = pd.read_csv(cfg.data[''] + '/humsavar.txt', delimiter=r"\s+",
                        usecols=['genename', 'AC', 'FTId', 'change', 'category', 'dbSNP', 'Diseasename'])
    # fix proteins that don't have gene name and acc and other cols right in place
    hsvar['change'] = hsvar.change.str.split(pat='p.').str[1]
    hsvar['orig_aa'] = hsvar.change.str[:3]
    hsvar['var_aa'] = hsvar.change.str[-3:]
    hsvar['position'] = hsvar.change.str[3:-3]

    return hsvar


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
    return list(set(transformed_regions))


def mutidr_bool_array_maker(input_df):
    # input df is merged df of mobidb and mutation positions from mutinfo
    ## checks if mutation position is in startend disorder region of mobidb or not
    array_is_in = []  # will be filled with boolean of 0,1
    for index, row in input_df.iterrows():
        disorder_region_lst = []
        disorder_region_lst = expand_regions(row.startend)  # temp set of data, convert each startend lst to a set,
        # write in report
        if row.position in disorder_region_lst:
            array_is_in.append('1')
        else:
            array_is_in.append('0')
    return array_is_in


def isin_idr_col_adder(df):
    # this gets list of dataframes to be check if they have mut positions in IDR and adds a bool array to them
    array = mutidr_bool_array_maker(df)
    ## add bool array to the df
    df['isin_idr'] = array
    return df


if __name__ == '__main__':
    hmvar = humsavar_pos_adder()
    hmvar.to_csv(cfg.data[''] + '/new-humsavar.csv')
    mobidb = pd.read_csv(cfg.data['phens'] + '/mobidb-results+ndd-tsv-damiano-shared.tsv')
    mobidb = mobidb.rename(columns={'start..end' : 'startend'})
    # hmvar_mrg = pd.merge(mobidb, hmvar, left_on='acc', right_on='AC')
    hmvar_mrg = pd.read_csv(cfg.data['vars'] + '/humsavars-mobidb-merged.csv')

    ## adding isin_idr bool array
    # vars_checked_df = isin_idr_col_adder(hmvar_mrg)
    vars_checked_df = pd.read_csv(cfg.data['vars'] + '/humsavars-mobidb-merged-checked.csv')
    vars_mobilite_df = vars_checked_df.loc[vars_checked_df['feature'] == 'prediction-disorder-mobidb_lite']
    vars_mobilite_inidr_df = vars_mobilite_df.loc[vars_mobilite_df['isin_idr'] == 1]

    # todo categorize based on benign or pathologic