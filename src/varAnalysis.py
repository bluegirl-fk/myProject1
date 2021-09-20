# this is to perform previous analysis which was done in mobibool and protanalyse but for 'variants in IDR' dataset
# Sep 17th
import pandas as pd
from pandas import DataFrame
import config as cfg
import brain as brn


def vars_df_generator(input_df):
    var_pr_lst = input_df['acc'].unique().tolist()
    ## ndds in variants
    ndd_in_idr = list(set(ndd_pr_lst).intersection(var_pr_lst))
    ndd_in_idr_subdf = ndd_subdf[ndd_subdf.acc.isin(ndd_in_idr)]
    del ndd_in_idr_subdf['Unnamed: 0']
    ## Brain
    brn_in_idr_lst = list(set(brain_prot_lst).intersection(var_pr_lst))
    brain_inidr_subdf = DataFrame(brn_in_idr_lst, columns=['acc'])

    return input_df, brain_inidr_subdf, ndd_in_idr_subdf


def all_vars_or_vars_inidr(input):
    vars_in_idr_df = pd.read_csv(cfg.data['vars'] + '/all-vars-all-features.csv')
    all_vars = pd.read_csv(cfg.data['vars'] + '/all-variants-mobidb-merged.csv')
    del all_vars['Unnamed: 0_x']
    del all_vars['Unnamed: 0_y']
    if input == 'all':
        return vars_df_generator(all_vars)
    elif input == 'idr':
        return vars_df_generator(vars_in_idr_df)


## NDD and brain original import
ndd_subdf = pd.read_csv(cfg.data['phens-fdr'] + '/acc-phen-5percentFDR.csv')
ndd_subdf = ndd_subdf.drop_duplicates()  # (4531, 3)
ndd_pr_lst = ndd_subdf['acc'].unique().tolist()  # 1308 proteins
brain_prot_lst = brn.brain_pr_lst_generator()  # n: 8320






