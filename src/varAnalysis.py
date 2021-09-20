# this is to perform previous analysis which was done in mobibool and protanalyse but for 'variants in IDR' dataset
# Sep 17th
import pandas as pd
from pandas import DataFrame
import config as cfg
import brain as brn
# TODO: insert each of this in mobibool instead of mobidb, but you need to also delete ndd and brain that are not in
#  this list maybe use left merge!


vars_in_idr_df = pd.read_csv(cfg.data['vars'] + '/all-vars-all-features.csv')
## NDD
ndd_subdf = pd.read_csv(cfg.data['phens-fdr'] + '/acc-phen-5percentFDR.csv')
ndd_subdf = ndd_subdf.drop_duplicates()  # (4531, 3)
ndd_pr_lst = ndd_subdf['acc'].unique().tolist()  # 1308 proteins
## ndds in variants
ndd_in_idr = list(set(ndd_pr_lst).intersection(vars_in_idr_df))
## Brain
brain_prot_lst = brn.brain_pr_lst_generator()  # n: 8320
brn_in_idr = list(set(brain_prot_lst).intersection(vars_in_idr_df))



