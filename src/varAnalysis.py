# this is to perform previous analysis which was done in mobibool and protanalyse but for 'variants in IDR' dataset
# Sep 17th
import pandas as pd
from pandas import DataFrame
import config as cfg
import brain as brn
# TODO: insert each of this in mobibool instead of mobidb, but you need to also delete ndd and brain that are not in
#  this list maybe use left merge!

mobidblite_inidr_df = pd.read_csv(cfg.data['vars'] + '/muts-in-IDR-mobidb_lite.csv')
dismajority_inidr_df = pd.read_csv(cfg.data['vars'] + '/muts-in-IDR-disorder_majority.csv')
variants_df = pd.read_csv(cfg.data['xml-p'] + '/protein-vars.csv')
variants_lst = variants_df['acc'].unique().tolist()
# NDD
ndd_subdf = pd.read_csv(cfg.data['phens-fdr'] + '/acc-phen-5percentFDR.csv')
ndd_subdf = ndd_subdf.drop_duplicates()  # (4531, 3)
ndd_pr_lst = ndd_subdf['acc'].unique().tolist()  # 1308 proteins
## ndd not in variants
ndd_notin_vars_lst = list(set(ndd_pr_lst).difference(variants_lst))
## ndd prs (not) in mobidblite_inidr: (ndds that are in IDR = 189 , not in IDR = 1119)
# (to be used to make new ndd_subdf based on these proteins only)
mobidblite_inidr_pr_lst = mobidblite_inidr_df['acc'].unique().tolist()
ndds_inidr_mobilite_lst = list(set(ndd_pr_lst).intersection(mobidblite_inidr_pr_lst))  # n: 189
# ndd_notin_mobilite_idr = list(set(ndd_pr_lst).difference(mobidblite_inidr_pr_lst))  # n: 1119
## ndd prs (not) in dismajority_inidr:
dismajority_inidr_pr_lst = dismajority_inidr_df['acc'].unique().tolist()
ndds_inidr_dismaj_lst = list(set(ndd_pr_lst).intersection(dismajority_inidr_pr_lst))  # n: 284
# ndd_notin_dismajority_idr = list(set(ndd_pr_lst).difference(dismajority_inidr_pr_lst))  # n: 1024
## Brain
brain_prot_lst = brn.brain_pr_lst_generator()  # n: 8320




