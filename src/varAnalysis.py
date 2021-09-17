# this is to perform previous analysis which was done in mobibool and protanalyse but for 'variants in IDR' dataset
# Sep 17th
import pandas as pd
import config as cfg
# TODO: insert each of this in mobibool instead of mobidb, but you need to also delete ndd and brain that are not in
#  this list maybe use left merge!

mobidblite_inidr_df = pd.read_csv(cfg.data['vars'] + '/muts-in-IDR-mobidb_lite.csv')
dismajority_inidr_df = pd.read_csv(cfg.data['vars'] + '/muts-in-IDR-disorder_majority.csv')
