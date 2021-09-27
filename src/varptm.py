import pandas as pd
import config as cfg


ptms_df = pd.read_csv(cfg.data['ptm'] + '/mobidb-parsed-ptm-info.csv', usecols=['acc', 'ptm_pos'])


disorder_maj = pd.read_csv(cfg.data['vars'] + '/disorder-majority-inout-idr-vars-count-normalized.csv')
dismaj_ptm_df = pd.merge(disorder_maj, ptms_df, on='acc')
del dismaj_ptm_df['Unnamed: 0']

for row in dismaj_ptm_df.iterrows():
    for i in :
        if int(row['position']) == i:
            row['ptm_var'] = 1
        else:
            row['ptm_var'] = 0

