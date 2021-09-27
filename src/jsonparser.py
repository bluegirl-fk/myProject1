# json file from mobidb to get ptm positions
# September 27th of 2021
import json
import pandas as pd
import config as cfg

acc_ptm_info_lst = []
with open(cfg.data['jsn'] + '/mobidb_result_2021-09-27T09_48_24.861Z.json', 'r') as f:
    json_data = f.read()
doc = json.loads(json_data)
for item in doc:
    acc = item.get('acc')
    ptm_type = item.get('ptms').get('names')
    ptm_pos = item.get('ptms').get('positions')
    acc_ptm_info_lst.append([acc, ptm_type, ptm_pos])

ptms_df = pd.DataFrame(acc_ptm_info_lst, columns=['acc', 'ptm_type', 'ptm_pos'])
ptms_df.to_csv(cfg.data['ptm'] + '/mobidb-parsed-ptm-info.csv')
