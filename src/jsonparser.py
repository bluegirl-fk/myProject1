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
    if item.get('ptms') is not None:
        # convert to a flat list for easier var analysis
        ptm_pos = [item for sublist in item.get('ptms').get('positions') for item in sublist]
        # convert list to 1 str, these 2 lines could be replaced by item.get('ptms').get('position') to reverse changes
        ptm_pos = ','.join(str(v) for v in ptm_pos)
        acc_ptm_info_lst.append([item.get('acc'), item.get('ptms').get('names'), ptm_pos])
    elif item.get('ptms') is None:
        acc_ptm_info_lst.append([item.get('acc'), None, None])

ptms_df = pd.DataFrame(acc_ptm_info_lst, columns=['acc', 'ptm_type', 'ptm_pos'])
ptms_df.to_csv(cfg.data['ptm'] + '/mobidb-parsed-ptm-info.csv')
