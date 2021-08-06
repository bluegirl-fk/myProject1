# this file is to read the denovo-db.csv and then see what happens next! Aug 5th
## version=denovo-db.v.1.6.1

import pandas as pd
import config as cfg

dndb = pd.read_csv(cfg.data['dndb'] + '/denovo-db.ssc-samples.variants.v.1.6.1.tsv',
                   skiprows=[0], sep='\t', low_memory=False)
dndb_rsID_True = dndb.loc[dndb['rsID'] != 0]  # (27169, 31)
dndb_valid = dndb.loc[dndb['Validation'] == 'yes']  # (9419,31)
dndb_valid_rsID = dndb_rsID_True.loc[dndb_rsID_True['Validation'] == 'yes']  # (462, 31)

rsID_lst = dndb_rsID_True['rsID'].unique().tolist()
refseq_lst = dndb['Transcript'].unique().tolist()

with open(cfg.data['dndb']+'/denovodb-rsid.txt', 'w') as f:
    for item in rsID_lst:
        f.write("%s\n" % item)

with open(cfg.data['dndb']+'/denovodb-refseq.txt', 'w') as f:
    for item in refseq_lst:
        f.write("%s\n" % item)

# also mapped to 14000 proteins!!! (put refseq into uniprot anmd idmap)