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
valid_rsID_lst = dndb_valid_rsID['Transcript'].unique().tolist()
with open(cfg.data['dndb']+'/denovodb-valid-rseq.txt', 'w') as f:  # maps to 249 uniprot IDs
    for item in valid_rsID_lst:
        f.write("%s\n" % item)


valid_acc_df = pd.read_csv(cfg.data['dndb'] + '/valid-acc-from-rseq.tab')
valid_pr_lst = valid_acc_df['Entry'].unique().tolist()
phen_df = pd.read_csv(cfg.data['gene4'] + '/positive_cand_g4mobi_concat.csv', usecols=["acc_x", "Phenotype"])
phens_acc_lst = phen_df['acc_x'].unique().tolist()
l1 = [i for i in phens_acc_lst if i in valid_pr_lst]  # only 5 mutual proteins among my list and denovodb! maybe I could merge them!


