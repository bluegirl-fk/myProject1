# this class is made to handle new mobidb with the data for ndd proteins Damiano provided, these are not present in
# the original mobidb and I'm going to append it to the previous mobidb, also discard brain/ndd prs that are not
# in mobidb even with this applied changes  Sep 14th
import config as cfg
import pandas as pd
import re
import brain as brn

# TODO: check with Damiano what he meant by taking the numbers into account for stats although not in mobidb
# probably should use refrence list in the beginning, and explain your exclution and limitation criteria
# e.g for ndds there was 95% fdr, for human, brain and ndd only kept the ones existing in mobidb using the new mobidb
# (confirm from Damiano)
mobidb1 = pd.read_csv(cfg.data[''] + '/mobidb_result.tsv', sep='\t')
# new tsv that Damiano sent based on NDDs that were not in mobidb and now they are gonna be added!
mobidb2 = pd.read_csv(cfg.data['phens'] + '/ndd-not-in-mobidb-byDamiano.tsv', sep='\t')
mobidb = mobidb1.append(mobidb2, ignore_index=True)
mobidb.to_csv(cfg.data['phens'] + '/mobidb-results+ndd-tsv-damiano-shared.tsv')
mobidb_pr_lst = mobidb['acc'].unique().tolist()

reference_lst = []
with open(cfg.data['phens']+'/uniprot-organism__homo+sapiens_-filtered-proteome_UP000005640+AND+organi--.list', 'r') as f:
    for line in f:
        # add item to the list
        line = re.sub(r'\n', '', line)
        reference_lst.append(line)

# NDD
ndd_subdf = pd.read_csv(cfg.data['phens-fdr'] + '/acc-phen-5percentFDR.csv')
ndd_subdf = ndd_subdf.drop_duplicates()  # (4531, 3)
ndd_pr_lst = ndd_subdf['acc'].unique().tolist()  # 1308 proteins  => these are all, not the selected phens
## brain
brain_prot_lst = brn.brain_pr_lst_generator()  # n: 8320
## lists
ndd_not_in_mobidb = list(set(ndd_pr_lst).difference(mobidb_pr_lst))  # n = 11
brain_not_in_mobidb = list(set(brain_prot_lst).difference(mobidb_pr_lst))  # n = 26
##
reference_lst = reference_lst + ndd_pr_lst + brain_prot_lst
reference_lst = list(set(reference_lst))
l50 = list(set(mobidb_pr_lst).difference(reference_lst))  # mobidb has 573 proteins which is not even in reflist
# (the one damiano sent plus ndd plus brain!)
reference_lst = reference_lst + list(set(mobidb_pr_lst).difference(reference_lst))



