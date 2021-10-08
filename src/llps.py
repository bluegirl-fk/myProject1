import pandas as pd
import config as cfg

disorder_vars = pd.read_csv(cfg.data['vars'] + '/disorder-majority-inout-idr-vars-count-normalized.csv', usecols=
['acc', 'var_id', 'orig_aa', 'var_aa', 'position', 'isin_idr', 'total_vars', 'in_idr_vars', 'out_idr_vars'])
phasepro = pd.read_csv(cfg.data['llps'] + '/phasepro.tsv', sep='\t')
# phasepro.columns = ['Gene', 'Name', 'acc', 'Organism', 'Sequence', 'Gene_name', 'subgroup', 'index',
#                     'Membrane cluster', 'Partner dependent', 'RNA dependancy', 'PTM dependancy',
#                     'Domain-motif interaction', 'Discrete oligomerization']
