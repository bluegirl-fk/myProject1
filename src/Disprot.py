# this is to check my hotspot mutations whitin idr in disprot to find the subregions and structural characteristics
import pandas as pd
import config as cfg
import varAnalysis as var

disprot = pd.read_csv(cfg.data['dis'] + '/DisProt release_2021_08 with_ambiguous_evidences.tsv', sep='\t')
variant_count = pd.read_csv(cfg.data['vars'] + '/variants-count-and-hotspots.csv')
variant_count = variant_count.loc[variant_count['percentage'] > 50]
variant_count = variant_count.rename(columns={'Unnamed: 0': 'acc'})

var_disprot = pd.merge(disprot, variant_count, on='acc')

ndd_var_lst_50 = var.ndd_var_lst()
ndd_disprot = disprot.loc[disprot.acc.isin(ndd_var_lst_50)]
