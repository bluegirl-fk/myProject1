import pandas as pd
import config as cfg

ndd_clinvar = pd.read_csv(cfg.data['clin'] + '/ndd-clin-pos-merged.tsv')