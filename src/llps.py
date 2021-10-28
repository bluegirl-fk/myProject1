import pandas as pd
import config as cfg

disorder_vars = pd.read_csv(cfg.data['vars'] + '/disorder-majority-inout-idr-vars-count-normalized.csv', usecols=
['acc', 'var_id', 'orig_aa', 'var_aa', 'position', 'isin_idr', 'total_vars', 'in_idr_vars', 'out_idr_vars'])
disorder_vars = disorder_vars.loc[disorder_vars['isin_idr'] == 1]
phasepro = pd.read_csv(cfg.data['fp'] + '/phasepro.tsv', sep='\t')
phasepro.columns = ['common_name', 'name', 'acc', 'organism', 'sequence', 'gene', 'taxon', 'id', 'segment',
                    'boundaries', 'region', 'partners', 'determinants', 'forms', 'organelles', 'pmids', 'description',
                    'experiment_llps', 'in_vivo', 'in_vitro', 'experiment_state', 'rna_req', 'ptm_affect', 'disease',
                    'splice', 'interaction', 'membrane_cluster', 'partner_dep', 'rna_dep', 'ptm_dep',
                    'domain-motif_interaction', 'discrete_oligo', 'under_annote', 'annotator', 'functional_class',
                    'date']
disorder_vars_lst = disorder_vars['acc'].unique().tolist()
phasepro_lst = phasepro['acc'].unique().tolist()
# todo: organzie this, for phsepro phasep, retrieve both "all variants" and omly "ndd variants", decide on var_in_idr
a = pd.merge(disorder_vars, phasepro, on='acc')
vars_phase_pro_lst = a['acc'].unique().tolist()
ndd_subdf = pd.read_csv(cfg.data['phens-fdr'] + '/acc-phen-5percentFDR.csv')
b = pd.merge(ndd_subdf, phasepro, on='acc')
b_lst = b['acc'].unique().tolist()
# print(','.join(b_lst))

phasep = pd.read_excel(cfg.data['fs'] + '/High throughput Data V1.3.xlsx', engine='openpyxl')
phasep_ndd = pd.merge(phasep, ndd_subdf, left_on='UniprotEntry', right_on='acc')
phasep_ndd_prs = phasep_ndd['acc'].unique().tolist()
print(','.join(phasep_ndd_prs))