import pandas as pd
import config as cfg


def var_llps_df_merger(source):
    # df with checked vars for being in idr or not (based on disorder majority consensus)
    disorder_majority = pd.read_csv(cfg.data['vars'] + '/disorder-majority-inout-idr-vars-count-normalized.csv',
                                    usecols=
                                    ['acc', 'var_id', 'orig_aa', 'var_aa', 'position', 'isin_idr', 'total_vars',
                                     'in_idr_vars', 'out_idr_vars'])
    in_idr_var_df = disorder_majority.loc[disorder_majority['isin_idr'] == 1]
    in_idr_var_lst = in_idr_var_df['acc'].unique().tolist()
    # NDD (in IDR vaqriants df but just for NDD protein)
    ndd_invar_df = ndd_var_df()
    # LLPS
    phasepro_df, phasepro_lst = phase_pro()
    phasep_df, phasep_lst = phasep()
    # merging
    if source == 'phasepro':
        var_phasepro_mrg = pd.merge(in_idr_var_df, phasepro_df, on='acc')
        ndd_var_phasepro_mrg = pd.merge(ndd_invar_df, phasepro_df, on='acc')
        return var_phasepro_mrg, ndd_var_phasepro_mrg
    elif source == 'phasep':
        var_phasep_mrg = pd.merge(in_idr_var_df, phasep_df, left_on='acc', right_on='UniprotEntry')
        ndd_var_phasep_mrg = pd.merge(ndd_invar_df, phasep_df, left_on='acc', right_on='UniprotEntry')
        return var_phasep_mrg, ndd_var_phasep_mrg


def phase_pro():
    phasepro = pd.read_csv(cfg.data['fp'] + '/phasepro.tsv', sep='\t')
    phasepro.columns = ['common_name', 'name', 'acc', 'organism', 'sequence', 'gene', 'taxon', 'id', 'segment',
                        'boundaries', 'region', 'partners', 'determinants', 'forms', 'organelles', 'pmids',
                        'description',
                        'experiment_llps', 'in_vivo', 'in_vitro', 'experiment_state', 'rna_req', 'ptm_affect',
                        'disease',
                        'splice', 'interaction', 'membrane_cluster', 'partner_dep', 'rna_dep', 'ptm_dep',
                        'domain-motif_interaction', 'discrete_oligo', 'under_annote', 'annotator', 'functional_class',
                        'date']
    phasepro_lst = phasepro['acc'].unique().tolist()
    return phasepro, phasepro_lst


def phasep():
    phasep = pd.read_excel(cfg.data['fs'] + '/High throughput Data V1.3.xlsx', engine='openpyxl')
    phasep = phasep.drop(['No', 'protein material states', 'Mutation/disease', 'Organism'], axis=1)
    phasep_lst = phasep['UniprotEntry'].unique().tolist()
    return phasep, phasep_lst


def ndd_var_df():
    disorder_majority = pd.read_csv(cfg.data['vars'] + '/disorder-majority-inout-idr-vars-count-normalized.csv',
                                    usecols=
                                    ['acc', 'var_id', 'orig_aa', 'var_aa', 'position', 'isin_idr', 'total_vars',
                                     'in_idr_vars', 'out_idr_vars'])
    in_idr_var_df = disorder_majority.loc[disorder_majority['isin_idr'] == 1]
    ndd_subdf = pd.read_csv(cfg.data['phens-fdr'] + '/acc-phen-5percentFDR.csv')
    del ndd_subdf['Unnamed: 0']
    ndd_inidr_var_mrg = pd.merge(in_idr_var_df, ndd_subdf, on='acc')
    return ndd_inidr_var_mrg


if __name__ == '__main__':
    varin_phasep_df, ndd_varin_phasep_df = var_llps_df_merger('phasep')  # 2296 # 840
    varin_phasepro_df, ndd_varin_phasepro_df = var_llps_df_merger('phasepro')  # 134 # 30
    # lst of NDD proteins with var in IDR and involved in LLPS
    ndd_varin_phasep_lst = ndd_varin_phasep_df['acc'].unique().tolist()  # 41
    # ndd_varin_phasepro_lst = ndd_varin_phasepro_df['acc'].unique().tolist()  # 5
    print('\n'.join(ndd_varin_phasep_lst))