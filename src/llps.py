import pandas as pd
import config as cfg
import matplotlib.pyplot as plt
import seaborn as sns

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
    ## the two llps datasets have 28 proteins in common

    varin_phasep_df, ndd_varin_phasep_df = var_llps_df_merger('phasep')  # 2296 # 840
    varin_phasepro_df, ndd_varin_phasepro_df = var_llps_df_merger('phasepro')  # 134 # 30
    # lst of NDD proteins with var in IDR and involved in LLPS
    hs_varin_phasep_lst = varin_phasep_df['acc'].unique().tolist()  # 640
    hs_varin_phasepro_lst = varin_phasepro_df['acc'].unique().tolist()  # 23
    ndd_varin_phasep_lst = ndd_varin_phasep_df['acc'].unique().tolist()  # 41
    ndd_varin_phasepro_lst = ndd_varin_phasepro_df['acc'].unique().tolist()  # 5

    # print('\n'.join(ndd_varin_phasepro_lst))
    # print('\n'.join(ndd_varin_phasep_lst))
    # check if llps proteins have more disordered regions (maybe), more low complexity
    # proteins high-confidence associate to MLOs have more PTMs
    # see your disease proteins are mostly in which llps category, (driver/scaffold, other types ...)
    # this paper : https://www.sciencedirect.com/science/article/pii/S2001037021002804
    # see this protein opened in mobidb in your df and compare
    fpro_df, fpro_lst = phase_pro()  # 120 proteins
    fsep_df, fsep_lst = phasep()  # 2572
    # see where this proteins are involves, the ones related to ndds, for fpro they are 5 and not in a patway together,
    # but all important and since are not assigned to a disease in the fpro, this can gives ideas where to look.
    # the ndd proteins var in idr in llps based on fsep are 30-40 and in a network, so even better!
    # in the end you can merge them in a list or something, giveout as a database
    # or try the one which is merged, form the above like if they fix the server!

    mlo_disease = pd.read_csv(cfg.data['mlo-d'] + '/mlodisdb_components.csv')
    # from https://academic.oup.com/bib/article/22/4/bbaa271/5943794
    # mlo-dis database not so powerful, could be done better
    mlo_dis = pd.read_excel(cfg.data['mlo-d'] + '/mlodisdb_relations_all.xlsx', engine='openpyxl')

    ## new mlo dis based on my ndds, var in idr ndds, merged two llps dfs
    ndd_llps_merged_lst = ndd_varin_phasep_lst + ndd_varin_phasepro_lst
    print('\n'.join(ndd_llps_merged_lst))
    mlo_disease_ndd = mlo_disease.loc[mlo_disease.Entry.isin(ndd_llps_merged_lst)]

    ## file from this paper:https://www.sciencedirect.com/science/article/pii/S2001037021002804
    ## challenge, open this file correctly, for now not very necessary
    # mlo = pd.read_excel(cfg.data['mlo']+'/complementary-data-paper-S2001037021002804.xlsx', engine='openpyxl')
    ## STRING analysis of ndd_llps_merged_lst with 46 total proteins
    molecular_function = pd.read_csv(cfg.data['fs-str'] + '/enrichment.Function.tsv', sep='\t')
    cell_process = pd.read_csv(cfg.data['fs-str'] + '/enrichment.Process.tsv', sep='\t')
    # count number of variants in llps dataset
    disorder_majority = pd.read_csv(cfg.data['vars'] + '/disorder-majority-inout-idr-vars-count-normalized.csv', usecols=['acc', 'in_idr_vars_perc'])
    llps_ndd_dismaj = disorder_majority.loc[disorder_majority.acc.isin(ndd_llps_merged_lst)]
    llps_ndd_dismaj = llps_ndd_dismaj.drop_duplicates()
    llps_ndd_dismaj = llps_ndd_dismaj.rename(columns={"in_idr_vars_perc": "IDR variant fraction (NDD-associated Proteins with LLPS roles)"})

    ndd_subdf = pd.read_csv(cfg.data['phens-fdr'] + '/acc-phen-5percentFDR.csv')
    ndd_pr_lst = ndd_subdf['acc'].unique().tolist()  # 1308 proteins
    ndd_dismaj = disorder_majority.loc[disorder_majority.acc.isin(ndd_pr_lst)]
    ndd_dismaj = ndd_dismaj.drop_duplicates()
    ndd_dismaj = ndd_dismaj.rename(columns={"in_idr_vars_perc": "IDR variant fraction (All NDD-associated proteins)"})
    g = sns.violinplot(x=ndd_dismaj['IDR variant fraction (All NDD-associated proteins)'])
    g.set_xlim(0, 1)
    plt.show()


