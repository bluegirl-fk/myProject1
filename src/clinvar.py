from io import StringIO
import pandas as pd
import config as cfg

# ndd_clinvar = pd.read_csv(cfg.data['clin'] + '/ndd-clin-pos-merged.tsv') annovar = pd.read_csv(cfg.data['clin'] +
# '/query.output.exome_summary.csv') annovar['gnomAD_exome_ALL'] = annovar.loc[(annovar['gnomAD_exome_ALL'] != '.') &
# (annovar['gnomAD_exome_ALL'] != 'nan'), 'gnomAD_exome_ALL'].astype(float).apply(lambda x: '%.15f' % x)
# annovar.to_csv(cfg.data['clin'] + '/annovar-ndd-vars.tsv')

## from here started at Monday Aug 8th
def expand_regions(region_ranges_lst):
    transformed_regions = []
    region_ranges_lst = region_ranges_lst.split(',')
    # region_ranges_lst = list(region_ranges_lst)
    for region in region_ranges_lst:
        start = int(region.split('..')[0])
        end = int(region.split('..')[1])
        while start <= end:
            transformed_regions.append(start)
            start += 1
    return set(transformed_regions)


def mutidr_bool_array_maker(input_df, position_col_name):
    # input df is merged df of mobidb and mutation positions from mutinfo
    ## checks if mutation position is in startend disorder region of mobidb or not
    array_is_in = []  # will be filled with boolean of 0,1
    for index, row in input_df.iterrows():
        set_disorder_region = expand_regions(row.startend)  # temp set of data, convert each startend lst to a set,
        # write in report
        if row.position_col_name in set_disorder_region:
            array_is_in.append('1')
        else:
            array_is_in.append('0')
    input_df['isin_idr'] = array_is_in
    return array_is_in


## clinvar version: variant summary.txt 135,193,289	2022-08-01 17:10:07	2022-08-01 17:10:07
clinvar = pd.read_csv(cfg.data['clin'] + '/variant_summary.txt', sep='\t', low_memory=False)
## Gene4denovo candidate genes   Release version 1.2 (Updated: 07/08/2022)
candidates = pd.read_excel(cfg.data['gene4'] + '/Candidate_gene_1.2.xlsx', engine='openpyxl', usecols=['Groups', 'Gene_symbol', 'FDR'])
candidates = candidates.loc[candidates['FDR'] < 0.05]
cand_genes = candidates['Gene_symbol'].unique().tolist()
print(', '.join(cand_genes))
## only clinvar row that are representing our candidate genes
# clinvar_ndd = clinvar.loc[clinvar.GeneSymbol.isin(cand_genes)]
# clinvar_ndd.to_csv(cfg.data['clin'] + '/clinvar-ndd-candidate-genes.csv')
# clinvar_ndd = pd.read_csv(cfg.data['clin'] + '/clinvar-ndd-candidate-genes.csv')
# clinvar_ndd = clinvar_ndd.loc[clinvar_ndd['Type'] == 'single nucleotide variant']
# clinvar_ndd[['Name', 'pr_change']] = clinvar_ndd['Name'].str.split('p.', 1, expand=True)
# clinvar_ndd = clinvar_ndd.dropna(subset=['pr_change'])
# clinvar_ndd['pr_change'] = clinvar_ndd['pr_change'].str.replace(r')', '')
# clinvar_ndd = clinvar_ndd.loc[~clinvar_ndd.pr_change.str.contains('=')]
# clinvar_ndd['ref_aa'] = clinvar_ndd['pr_change'].str[:3]
# clinvar_ndd['alt_aa'] = clinvar_ndd['pr_change'].str[-3:]
# clinvar_ndd['pr_pos'] = clinvar_ndd.pr_change.str.extract('(\d+)')
# ## Mobidb downloaded on Aug 8th 2022, 78,106 entries
# mobidb = pd.read_csv(cfg.data['clin'] + '/mobidb_result_2022-08-08T13_12_39.379Z.tsv', sep='\t')
# mobidb = mobidb.rename(columns={'start..end': 'startend'})
#
# mutidr_bool_array_maker()