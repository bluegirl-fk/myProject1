#!/usr/bin/env python3

import config as cfg
import pandas as pd
import sys
import limitedmutations as lmut


def main():
    # generate_mutation_file()
    generate_mutation_file2()


def prep_orig_df(input_df):  # this gets df, renames some columns and deletes extra columns
    del input_df['Unnamed: 0']  # TODO:delete the rest of columns here instead of doing it in the end, easier processing
    input_df = input_df.reset_index()
    input_df = input_df.rename(columns={'index': 'idx1', 'AAChange.refGene': 'AAChange_refGene'})
    return input_df


def pr_mut_subdf_handler(input_org_df, desired_col, new_cols_count, new_cols_lst):  # gets one column and turns it into
    # new df with several columns, aim is to get mutation position and other info in separate columns
    # stacked the refseq mut positions, now have repeated proteins with possible dif mut positions per each
    subdf1 = input_org_df[[desired_col]]
    subdf1 = subdf1[desired_col].str.split(',', expand=True).stack().to_frame(desired_col)
    subdf2 = subdf1[desired_col].str.split(pat=':', expand=True)
    subdf3 = subdf2[subdf2.columns[:new_cols_count]]  # del rest of None columns
    subdf3.columns = new_cols_lst

    subdf3['mutPr'] = subdf3.AAChange_refGene.str.split(pat='fs*', expand=True, )
    subdf3['aa1'] = subdf3['mutPr'].str[2]
    subdf3['aa2'] = subdf3['mutPr'].str[-1]
    subdf3['position'] = subdf3['mutPr'].str.replace(r'\D', '')
    subdf3['frameshift'] = subdf3['AAChange_refGene'].str.split('fs', 1).str[1]
    del subdf3['mutPr']
    subdf3.to_csv(cfg.data['gene4'] + '/subdf-mut-beforeACC.csv')
    # TODO: consider reset_index here, fist check if you are not using it anywhere
    return subdf3


def merge_dfs_on_index(df1, df2, mutual_col, saving_route):
    # changes col names, merge based on index, cuz mutinfo was extracted from g4dn col originally
    # and should merge back correctly
    df1 = pd.read_csv(cfg.data['gene4'] + '/subdf-mut-beforeACC.csv')  # (201372, 11)
    df1 = df1.rename(columns={'Unnamed: 0': 'idx1', 'Unnamed: 1': 'sub-idx1'})
    merged_df = pd.merge(df1, df2, on=mutual_col)
    merged_df.to_csv(cfg.data['gene4'] + saving_route)
    return merged_df


def refseq_acc_df_handler(input_df, file_name):
    input_df.columns = ['refseq_id', 'isoforms', 'acc', 'organism', 'Length', 'Gene names']
    del input_df['organism']  # (50930, 5)
    df2 = input_df['refseq_id'].str.split(',', expand=True).stack()
    idx_tmp = df2.index._get_level_values(0)
    df3 = input_df.iloc[idx_tmp].copy()
    df3['refSeq'] = df2.values
    del df3['refseq_id']  # (93949, 5)
    df3.to_csv(cfg.data['gene4'] + file_name)
    return df3


def g4dn_mut_acc_merger(df1, df2, merger_column, file_name):
    mrg_df = pd.merge(df1, df2, on=merger_column)
    del mrg_df['Gene names']
    del mrg_df['idx1']
    del mrg_df['sub-idx1']
    mrg_df['position'] = mrg_df['position'].fillna(0).astype(int)  # (551773, 169)
    mrg_df = mrg_df.drop_duplicates(ignore_index=True)
    mrg_df.reset_index(level=0, inplace=True)  # (236699, 169)
    mrg_df.to_csv(cfg.data['gene4'] + file_name)
    return mrg_df


def expand_regions(region_ranges_lst):
    transformed_regions = []
    for reg in region_ranges_lst:
        start = int(reg.split('..')[0])
        end = int(reg.split('..')[1])
        while start <= end:
            transformed_regions.append(start)
            start += 1
    return set(transformed_regions)


def mutidr_bool_array_maker(input_df):
    # input df is merged df of mobidb and mutation positions from mutinfo
    ## checks if mutation position is in startend disorder region of mobidb or not
    array_is_in = []  # will be filled with boolean of 0,1
    for index, row in input_df.iterrows():
        set_disorder_region = expand_regions(row.startend)  # temp set of data, convert each startend lst to a set,
        # write in report
        if row.position in set_disorder_region:
            array_is_in.append('1')
        else:
            array_is_in.append('0')
    return array_is_in


# made a column with bool values to check if it is ndd protein or not
def mobi_mut_inidr_checker(mobi_df, mutinfo_df, filename):
    # mutinfo_df is mut_acc_mrg_df and is needed to get mut position from to check if in idr based on mobidb startend
    # converting disorder content ranges in each cell to a list
    mobi_df['startend'] = mobi_df['startend'].str.split(',')
    # subdf of mut pos to be merged with mobidb
    mutinfo_subdf = mutinfo_df[['index', 'acc', 'position', 'Phenotype', 'Length']]
    mutinfo_subdf['ndd'] = 1

    ##  lots of rows cuz accs are repeated in both databases with dif features or mutation per each ACC
    # (maybe it could be done with pivot table already and have a matrix instead, like before with the heatmaps)
    mobi_mutpos_df = pd.merge(mobi_df, mutinfo_subdf, on='acc')
    mobi_mutpos_df['ndd'] = mobi_mutpos_df['ndd'].fillna(0)

    array_is_in = mutidr_bool_array_maker(mobi_mutpos_df)
    ## add bool array to the df
    mobi_mutpos_df['is_in_startend'] = array_is_in
    mobi_mutpos_df.to_csv(cfg.data['gene4'] + filename)  # (4258689, 10)
    return mobi_mutpos_df


def mobi_mut_in_df_generator(input_df):
    # for mutations inside IDR
    # input_df = mobi_mutpos_checked_df # output = 6 dataframes
    # this will generate different possible dataframes based on:
    # cf_or_cc: if we have the pivot table, we can have either content_fraction or content_count as df values
    # in_or_out: if mutation is in idr or not
    mobi_mut_in_df = input_df[input_df['is_in_startend'] == 1]  # mutations inside IDR# ()
    mobi_mut_in_df.to_csv(cfg.data['gene4'] + '/mobidb-mut-pos-true.csv')
    # if needed not-pivoted(original mobidb to be merged), here merge mut_acc_mrg_df with mobi_mut_in_df
    mut_in_cf_df = mobi_mut_in_df.pivot_table(index=['acc'], columns=['feature'],
                                              values='content_fraction').fillna(0)
    mut_in_cf_df.to_csv(cfg.data['gene4'] + '/mobidb-pivot-cf-mut-true.csv')
    mut_in_cc_df = mobi_mut_in_df.pivot_table(index=['acc'], columns=['feature'],
                                              values='content_count').fillna(0)
    mut_in_cc_df.to_csv(cfg.data['gene4'] + '/mobidb-pivot-cc-mut-true.csv')
    return mobi_mut_in_df, mut_in_cf_df, mut_in_cc_df


def mobi_mut_out_df_generator(input_df):
    # for mutations not inside IDR
    mobi_mut_out_df = input_df[input_df['is_in_startend'] == 0]  # ()
    mobi_mut_out_df.to_csv(cfg.data['gene4'] + '/mobidb-mut-pos-false.csv')
    # if needed not-pivoted(original mobidb to be merged), here merge mut_acc_mrg_df with mobi_mut_out_df
    mut_out_cf_df = mobi_mut_out_df.pivot_table(index=['acc'], columns=['feature'],
                                                values='content_fraction').fillna(0)
    mut_out_cf_df.to_csv(cfg.data['gene4'] + '/mobidb-pivot-cf-mut-false.csv')
    mut_out_cc_df = mobi_mut_out_df.pivot_table(index=['acc'], columns=['feature'],
                                                values='content_count').fillna(0)
    mut_out_cc_df.to_csv(cfg.data['gene4'] + '/mobidb-pivot-cc-mut-flase.csv')
    return mobi_mut_out_df, mut_out_cf_df, mut_out_cc_df


# TODO: this should merge based on index, not acc!
def mobi_g4dn_merger(df1, df2, df3, df4, df5, df6, g4dn_mutinfo_acc_df):
    # inputs: dif categorized mobidb dfs + mut_acc_mrg_df
    df1_g4dn_mrgdf = pd.merge(df1, g4dn_mutinfo_acc_df, on='acc')
    df2_g4dn_mrgdf = pd.merge(df2, g4dn_mutinfo_acc_df, on='acc')
    df3_g4dn_mrgdf = pd.merge(df3, g4dn_mutinfo_acc_df, on='acc')
    df4_g4dn_mrgdf = pd.merge(df4, g4dn_mutinfo_acc_df, on='acc')
    df5_g4dn_mrgdf = pd.merge(df5, g4dn_mutinfo_acc_df, on='acc')
    df6_g4dn_mrgdf = pd.merge(df6, g4dn_mutinfo_acc_df, on='acc')
    # write to_csv
    df1_g4dn_mrgdf.to_csv(cfg.data['gene4'] + '/' + df1 + '-mrg.csv')
    df2_g4dn_mrgdf.to_csv(cfg.data['gene4'] + '/' + df2 + '-mrg.csv')
    df3_g4dn_mrgdf.to_csv(cfg.data['gene4'] + '/' + df3 + '-mrg.csv')
    df4_g4dn_mrgdf.to_csv(cfg.data['gene4'] + '/' + df4 + '-mrg.csv')
    df5_g4dn_mrgdf.to_csv(cfg.data['gene4'] + '/' + df5 + '-mrg.csv')
    df6_g4dn_mrgdf.to_csv(cfg.data['gene4'] + '/' + df6 + '-mrg.csv')

    return df1_g4dn_mrgdf, df2_g4dn_mrgdf, df3_g4dn_mrgdf, df4_g4dn_mrgdf, df5_g4dn_mrgdf, df6_g4dn_mrgdf


def generate_mutation_file():
    # ### G4dn code more ## Gene4denovo ## only exonic mutations # g4dn_exonic_df = pd.read_csv(
    # 'data/gene4denovo/exonic-df.csv')  # (70879, 156) # stacked the refseq mut positions, now have repeated
    # ...............

    df = pd.read_csv(cfg.data['gene4'] + '/yourfile.csv', sep='\t')
    # ...............


def generate_mutation_file2():
    # ### G4dn code more ## Gene4denovo ## only exonic mutations # g4dn_exonic_df = pd.read_csv(
    # 'data/gene4denovo/exonic-df.csv')  # (70879, 156) # stacked the refseq mut positions, now have repeated
    # ...............

    df = pd.read_csv(cfg.data['gene4'] + '/yourfile.csv', sep='\t')
    # ...............

    # if __name__ == '__main__':
    g4dn_exonic_df = pd.read_csv(cfg.data['gene4'] + '/exonic-df.csv')
    g4dn_exonic_df = prep_orig_df(g4dn_exonic_df)  # (70879, 156)
    refseq_mut_subdf = pr_mut_subdf_handler(g4dn_exonic_df, 'AAChange_refGene', 10,
                                            ['Gene_refGene', 'refSeq', 'exon#', 'mutNA', 'AAChange_refGene', 'aa1',
                                             'aa2', 'position', 'frameshift', 'mutPr'])  # (201372, 9)
    # (201372, 166)
    g4dn_exo_mutinfo_df = merge_dfs_on_index(refseq_mut_subdf, g4dn_exonic_df, 'idx1', '/exonic-mutinfo.csv')
    ## positive candidate genes list
    pos_candidate_gene_lst, ctrl_candidate_genes_lst = lmut.genes_lst_maker()  # 181
    # (24387, 166)
    g4dn_exo_pos_cand_df = g4dn_exo_mutinfo_df[g4dn_exo_mutinfo_df.Gene_refGene.isin(pos_candidate_gene_lst)]

    # * Got refseq_ids from refseq_mut_subdf['refSeq'] and wrote this list to txt, retrieved ACCs from uniprot
    # (split my text file using bash : split -l 70000 refseq-gene4dn.txt, the 7000 is number of the lines)

    ## g4dn mutInfo + uniprot ACCs file (merge)
    refseq_acc_df = pd.read_csv(cfg.data['rseq'] + '/refseq-acc.tab', sep='\t')  # from Uniprot
    refseq_acc_modified_df = refseq_acc_df_handler(refseq_acc_df, '/refseg-acc-modified.csv')  # (93949, 5)

    ## merge g4dn exonic mutInfo with Uniprot ACC # (22858, 168)
    mut_acc_mrg_df = g4dn_mut_acc_merger(refseq_acc_modified_df, g4dn_exo_pos_cand_df, 'refSeq',
                                         '/mut-acc-mrg-df100.csv')
    # TODO works fine till here

    ## mobidb
    ## delete this line later


if __name__ == '__main__':
    mut_acc_mrg_df = pd.read_csv(cfg.data['gene4'] + '/mut-acc-mrg-df100.csv', low_memory=False)

    mobidb_original_df = pd.read_csv(cfg.data[''] + '/mobidb_result.tsv', sep='\t')  # (1212280,6)
    mobidb_original_df.columns = ['acc', 'feature', 'startend', 'content_fraction', 'content_count', 'length']
    # (1212282, 11)
    mobi_mutpos_checked_df = mobi_mut_inidr_checker(mobidb_original_df, mut_acc_mrg_df, '/mut-pos-mobi100.csv')
    # the problem with using join instead of concat is that it does not contain all mobidb rows
    sys.exit()

# maybe make another file to process that one

# merged_mobidbp_g4dn_df = pd.read_csv(
# 'data/gene4denovo/final-merged-with-mobidb-pivot.csv', low_memory=False)
#       (180315, 245) phenotypes_lst = [
# 'ASD', 'EE', 'ID', 'CMS', 'SCZ', 'NDDs'] # (41081, 245) phens_mobip_g4dn_muttrue_df = merged_mobidbp_g4dn_df[
# merged_mobidbp_g4dn_df.Phenotype.isin(phenotypes_lst)] phens_mobip_g4dn_limited_df =
# phens_mobip_g4dn_muttrue_df.drop( columns=['Unnamed: 0', 'index', 'Unnamed: 0.1', 'mutNA',
# 'AAChange_refGene_x', 'Rare_or_Common', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene',
# 'AAChange_refGene_y', 'GeneFullName.refGene', 'GeneFullName.ensGene', 'GeneFunction.ensGene',
# 'GeneExpressionTissue.ensGene', 'GeneDisease.ensGene', 'OMIM.ensGene', 'MGI.ensGene', 'RVIS.ensGene',
# 'LoFtool.ensGene', 'GDI.ensGene', 'Episcore.ensGene', 'Aggarwala.ensGene', 'pLi_EXAC.ensGene',
# 'HIPred.ensGene'])  # (41081, 221) phen_asd = ['ASD'] ASD_mobip_g4dn_limited_df = phens_mobip_g4dn_limited_df[
# phens_mobip_g4dn_limited_df.Phenotype.isin(phen_asd)]  # (28486, 221)
#
# g4dn_llps_subdf = phens_mobip_g4dn_limited_df.loc[phens_mobip_g4dn_limited_df['curated-phase_separation-merge'] !=
#                                                   0.0, ('acc', 'curated-phase_separation-merge')]
# g4dn_llps_subdf = g4dn_llps_subdf.drop_duplicates()
#
# g4dn_llps_subdf.to_csv(r'data/gene4denovo-llps-mrg-Pr-list.csv')
#
# ## lst of column names
# mrg_mobip_d4dn_cols_lst = list(merged_mobidbp_g4dn_df.columns)
#
# # percentage of idr mutations (merged_mobidbp_g4dn_df/all(mut_acc_mrg_df) )=> 180315/236699 = 76.17 % in idrs #
# print(ASD_mobip_g4dn_limited_df['ExonicFunc.refGene'].value_counts().keys()[0:1])  # => ['nonsynonymous SNV'] #
# print(ASD_mobip_g4dn_limited_df['Chr'].value_counts().keys()[0:5])  # => ['2', '1', '3', '19', '12'] # print(
# ASD_mobip_g4dn_limited_df['GeneFunction.refGene'].value_counts().keys()[0:4])   # => # [1- May be involved in
# transcriptional regulation., 2- Key component in the assembly and functioning of vertebrate # striated muscles.
# By providing connections at the level of individual microfilaments, it contributes to the fine # balance of
# forces between the two halves of the sarcomere. The size and extensibility of the cross-links are the # main
# determinants of sarcomere extensibility properties of muscle. In non-muscle cells, seems to play a role in #
# chromosome condensation and chromosome segregation during mitosis. Might link the lamina network to chromatin
# or # nuclear actin, or both during interphase. {ECO:0000269|PubMed:9804419}. # 3- Voltage-sensitive calcium
# channels (VSCC) mediate the entry of calcium ions into excitable cells and are also # involved in a variety of
# calcium-dependent processes, including muscle contraction, hormone or neurotransmitter # release,
# gene expression, cell motility, cell division and cell death. The isoform alpha-1C gives rise to L-type #
# calcium currents. Long-lasting (L-type) calcium channels belong to the ' high - voltage activated ' (HVA)
# group. # They are blocked by dihydropyridines (DHP), phenylalkylamines, benzothiazepines,
# and by omega-agatoxin-IIIA # (omega-Aga-IIIA). They are however insensitive to omega-conotoxin- GVIA (
# omega-CTx-GVIA) and omega-agatoxin-IVA # (omega-Aga-IVA). Calcium channels containing the alpha-1C subunit play
# an important role in excitation-contraction # coupling in the heart. The various isoforms display marked
# differences in the sensitivity to DHP compounds. # Binding of calmodulin or CABP1 at the same regulatory sites
# results in an opposit effects on the channel function. # {ECO:0000269|PubMed:12176756,
# ECO:0000269|PubMed:17071743, ECO:0000269|PubMed:7737988, # ECO:0000269|PubMed:8392192,
# ECO:0000269|PubMed:9013606, ECO:0000269|PubMed:9607315}., '], dtype = 'object') print(
# ASD_mobip_g4dn_limited_df['ExonicFunc.refGene'].value_counts(dropna=False))
#
# # filtered with mutations out of idrs
# # mobidb_mut_check_df = pd.read_csv('data/mutations-position-mobidb-all.csv')
# # filtered_mut_pos_false_df = mobidb_mut_check_df[mobidb_mut_check_df['is_in_startend'] == 0]
# # filtered_mut_pos_false_df.to_csv(r'data/gene4denovo/mobidb-mut-pos-false.csv')  # (3255439, 10)
# filt_mut_pos_false_df = pd.read_csv('data/gene4denovo/mobidb-mut-pos-false.csv')
# mobidb_mutfalse_cf_pivot_df = filt_mut_pos_false_df.pivot_table(index=['acc'], columns=['feature'],
#                                                                 values='content_fraction').fillna(0)
# mobidb_mutfalse_cc_pivot_df = filt_mut_pos_false_df.pivot_table(index=['acc'], columns=['feature'],
#                                                                 values='content_count').fillna(0)
# # TODO: also the content count for muttrue df
#
# ## merged mobidb_muttrue(pivot df) with (g4dn+acc)
# mobidbp_g4dn_cf_mutfalse_df = pd.merge(mobidb_mutfalse_cf_pivot_df, mut_acc_mrg_df, on='acc')
# mobidbp_g4dn_cf_mutfalse_df.to_csv(r'data/gene4denovo/merged-with-mobidb-pivot-mutfalse.csv')
# merged_mobidbp_g4dn_mutfalse_df = pd.read_csv('data/gene4denovo/merged-with-mobidb-pivot-mutfalse.csv',
#                                               low_memory=False)
# phenotypes_lst = ['ASD', 'EE', 'ID', 'CMS', 'SCZ', 'NDDs']
# phens_mobip_g4dn_mutfalse_df = merged_mobidbp_g4dn_mutfalse_df[
#     merged_mobidbp_g4dn_mutfalse_df.Phenotype.isin(phenotypes_lst)]
# phens_mobip_g4dn_mutfalse_df = phens_mobip_g4dn_mutfalse_df

sys.exit(main())
