#!/usr/bin/env python3

import config as cfg
import pandas as pd
import sys


def main():
    # generate_mutation_file()
    generate_mutation_file2()


def prepare_df(input_df):
    del input_df['Unnamed: 0']
    input_df = input_df.reset_index()
    input_df = input_df.rename(columns={'index': 'idx', 'AAChange.refGene': 'AAChange_refGene'})
    return input_df


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


if __name__ == '__main__':
    g4dn_exonic_df = pd.read_csv('data/gene4denovo/exonic-df.csv')  # (70879, 156)
    g4dn_exonic_df = prepare_df(g4dn_exonic_df)

    # stacked the refseq mut positions, now have repeated
    # proteins but with possible dif mut positions per each protein # (several possible rows per protein) #
    # aachange_g4dn_subdf3.to_csv(r'data/gene4denovo/subdf-mut-beforeACC.csv')  # (201372, 9) # Then made a list of
    # refSeq ids from this df, wrote to .txt, retrived ACCs from uniprot # splited my text file using bash : split -l
    # 70000 refseq-gene4dn.txt, the 7000 is number of the lines
    #
    # ## mut position df import
    # # rseq_mutinfo_df = pd.read_csv('data/gene4denovo/subdf-mut-beforeACC.csv')  # (201372, 11)
    # # # merge with gene4dn exonic, now original exonic g4dn file + mutInfo e.g position
    # # g4dn_exonic_mutinfo_df = pd.merge(rseq_mutinfo_df, g4dn_exonic_df, on='idx')  # (201372, 166)
    # # g4dn_exonic_mutinfo_df.to_csv(r'data/gene4denovo/exonic-mutinfo.csv')
    #
    # # ## g4dn mutInfo + uniprot ACCs file
    # # g4dn_exonic_mutinfo_df = pd.read_csv('data/gene4denovo/exonic-mutinfo.csv')  # (201372, 166)
    # # refseq_acc_df1 = pd.read_csv('data/refseq/refseq-acc.tab', sep='\t')  # from Uniprot
    #
    # # # * merge g4dn exonic mutInfo with Uniprot ACC
    # mut_acc_mrg_df = pd.read_csv('data/mut-acc-mrg-df.csv')  # (236699, 169)
    #
    # ## mobidb
    # # mobidb_original_df = pd.read_csv('data/mobidb_result.tsv', sep='\t')  # (1212280,6)
    # # mobidb_original_df.columns = ['acc', 'feature', 'startend', 'content_fraction', 'content_count', 'length']
    # # # converting disorder content ranges in each cell to a list
    # # mobidb_original_df['startend'] = mobidb_original_df['startend'].str.split(',')
    # ## subdf of mut pos to be merged with mobidb
    # # mutinfo_subdf = mut_acc_mrg_df[['index', 'acc', 'position']]  # (236699, 3)
    #
    # ##  lots of rows cuz accs are repeated in both databases with dif features or mutation per each ACC,
    # # (this can be merged with d4dn based on idx)
    # # mobidb_mutpos_df = pd.merge(mobidb_original_df, mutinfo_subdf, on='acc')  # (4258689, 8)
    #
    # ## check if mutation position is in startend disorder region of mobidb or not
    # # array_is_in = []  # will be filled with boolean of 0,1 for pos in startend or not
    # # for index, row in mobidb_mutpos_df.iterrows():
    # #     set_disorder_region = expand_regions(row.startend)  # temp set of data, convert each startend lst to a set,
    # #     # write in report
    # #     if row.position in set_disorder_region:
    # #         print('yes')
    # #         array_is_in.append('1')
    # #     else:
    # #         print('no')
    # #         array_is_in.append('0')
    #
    # ## add bool array to the df
    # # mobidb_mutpos_df['is_in_startend'] = array_is_in
    # # mobidb_mutpos_df.to_csv(r'data/mutations-position-mobidb-all.csv')
    # # final_mut_check_df = pd.read_csv('data/mutations-position-mobidb-all.csv')  # (4258689, 10)
    #
    # # filtered with mutations inside IDRs
    # # filtered_mut_pos_df = final_mut_check_df[final_mut_check_df['is_in_startend'] == 1]  # (1003250, 10)
    # # filtered_mut_pos_df.to_csv(r'data/gene4denovo/mobidb-mut-pos-true.csv')
    # # mobidb_mutpos_true_df = pd.read_csv('data/gene4denovo/mobidb-mut-pos-true.csv')
    #
    # # mobidbp_muttrue_cf_df = mobidb_mutpos_true_df.pivot_table(
    # #     index=['acc'], columns=['feature'], values='content_fraction').fillna(0)  # also content_count could be used
    # ## merged mobidb_muttrue(normal df) with (g4dn+acc)
    # # # merged_filtered_mobidb_d4dn_df = pd.merge(filtered_mut_pos_df, mut_acc_mrg_df, on='index')
    # # # merged_filtered_mobidb_d4dn_df.to_csv(r'data/gene4denovo/final-merged-mobi-g4dn-true.csv')
    #
    # ## merged mobidb_muttrue(pivot df) with (g4dn+acc) # merged_mobidbp_g4dn_df = pd.merge(mobidbp_muttrue_cf_df,
    # mut_acc_mrg_df, on='acc') # merged_mobidbp_g4dn_df.to_csv(
    # r'data/gene4denovo/final-merged-with-mobidb-pivot.csv') merged_mobidbp_g4dn_df = pd.read_csv(
    # 'data/gene4denovo/final-merged-with-mobidb-pivot.csv', low_memory=False)  # (180315, 245) phenotypes_lst = [
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
    # # percentage of IDR mutations (merged_mobidbp_g4dn_df/all(mut_acc_mrg_df) )=> 180315/236699 = 76.17 % in IDRs #
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
    # # filtered with mutations out of IDRs
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