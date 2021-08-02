# this file is supposed to import the mutations inside and outside IDRs separately and perform different operations:
# a df with mutations filtered based on mobidb-lite feature
# disorder residue length + complement
# July 27th
import pandas as pd
import config as cfg
import sys
import limitedmutations as lmut


def main():
    print("whatever")
    return


if __name__ == '__main__':
    # mobi_mut_in_idr_df = pd.read_csv(cfg.data['gene4'] + '/mobidb-mut-pos-true.csv')  # (994383, 11)
    # mobi_mut_out_idr_df = pd.read_csv(cfg.data['gene4'] + '/mobidb-mut-pos-false.csv')  # (3238392, 11)
    #
    # # (25961, 11)
    # mut_in_mobi_lite_df = mobi_mut_in_idr_df.loc[mobi_mut_in_idr_df['feature'] == 'prediction-disorder-mobidb_lite']
    # # (102495, 11)
    # mut_out_mobi_lite_df = mobi_mut_out_idr_df.loc[mobi_mut_out_idr_df['feature']=='prediction-disorder-mobidb_lite']

    ## merge with g4dn (proteins mutation position checked)

    in_mobi_g4dn_df = pd.read_csv(cfg.data['gene4'] + '/in-mobi-newG4DN.csv', low_memory=False)
    out_mobi_g4dn_df = pd.read_csv(cfg.data['gene4'] + '/out-mobi-newG4DN.csv', low_memory=False)

    # (25961, 27)
    in_mobi_g4dn_df0 = in_mobi_g4dn_df[['index', 'acc_x', 'length', 'position_x', 'aa1', 'aa2', 'content_fraction',
                                        'content_count', 'refSeq', 'Gene_refGene', 'exon#', 'frameshift', 'Extreme',
                                        'Chr', 'Start', 'End', 'Ref', 'Alt', 'ExonicFunc.refGene',
                                        'GeneFunction.refGene', 'GeneExpressionTissue.refGene', 'GeneDisease.refGene',
                                        'InterVar_automated', 'Phenotype', 'Platform', 'Study', 'PubMed ID']]
    # (102495, 27)
    out_mobi_g4dn_df0 = out_mobi_g4dn_df[['index', 'acc_x', 'length', 'position_x', 'aa1', 'aa2', 'content_fraction',
                                          'content_count', 'refSeq', 'Gene_refGene', 'exon#', 'frameshift', 'Extreme',
                                          'Chr', 'Start', 'End', 'Ref', 'Alt', 'ExonicFunc.refGene',
                                          'GeneFunction.refGene', 'GeneExpressionTissue.refGene', 'GeneDisease.refGene',
                                          'InterVar_automated', 'Phenotype', 'Platform', 'Study', 'PubMed ID']]

    # (17329, 27)
    in_mobi_g4dn_df1 = in_mobi_g4dn_df0.drop_duplicates(subset=['acc_x', 'length', 'position_x', 'aa1', 'aa2',
                                                                'content_fraction', 'content_count',
                                                                'Gene_refGene', 'exon#', 'frameshift', 'Extreme',
                                                                'Chr', 'Ref', 'Alt',
                                                                'ExonicFunc.refGene', 'GeneFunction.refGene',
                                                                'GeneExpressionTissue.refGene', 'GeneDisease.refGene',
                                                                'InterVar_automated', 'Phenotype', 'Platform'])
    # (67822, 27)
    out_mobi_g4dn_df1 = out_mobi_g4dn_df0.drop_duplicates(subset=['acc_x', 'length', 'position_x', 'aa1', 'aa2',
                                                                  'content_fraction', 'content_count', 'Gene_refGene',
                                                                  'exon#', 'frameshift', 'Extreme', 'Chr', 'Ref', 'Alt', 'ExonicFunc.refGene',
                                                                  'GeneFunction.refGene',
                                                                  'GeneExpressionTissue.refGene',
                                                                  'GeneDisease.refGene', 'InterVar_automated',
                                                                  'Phenotype', 'Platform'])
    ## number of unique proteins
    in_mobi_uniprotid_lst = in_mobi_g4dn_df1['acc_x'].unique().tolist()  # 5047
    out_mobi_uniprotid_lst = out_mobi_g4dn_df1['acc_x'].unique().tolist()  # 9397
    all_g4dn_acc_lst = in_mobi_uniprotid_lst + out_mobi_uniprotid_lst  # 14444

    ##  if first dataset of proteins (1090) is inside g4dn proteins
    first_proteins_df = pd.read_csv(cfg.data[''] + '/allUniqueEntry.tab', sep='\t')
    first_proteins_lst = first_proteins_df['Entry'].tolist()

    mapped_pr_lst = []
    unmapped_pr_lst = []
    for i in first_proteins_lst:
        if i in all_g4dn_acc_lst:
            mapped_pr_lst.append(i)
        else:
            unmapped_pr_lst.append(i)

    ## check if positive candidate genes in mobi_g4dn merged files (both inside IDRs and outside)
    # positive genes imported from limitedmutations.py method genes_lst_maker()
    pos_candidate_gene_lst, ctrl_candidate_genes_lst = lmut.genes_lst_maker()  # 181

    ## this two new dfs could have redundant proteins (several prs from one gene can have muts in or out idr or both!)
    # (2143, 27) protein mutations in IDR
    pos_cand_idr_g4mobi_df = in_mobi_g4dn_df1[in_mobi_g4dn_df1.Gene_refGene.isin(pos_candidate_gene_lst)]
    pos_cand_idr_g4mobi_df = pos_cand_idr_g4mobi_df.drop_duplicates(
        subset=pos_cand_idr_g4mobi_df.columns.difference(['index']))
    # (6853, 27) protein mutations out IDR
    pos_cand_oidr_g4mobi_df = out_mobi_g4dn_df1[out_mobi_g4dn_df1.Gene_refGene.isin(pos_candidate_gene_lst)]
    pos_cand_oidr_g4mobi_df = pos_cand_oidr_g4mobi_df.drop_duplicates(
        subset=pos_cand_oidr_g4mobi_df.columns.difference(['index']))

    # control df
    # (first concat mobi_g4dn in and out dfs, then put ctrl in that df ), also same for pos_cand_genes_lst

    ## Unique proteins, total: 304, in:131, out: 173
    pos_cand_idr_pr_lst = pos_cand_idr_g4mobi_df['acc_x'].unique().tolist()
    pos_cand_oidr_pr_lst = pos_cand_oidr_g4mobi_df['acc_x'].unique().tolist()

    # TODO: change the names first, then make merged dfs with phens_pos_genes_df and g4mobi
sys.exit(main())
