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


def df_col_selecter(inputdf):
    inputdf1 = inputdf[['index', 'acc_x', 'length', 'position_x', 'aa1', 'aa2', 'content_fraction',
                                        'content_count', 'refSeq', 'Gene_refGene', 'exon#', 'frameshift', 'Extreme',
                                        'Chr', 'Start', 'End', 'Ref', 'Alt', 'ExonicFunc.refGene',
                                        'GeneFunction.refGene', 'GeneExpressionTissue.refGene', 'GeneDisease.refGene',
                                        'InterVar_automated', 'Phenotype', 'Platform', 'Study', 'PubMed ID']]
    return inputdf1


def col_based_drop_duplicate(inputdf):
    inputdf1 = inputdf.drop_duplicates(subset=['acc_x', 'length', 'position_x', 'aa1', 'aa2',
                                                                'content_fraction', 'content_count', 'Gene_refGene',
                                                                'exon#', 'frameshift', 'Extreme', 'Chr', 'Ref', 'Alt',
                                                                'ExonicFunc.refGene', 'GeneFunction.refGene',
                                                                'GeneExpressionTissue.refGene', 'GeneDisease.refGene',
                                                                'InterVar_automated', 'Phenotype', 'Platform'])
    return inputdf1


if __name__ == '__main__':
    # mobi_mut_in_idr_df = pd.read_csv(cfg.data['gene4'] + '/mobidb-mut-pos-true.csv')  # (994383, 11)
    # mobi_mut_out_idr_df = pd.read_csv(cfg.data['gene4'] + '/mobidb-mut-pos-false.csv')  # (3238392, 11)
    # # (25961, 11)
    # mut_in_mobi_lite_df = mobi_mut_in_idr_df.loc[mobi_mut_in_idr_df['feature'] == 'prediction-disorder-mobidb_lite']
    # # (102495, 11)
    # mut_out_mobi_lite_df = mobi_mut_out_idr_df.loc[mobi_mut_out_idr_df['feature']=='prediction-disorder-mobidb_lite']

    ## merge with g4dn (proteins mutation position checked)

    in_mobi_g4dn_df = pd.read_csv(cfg.data['gene4'] + '/in-mobi-newG4DN.csv', low_memory=False)
    out_mobi_g4dn_df = pd.read_csv(cfg.data['gene4'] + '/out-mobi-newG4DN.csv', low_memory=False)
    ## limiting the columns
    in_mobi_g4dn_df0 = df_col_selecter(in_mobi_g4dn_df)  # (25961, 27)
    out_mobi_g4dn_df0 = df_col_selecter(out_mobi_g4dn_df)  # (102495, 27)
    ## drop duplicates except columns: index, refSeq, Start, End, Study, PubmedID
    # (17329, 27)
    in_mobi_g4dn_df1 = col_based_drop_duplicate(in_mobi_g4dn_df0)
    # (67822, 27)
    out_mobi_g4dn_df1 = col_based_drop_duplicate(out_mobi_g4dn_df0)

    ## add a new column with boolean values for in IDR or not
    in_mobi_g4dn_df1['is_in_idr'] = 1  # (17329, 28)
    out_mobi_g4dn_df1['is_in_idr'] = 0  # (67822, 28)

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
    # positive candidate genes are the ones with FDR =< 0.05
    pos_candidate_gene_lst = lmut.genes_lst_maker(0.05)  # 181

    ## this two new dfs could have redundant proteins (several prs from one gene can have muts in or out idr or both!)
    # (2143, 27) protein mutations in IDR
    pos_cand_idr_g4mobi_df = in_mobi_g4dn_df1[in_mobi_g4dn_df1.Gene_refGene.isin(pos_candidate_gene_lst)]
    pos_cand_idr_g4mobi_df = pos_cand_idr_g4mobi_df.drop_duplicates(
        subset=pos_cand_idr_g4mobi_df.columns.difference(['index']))
    # (6853, 27) protein mutations out IDR
    pos_cand_oidr_g4mobi_df = out_mobi_g4dn_df1[out_mobi_g4dn_df1.Gene_refGene.isin(pos_candidate_gene_lst)]
    pos_cand_oidr_g4mobi_df = pos_cand_oidr_g4mobi_df.drop_duplicates(
        subset=pos_cand_oidr_g4mobi_df.columns.difference(['index']))

    ## concat pos_cand_idr_g4mobi_df and pos_cand_oidr_g4mobi_df
    positive_cand_g4mobi_df = pos_cand_idr_g4mobi_df.append(pos_cand_oidr_g4mobi_df, sort=False)  # (8996, 28)
    positive_cand_g4mobi_df.to_csv(cfg.data['gene4'] + '/positive_cand_g4mobi_concat.csv')

    # control df
    # (first concat mobi_g4dn in and out dfs, then put ctrl in that df ), also same for pos_cand_genes_lst

    ## Unique proteins, total: 304, in:131, out: 173
    pos_cand_idr_pr_lst = pos_cand_idr_g4mobi_df['acc_x'].unique().tolist()
    pos_cand_oidr_pr_lst = pos_cand_oidr_g4mobi_df['acc_x'].unique().tolist()


    sys.exit(main())
