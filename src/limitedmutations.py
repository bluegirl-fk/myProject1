# the aim of this file is to process candidate genes_1.2 file(7mb) and all de novo mutations_1.2 (59.6mb) from gene4dn
# Release version 1.1 (Updated: 09/09/2020)
# start date: July 30th 2021

import pandas as pd
import config as cfg
import sys


def main():
    print("Hello")
    return


def genes_lst_maker(fdr_limit):  # this is candidate genes from download page
    all_candidate_genes_df = pd.read_csv(cfg.data['gene4'] + '/limitedmut/Candidate_gene_1.2.txt',
                                         sep='\t')  # (124137, 12)

    positive_genes_df = all_candidate_genes_df.loc[all_candidate_genes_df['FDR'] <= fdr_limit]
    positive_genes_lst = positive_genes_df['Gene symbol'].unique().tolist()

    return positive_genes_lst


if __name__ == '__main__':
    ## Candidate genes from excel file (list is not complete)
    # genes_df = pd.read_excel(cfg.data['gene4'] + '/Gene.xlsx', engine='openpyxl')  # (8271, 13)
    # lst_1 = genes_lst_maker(0.05)  # 679 genes
    # lst_2 = genes_lst_maker(0.1)  # 0.02 = 1144, 0.1 = 842
    sys.exit()
