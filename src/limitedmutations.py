# the aim of this file is to process candidate genes_1.2 file(7mb) and all de novo mutations_1.2 (59.6mb) from gene4dn
# Release version 1.1 (Updated: 09/09/2020)
# start date: July 30th 2021

import pandas as pd
import config as cfg
import sys


def main():
    print("Hello darkness my old friend..")
    return


def genes_lst_maker():  # this is candidate genes from download page
    all_candidate_genes_df = pd.read_csv(cfg.data['gene4'] + '/limitedmut/Candidate_gene_1.2.txt',
                                         sep='\t')  # (124137, 12)
    phens_lst = ['ASD', 'ID', 'SCZ', 'EE']
    candidate_gens_df = all_candidate_genes_df[all_candidate_genes_df.Groups.isin(phens_lst)]  # (27844, 12)
    ctrl_lst = ['Control']
    control_df = all_candidate_genes_df[all_candidate_genes_df.Groups.isin(ctrl_lst)]  # (20907, 12)

    more_accurate_genes_df = candidate_gens_df.loc[candidate_gens_df['FDR'] <= 0.05]
    more_accurate_genes_lst = more_accurate_genes_df['Gene symbol'].unique().tolist()  # 181
    return more_accurate_genes_lst


if __name__ == '__main__':

    ## Candidate genes from excel file
    genes_df = pd.read_excel(cfg.data['gene4'] + '/Gene.xlsx', engine='openpyxl')  # (8271, 13)
    genes_df.to_csv(cfg.data['gene4'] + 'new_genes_file.csv')
    phens_lst = ['ASD', 'ID', 'SCZ', 'EE']
    phens_genes_df = genes_df[genes_df.group_name.isin(phens_lst)]  # (1815, 13)

    ctrl_lst = ['Control']
    control_genes_df = genes_df[genes_df.group_name.isin(ctrl_lst)]  # (1266, 13)

    ## unique proteins lists of phens and control
    fdr_based_genes_df = phens_genes_df.loc[phens_genes_df['FDR'] <= 0.05]
    # fdr_based_genes_lst = fdr_based_genes_df['gene_symbol'].unique.tolist()
    #
    # fdr_based_ctrl_df = control_genes_df.loc[phens_genes_df['FDR'] <= 0.05]
    # fdr_based_ctrl_lst = fdr_based_ctrl_df['gene_symbol'].unique().tolist()

    sys.exit()
