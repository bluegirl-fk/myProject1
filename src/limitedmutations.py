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

    positive_genes_df = candidate_gens_df.loc[candidate_gens_df['FDR'] <= 0.05]
    positive_genes_lst = positive_genes_df['Gene symbol'].unique().tolist()  # 181

    ctrl_lst = ['Control']
    control_df = all_candidate_genes_df[all_candidate_genes_df.Groups.isin(ctrl_lst)]  # (20907, 12)
    # ctrl_pos_df = control_df.loc[control_df['FDR'] <= 0.05] # this gives 0 results
    ctrl_pos_lst = control_df['Gene symbol'].unique().tolist()  # 20907
    return positive_genes_lst, ctrl_pos_lst


if __name__ == '__main__':

    ## Candidate genes from excel file (list is not complete)
    # genes_df = pd.read_excel(cfg.data['gene4'] + '/Gene.xlsx', engine='openpyxl')  # (8271, 13)


    sys.exit()
