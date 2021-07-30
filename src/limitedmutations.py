# the aim of this file is to process candidate genes_1.2 file(7mb) and all de novo mutations_1.2 (59.6mb) from gene4dn
# Release version 1.1 (Updated: 09/09/2020)
# start date: July 30th 2021

import pandas as pd
import config as cfg
import sys


def main():
    print("Hello darkness my old friend..")
    return


def genes_lst_maker():
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


    sys.exit()
