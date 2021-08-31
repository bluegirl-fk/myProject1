# this is to analyze proteins of Human, Brain and NDD phenotypes creating arrays of intersection and unition of these.
import numpy as np
import pandas as pd
import config as cfg
import brain as brn


def phens_acc_dict_maker(phen_lst, phens_df):  # TODO: later try to use these two methods inside eachother
    phen_dict = dict.fromkeys(phen_lst)  # this dict has phen_name as key and lst of corresponding ACCs as values
    for each_phen in phen_lst:
        tmp_lst = list(set(phens_df.loc[phens_df['Phenotype'] == each_phen, 'acc'].tolist()))
        phen_dict[each_phen] = tmp_lst
    return phen_dict


def human_brain_acc_adder(phens_dict):
    mobidb = pd.read_csv(cfg.data[''] + '/mobidb_result.tsv', sep='\t', usecols=['acc'])
    human_lst = list(set(mobidb['acc'].tolist()))
    brain_lst = brn.brain_pr_lst_generator()
    phens_dict['Human'] = human_lst
    phens_dict['Brain'] = brain_lst
    return phens_dict


if __name__ == '__main__':

    phens_lst = ['Human', 'Brain', 'ASD', 'EE', 'ID', 'DD', 'SCZ', 'NDDs', 'Mix', 'Control']
    ndd_subdf = pd.read_csv(cfg.data['phens-fdr'] + '/acc-phen-5percentFDR.csv')
    ndd_subdf = ndd_subdf.drop_duplicates()  # (4531, 3)
    # Dictionary with phen as key and their corresponding  list of ACCs as value
    phens_acc_dict = human_brain_acc_adder(phens_acc_dict_maker(phens_lst, ndd_subdf))





