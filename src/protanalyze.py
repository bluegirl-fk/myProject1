# this is to analyze proteins of Human, Brain and NDD phenotypes creating arrays of intersection and unition of these.
import pandas as pd
import config as cfg
import brain as brn

if __name__ == '__main__':

    phens_lst = ['Human', 'Brain', 'ASD', 'EE', 'ID', 'DD', 'SCZ', 'NDDs', 'Mix', 'Control']
    ndd_subdf = pd.read_csv(cfg.data['phens-fdr'] + '/acc-phen-5percentFDR.csv')
    ndd_subdf = ndd_subdf.drop_duplicates()  # (4531, 3)


    def phens_acc_lst_maker(phens_lst, phens_df):
        phen_dict = dict.fromkeys(phens_lst)  # this dict has phen_name as key and lst of corresponding ACCs as values
        for each_phen in phens_lst:
            tmp_lst = list(set(phens_df.loc[phens_df['Phenotype'] == each_phen, 'acc'].tolist()))
            phen_dict[each_phen] = tmp_lst
        return phen_dict
    def human_brain_acc_adder(phens_dict):
        mobidb = pd.read_csv(cfg.data[''] + '/mobidb_result.tsv', sep='\t', cols=['acc'])
    dict1 = phens_acc_lst_maker(phens_lst, ndd_subdf)
