# created at sep 8th to extract variant id from uniprot txt API, store and insert in expasy API to
# get the mutation position fromn there
import urllib.request
import re
import pandas as pd
import config as cfg
import brain as bd

## example link: https://www.uniprot.org/uniprot/Q14204.txt
# my_proteins_acc_df = pd.read_csv(cfg.data['phens-fdr'] + '/acc-phen-5percentFDR.csv', usecols=['acc'])
# my_proteins_acc_lst = set(my_proteins_acc_df['acc'].tolist())  # 1308
# my_proteins_acc_lst = list(dict.fromkeys(my_proteins_acc_lst))
# my_proteins_acc_lst = my_proteins_acc_lst[:3]

my_proteins_acc_lst = ['Q14204', 'P14867']
proteins_dict = {}
par_list = []

for acc in my_proteins_acc_lst:
    url = 'https://www.uniprot.org/uniprot/' + acc + '.txt'
    # proteins_dict[acc] = []
    disease_par_lst = []
    file_test = urllib.request.urlopen(url)
    for l in file_test:
        line = l.decode("utf-8")
        if re.search(r'^CC\s+-!-', line):
            if re.search(r'^CC\s+-!- DISEASE:', line):
                # Crerate new paragraph
                disease_par_lst.append(line)
            # Otherwise, break only if all DISEASE line have been terminated
            elif len(disease_par_lst) > 0:
                break
        elif len(disease_par_lst) > 0:
            # print(line)
            disease_par_lst[-1] += line

    for par in disease_par_lst:
        print('Starting paragraph: ')
        # Extract subsequence
        par = re.sub(r'[\n\r]+', ' ', par)
        par = re.split(r'^CC\s+-!- DISEASE: ', par, 1)
        # par = re.split( r'\[MIM:\d{6}\]\s+CC', par[1], 1)
        par = re.split(r'\[MIM:\d{6}\]', par[1], 1)
        print('org: '+par[0])
        par = par[0]
        par = re.split(r'')
        #print(par) TODO: delete CC as well and replace with space
        # Retrieve list
        par_list = proteins_dict.setdefault(acc, [])
        # Store paragraph
        par_list.append(par)

# TODO: debug
print(proteins_dict)
print(par_list)
