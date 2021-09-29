# This file is created to parse uniprot xml data file to get mut locations and other useful info (sep 14th)
import xml.etree.ElementTree as ET
import config as cfg
import pandas as pd
import dateutil


def var_position_parser(xml_file_name):
    NS = {'uniprot': 'http://uniprot.org/uniprot'}
    # tree = ET.parse(cfg.data['xml'] + '/uniprot_example.xml')
    tree = ET.parse(cfg.data['xml'] + xml_file_name)
    root = tree.getroot()
    # sequence_length = int(root.find("uniprot:entry/uniprot:sequence", NS).attrib["length"])
    entries = root.findall('uniprot:entry', NS)
    acc_protinfo = []
    acc_no_seq_var = []
    for entry in entries:
        acc = (entry.find('uniprot:accession', NS)).text
        for feature in entry.findall('uniprot:feature', NS):
            if feature.attrib['type'] == 'sequence variant':
                for pos in feature.findall('uniprot:location/uniprot:position', NS):
                    if 'description' in feature.attrib:
                        if feature.find('uniprot:original', NS) is not None:
                            print('inside the if')
                            acc_protinfo.append([acc, feature.attrib['id'], feature.attrib['description'],
                                                 (feature.find('uniprot:original', NS)).text,
                                                 (feature.find('uniprot:variation', NS)).text,
                                                 pos.attrib['position']])
                        else:
                            print('Else!!!')
                            acc_protinfo.append([acc, feature.attrib['id'], feature.attrib['description'],
                                                 None, None, pos.attrib['position']])
                    else:
                        if feature.find('uniprot:original', NS) is not None:
                            print('inside the if')
                            acc_protinfo.append([acc, feature.attrib['id'], None,
                                                 (feature.find('uniprot:original', NS)).text,
                                                 (feature.find('uniprot:variation', NS)).text,
                                                 pos.attrib['position']])
                        else:
                            print('Else!!!')
                            acc_protinfo.append([acc, feature.attrib['id'], None,
                                                 None, None, pos.attrib['position']])

    return acc_protinfo


# acc_prot_info_lst = var_position_parser('/uniprot-reflst-40001tilltheend.xml')
# df = pd.DataFrame(acc_prot_info_lst, columns=['acc', 'var_id', 'description', 'orig_aa', 'var_aa', 'position'])
# df.to_csv(cfg.data['xml-p'] + '/uniprot_variants-parsed-40001tillend.csv')
## merging two dataframes:
# pr_vars_1to40k_df = pd.read_csv(cfg.data['xml-p'] + '/uniprot_variants-parsed-1to40000.csv')
# pr_vars_40Ktoend_df = pd.read_csv(cfg.data['xml-p'] + '/uniprot_variants-parsed-40001tillend.csv')
# pr_vars_df = pr_vars_1to40k_df.append(pr_vars_40Ktoend_df, ignore_index=True)
# del pr_vars_df['Unnamed: 0']
# pr_vars_df.to_csv(cfg.data['xml-p'] + '/protein-vars.csv')

# todo: make dict of acc and corresponding diseases
