# This file is created to parse uniprot xml data file to get mut locations and other useful info (sep 14th)
import xml.etree.ElementTree as ET
import config as cfg
import itertools
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


## PTM parser
NS = {'uniprot': 'http://uniprot.org/uniprot'}
tree = ET.parse(cfg.data['xml'] + '/uniprot_example.xml')
root = tree.getroot()
entries = root.findall('uniprot:entry', NS)
ptm_types = ['modified residue', 'glycosylation site', 'cross-link']  # disulfide bond not working
acc_ptm_lst = []
for entry in entries:
    acc = (entry.find('uniprot:accession', NS)).text
    for feature in entry.findall('uniprot:feature', NS):
        for ptm_type in ptm_types:
            if feature.attrib['type'] == ptm_type:
                for pos in feature.findall('uniprot:location/uniprot:position', NS):
                    acc_ptm_lst.append([acc, feature.attrib['description'], pos.attrib['position'], ptm_type])
            elif feature.attrib['type'] == 'disulfide bond':
                a = zip(feature.findall('uniprot:location/uniprot:begin', NS), feature.findall('uniprot:location/uniprot:end', NS))
                for (begin, end) in a:
                    acc_ptm_lst.append([acc, None, begin.attrib['position']+'&'+end.attrib['position'], 'disulfide bond'])

## todo add the end positionafter the begin, try to do both in one for loop, and concatenate them have them ass one value, eg: 46..79
#todo put an if for ptms if they have description or no

# df = pd.DataFrame(acc_ptm_lst, columns=['acc', 'description', 'pos', 'ptm_type'])
# df.to_csv(cfg.data['ptm-u'] + '/uniprot_parsed-ptms-40001tilltheend.csv')
# now I have two dfs to be concatenated into 1, but still lack the disulfide bond
