# This file is created to parse uniprot xml data file to get mut locations and other useful info (sep 14th)
import xml.etree.ElementTree as ET
import config as cfg
import pandas as pd
import dateutil

NS = {'uniprot': 'http://uniprot.org/uniprot'}
# tree = ET.parse(cfg.data['xml'] + '/uniprot_example.xml')
tree = ET.parse(cfg.data['xml'] + '/uniprot-reflst-40001tilltheend.xml')
root = tree.getroot()
# sequence_length = int(root.find("uniprot:entry/uniprot:sequence", NS).attrib["length"])
entries = root.findall('uniprot:entry', NS)
acc_protinfo = []
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


df = pd.DataFrame(acc_protinfo, columns=['acc', 'id', 'description', 'aa', 'variation', 'pos'])
df.to_csv(cfg.data['phens'] + '/uniprot_variants-parsed-40001tillend.csv')
# TODO: check the final df against the ref lst and make sure all proteins are included, if not they should not be found in uniprot
# TODO: make dict of acc and corresponding diseases


