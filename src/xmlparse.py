# This file is created to parse uniprot xml data file to get mut locations and other useful info (sep 14th)
import xml.etree.ElementTree as ET
import config as cfg

acc_protinfo_dic = {}
NS = {'uniprot': 'http://uniprot.org/uniprot'}
tree = ET.parse(cfg.data['xml'] + '/uniprot_example.xml')
root = tree.getroot()
# TODO: keep the length and compare with mobidb as a second checking measure of data accuracy
sequence_length = int(root.find("uniprot:entry/uniprot:sequence", NS).attrib["length"])
entries = root.findall('uniprot:entry', NS)
features = root.findall('uniprot:entry/uniprot:feature', NS)

# maybe should use var id as keys of child dicts
# data to be stored in the dictionary will be:
# {acc: {var_id: {description: xxxx, evidence:xxxxx, position:x, orig_aa:x, var_aa:x}}}
# (find a solution for entries without variant sequence data to be skipped) , in the end do it in terminal not pycharm

var_tmp_dict = {}
for entry in entries:
    acc = (entry.find('uniprot:accession', NS)).text
    acc_protinfo_dic[acc] = {}  # this will be changes later
    for feature in features:
        if feature.attrib['type'] == 'sequence variant':
            feature_tmp_dic = {}  # stores vars of a Pr, this child dict will be added to acc_protinfo parent dict
            var_id = feature.attrib['id']
            feature_tmp_dic[var_id] = {}
            description = feature.attrib['description']
            feature_tmp_dic[var_id].append(description)
            # acc_protinfo_dic[acc] = feature.attrib
            aa_orig = feature.find('uniprot:original', NS).text
            feature_tmp_dic[var_id].append(aa_orig)
            aa_var = feature.find('uniprot:variation', NS).text
            feature_tmp_dic[var_id].append(aa_var)
            for pos in feature.findall('uniprot:location/uniprot:position', NS):
                position = pos.attrib['position']
                feature_tmp_dic[var_id].append(position)  # check int/str based on expand regions method and mobidb data
    acc_protinfo_dic[acc] = feature_tmp_dic

print(acc_protinfo_dic)
print(len(acc_protinfo_dic))
print(acc_protinfo_dic.keys())
# add disease, then turn into a dataframe
