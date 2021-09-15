# This file is created to parse uniprot xml data file to get mut locations and other useful info (sep 14th)
import xml.etree.ElementTree as ET
import config as cfg

NS = {'uniprot': 'http://uniprot.org/uniprot'}
tree = ET.parse(cfg.data['xml'] + '/uniprot_example.xml')
root = tree.getroot()
# sequence_length = int(root.find("uniprot:entry/uniprot:sequence", NS).attrib["length"])
entries = root.findall('uniprot:entry', NS)
# {acc: {var_id: {description: xxxx, evidence:xxxxx, position:x, orig_aa:x, var_aa:x}}}
# (find a solution for entries without variant sequence data to be skipped)
acc_protinfo_dic = {}
feature_tmp_dic = {}
var_tmp_dict = {}
for entry in entries:
    acc = (entry.find('uniprot:accession', NS)).text
    acc_protinfo_dic[acc] = {}
    for feature in entry.findall('uniprot:feature', NS):
        if feature.attrib['type'] == 'sequence variant':
            var_tmp_dict = {}
            var_id = feature.attrib['id']
            var_tmp_dict['description'] = feature.attrib['description']
            var_tmp_dict['aa_orig'] = feature.find('uniprot:original', NS).text
            var_tmp_dict['aa_var'] = feature.find('uniprot:variation', NS).text
            for pos in feature.findall('uniprot:location/uniprot:position', NS):
                var_tmp_dict['position'] = pos.attrib[
                    'position']  # check int/str based on expand regions method and mobidb data
                feature_tmp_dic[var_id] = var_tmp_dict
        elif feature.attrib['type'] != 'sequence variant' and len(feature_tmp_dic) > 0:
            acc_protinfo_dic[acc] = feature_tmp_dic
            feature_tmp_dic = {}

print(acc_protinfo_dic)
acc_pos_d = {}
acc_lst = list(acc_protinfo_dic.keys())
for acc in acc_lst:
    vars_per_each_acc_lst = list(acc_protinfo_dic[acc].keys())
    muts_pos_per_each_acc_lst = []
    for id in vars_per_each_acc_lst:
        mut_position = acc_protinfo_dic[acc][id]['position']
        muts_pos_per_each_acc_lst.append(mut_position)
    acc_pos_d[acc] = muts_pos_per_each_acc_lst





# TODO: make dict of acc and corresponding diseases, for the ones not having that section, put None.
#  then do the mut in IDR analysis again
# in the end do it in terminal not pycharm
