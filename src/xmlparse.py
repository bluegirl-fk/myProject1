# This file is created to parse uniprot xml data file to get mut locations and other useful info (sep 14th)
import xml.etree.ElementTree as ET
import config as cfg

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
acc_protinfo_dic = {}
feature_tmp_dic = {}
var_tmp_dict = {}
for entry in entries:
    print(entry)
    acc = (entry.find('uniprot:accession', NS)).text
    acc_protinfo_dic[acc] = {}
    print('start!')
    print(acc)
    for feature in features:
        if feature.attrib['type'] == 'sequence variant':
            var_tmp_dict = {}
            var_id = feature.attrib['id']
            print(var_id)
            var_tmp_dict['description'] = feature.attrib['description']
            var_tmp_dict['aa_orig'] = feature.find('uniprot:original', NS).text
            var_tmp_dict['aa_var'] = feature.find('uniprot:variation', NS).text
            for pos in feature.findall('uniprot:location/uniprot:position', NS):
                var_tmp_dict['position'] = pos.attrib['position']  # check int/str based on expand regions method and mobidb data
                feature_tmp_dic[var_id] = var_tmp_dict
            print(len(feature_tmp_dic))

        elif feature.attrib['type'] != 'sequence variant' and len(feature_tmp_dic) > 0:
            print('End!')
            acc_protinfo_dic[acc] = feature_tmp_dic
            print(len(feature_tmp_dic), ' ==length')
            feature_tmp_dic = {}





print(acc_protinfo_dic)
print(len(acc_protinfo_dic))
print(acc_protinfo_dic.keys())
print(len(acc_protinfo_dic['P01911'].keys()))
print(len(acc_protinfo_dic['P04439'].keys()))
# # add disease, then turn into a dataframe
