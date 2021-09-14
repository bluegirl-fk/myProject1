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

# TODO: this code is not correct cuz it needs to append several seq variant data to rach acc which is dict key, and in the end the mut position will be added to the dict
# maybe should use var id as keys of child dicts
# data to be stored in the dictionary will be:
# {acc: {var_id: {description: xxxx, evidence:xxxxx, position:x, orig_aa:x, var_aa:x}}}
for entry in entries:
    var_positions_lst = []
    acc = (entry.find('uniprot:accession', NS)).text
    acc_protinfo_dic[acc] = {}  # this will be changes later
    for feature in features:
        if feature.attrib['type'] == 'sequence variant':
            acc_protinfo_dic[acc] = feature.attrib
            for pos in feature.findall('uniprot:location/uniprot:position', NS):
                var_positions_lst.append(pos.attrib['position'])
            feature.attrib['var_position'] = var_positions_lst

print(acc_protinfo_dic)

# for ele in root.findall('uniprot:entry/uniprot:feature', NS):
#     if ele.attrib['type'] == 'sequence variant':
#         print(ele.attrib)
#         for ele2 in ele.findall('uniprot:location/uniprot:position', NS):
#             print(ele2.attrib)
