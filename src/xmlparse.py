# This file is created to parse uniprot xml data file to get mut locations and other useful info (sep 14th)
import xml.etree.ElementTree as ET
from inspect import getmembers, isclass, isfunction
import config as cfg

NS = {'uniprot': 'http://uniprot.org/uniprot'}
tree = ET.parse(cfg.data['xml'] + '/uniprot_example.xml')
root = tree.getroot()
sequence_length = int(root.find("uniprot:entry/uniprot:sequence", NS).attrib["length"])
# for ele in root.findall("uniprot:entry/uniprot:dbReference", NS):
#     if ele.attrib["type"] == "PDB":
#         print(ele.attrib)
#         for ele2 in ele.findall("uniprot:property", NS):
#             if ele2.attrib["type"] == "chains":
#                 print(ele2.attrib)

# for ele in root.findall('uniprot:entry/uniprot:feature', NS):  # try this with tree instead of NS
#     if ele.attrib['type'] == 'sequence variant':
#         print(ele.attrib)
#         for ele2 in ele.findall('uniprot:location/uniprot:position', NS):
#             print(ele2.attrib)
