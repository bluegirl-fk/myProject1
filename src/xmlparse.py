# This file is created to parse uniprot xml data file to get mut locations and other useful info (sep 14th)
import xml.etree.ElementTree as ET
import config as cfg

NS = {'uniprot': 'http://uniprot.org/uniprot'}
tree = ET.parse(cfg.data['xml'] + '/uniprot_example.xml')
root = tree.getroot()
sequence_length = int(root.find("uniprot:entry/uniprot:sequence", NS).attrib["length"])
for ele in root.findall("uniprot:entry/uniprot:dbReference", NS):
    if ele.attrib["type"] == "PDB":
        print(ele.attrib)
        for ele2 in ele.findall("uniprot:property", NS):
            if ele2.attrib["type"] == "chains":
                print(ele2.attrib)
