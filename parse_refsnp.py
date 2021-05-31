import bz2
import json
import os

with open("data/genebank.tsv", "w") as f_gb, open("data/uniprot.tsv", "w") as f_unp:
    for name in os.listdir("/home/fatemeh/Desktop/ext-1tb/refsnp/"):
        if (name.startswith("refsnp") and ".json.bz2" in name):

            with bz2.open("/home/fatemeh/Desktop/ext-1tb/refsnp/" + name) as f:
                for i, line in enumerate(f):
                    if i % 1000000 == 0:
                        print("   {} lines processed of {}".format(i, name))
                    doc = json.loads(line)
                    rs_ids = ",".join([ele.get("merged_rsid") for ele in doc.get("dbsnp1_merges", [])])
                    if rs_ids:
                        for allele in doc.get("primary_snapshot_data", {}).get("placements_with_allele", []):
                            if allele.get("placement_annot", {}).get("mol_type") == "protein":
                                for alleles in allele.get("alleles", []):
                                    var = alleles.get("allele").get("spdi", {})
                                    if var:
                                        f_gb.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(doc["refsnp_id"], rs_ids, var.get("seq_id"), var.get("position"), var.get("deleted_sequence"), var.get("inserted_sequence")))

                        for allele in doc.get("primary_snapshot_data", {}).get("allele_annotations", []):
                            for clinical in allele.get("clinical", []):
                                for variant_identifier in clinical.get("variant_identifiers", []):
                                    if variant_identifier.get("organization") == "UniProtKB":
                                        f_unp.write("{}\t{}\t{}\n".format(doc["refsnp_id"], rs_ids, variant_identifier["accession"]))


