from Bio import SeqIO
import json

class ReferenceParser():
    """
    Parses a genbank file containing multiple references and searches for
    annotated CDS features. Writes json file with the relevant information.
    """

    def __init__(self):
        self.reference_genbank = None
        self.reference_json = {"schema": "v1.1", 'features': {}, 'references': {}}

    def load_reference_genbank(self, filepath):
        # tests:
        # what if no file
        # simple case check
        # expected attributes (reference,sequence), (genes,(start,end,strand))
        self.reference_genbank = SeqIO.index(filepath, "genbank")

    def identify_features(self, num_features=None):
        for i in self.reference_genbank.keys():
            # print(i)
            locations = {}
            length = None
            for feature in self.reference_genbank[i].features:
                if feature.type == 'source':
                    length = int(feature.location.end)
                    continue
                if feature.type == 'gene':
                    continue
                if feature.type == 'CDS':
                    # print(feature)
                    if "gene" in feature.qualifiers:
                        gene_id = feature.qualifiers['gene'][0]
                        self.reference_json['features'][gene_id] = {
                            "name": gene_id.lower(),
                            "type": "CDS"
                        }
                        if "note" in feature.qualifiers:
                            self.reference_json['features'][gene_id]["description"] = feature.qualifiers['note'][0],
                        locations[feature.qualifiers['gene'][0]] = {
                            "start": int(feature.location.start),
                            "end": int(feature.location.end),
                            "strand": feature.location.strand
                        }
                # else:
                #    print(feature.type)
                #    print(feature)
            record = self.reference_genbank[i]
            self.reference_json['references'][i] = {
                'accession': i,
                'description': record.description,
                'length': length,
                'locations': locations,
                'sequence': str(record.seq)
            }

    def write_json(self, filepath):
        with open(filepath, 'w') as json_file:
            json.dump(self.reference_json, json_file)

    def run(self, infilepath, outfilepath, num_features=None):
        self.load_reference_genbank(infilepath)
        self.identify_features(num_features)
        self.write_json(outfilepath)