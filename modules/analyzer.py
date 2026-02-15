from Bio import SeqIO

def analyze_vector_resistance(record):
    found = []
    for feature in record.features:
        if feature.type in ["gene", "CDS", "misc_feature"]:
            qualifiers = feature.qualifiers
            for tag in ["label", "gene", "note", "product"]:
                if tag in qualifiers:
                    name = str(qualifiers[tag][0]).lower()
                    for res in ["amp", "kan", "puro", "neo", "hyg", "bla"]:
                        if res in name: found.append(res)
    return list(set(found))