import pandas as pd

input_vcf = "results/mendelian_annotated.vcf"
output_tsv = "data/processed/snpeff_annotations.tsv"

records = []

with open(input_vcf) as f:
    for line in f:

        if line.startswith("#"):
            continue

        cols = line.strip().split("\t")

        chrom = cols[0]
        pos = cols[1]
        ref = cols[3]
        alt = cols[4]
        info = cols[7]

        ann_field = None

        for item in info.split(";"):
            if item.startswith("ANN="):
                ann_field = item.replace("ANN=", "")
                break

        if ann_field is None:
            continue

        annotations = ann_field.split(",")

        for ann in annotations:

            ann_parts = ann.split("|")

            if len(ann_parts) < 5:
                continue

            allele = ann_parts[0]
            effect = ann_parts[1]
            impact = ann_parts[2]
            gene = ann_parts[3]

            records.append({
                "CHROM": chrom,
                "POS": pos,
                "REF": ref,
                "ALT": alt,
                "GENE": gene,
                "EFFECT": effect,
                "IMPACT": impact
            })

df = pd.DataFrame(records)

df.to_csv(output_tsv, sep="\t", index=False)

print(df.head())
print("\nTotal annotations:", len(df))