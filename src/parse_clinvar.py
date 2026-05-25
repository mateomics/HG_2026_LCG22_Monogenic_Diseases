import pandas as pd

input_vcf = "data/raw/pathogenic_only.vcf"

variants = []

with open(input_vcf) as f:
    for line in f:
        if line.startswith("#"):
            continue

        fields = line.strip().split("\t")

        chrom = fields[0]
        pos = fields[1]
        ref = fields[3]
        alt = fields[4]
        info = fields[7]

        info_dict = {}

        for item in info.split(";"):
            if "=" in item:
                key, value = item.split("=", 1)
                info_dict[key] = value

        variants.append({
            "CHROM": chrom,
            "POS": pos,
            "REF": ref,
            "ALT": alt,
            "GENE": info_dict.get("GENEINFO", "NA"),
            "DISEASE": info_dict.get("CLNDN", "NA"),
            "DISEASE_DB": info_dict.get("CLNDISDB", "NA"),
            "CLNSIG": info_dict.get("CLNSIG", "NA"),
            "MC": info_dict.get("MC", "NA")
        })

df = pd.DataFrame(variants)

df.to_csv(
    "data/processed/clinvar_pathogenic_parsed.tsv",
    sep="\t",
    index=False
)

print(df.head())
print("\nTotal variants:", len(df))