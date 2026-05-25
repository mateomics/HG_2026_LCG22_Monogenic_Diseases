import pandas as pd

df = pd.read_csv(
    "data/processed/clinvar_pathogenic_parsed.tsv",
    sep="\t"
)

mendelian = df[
    df["DISEASE_DB"].str.contains(
        "OMIM:",
        na=False
    )
]

print(mendelian.head())
print("\nTotal mendelian variants:", len(mendelian))

mendelian.to_csv(
    "data/processed/mendelian_variants.tsv",
    sep="\t",
    index=False
)