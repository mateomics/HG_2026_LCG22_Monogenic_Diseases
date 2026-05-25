import pandas as pd

# Leer anotaciones SnpEff
df = pd.read_csv(
  "./data/processed/snpeff_annotations.tsv", 
  sep="\t", dtype={"CHROM": str}, low_memory=False)

# Jerarquía de impacto
impact_priority = {
  "HIGH": 1,
  "MODERATE": 2,
  "LOW": 3,
  "MODIFIER": 4
}

# Col. para ordenar
df["impact_rank"] = df["IMPACT"].map(impact_priority)

# Separar múltiples efectos
df["EFFECT"] = df["EFFECT"].str.strip().str.split("&") # SnpEff usa "&" como delimitador en efectos
df = df.explode("EFFECT").copy() # Cada fila con un efecto diferente

# Set de efectos de cada región
# Regiones codificantes
coding_effects = {
  "missense_variant",
  "synonymous_variant",
  "stop_gained",
  "frameshift_variant",
  "stop_lost",
  "start_lost",
  "protein_altering_variant",
  "coding_sequence_variant"
}

# Regiones no codificantes
utr_effects = {"3_prime_UTR_variant", "5_prime_UTR_variant"}

# Regiones intrónicas
intronic_effects = {"intron_variant"}

# Splicing
splicing_effects = {
  "splice_acceptor_variant",
  "splice_donor_variant",
  "splice_region_variant"
}

# Intergénicas y regulatorias
intergenic_effects = {
  "intergenic_region",
  "upstream_gene_variant",
  "downstream_gene_variant"
}

# Clasificar cada efecto en una región genómica
def classify_region(effect):
  if effect in coding_effects:
    return "Coding"
  elif effect in utr_effects:
    return "UTR"
  elif effect in intronic_effects:
    return "Intronic"
  elif effect in splicing_effects:
    return "Splicing"
  elif effect in intergenic_effects:
    return "Intergenic/Regulatory"
  else:
    return "Other"

df["REGION"] = df["EFFECT"].apply(classify_region)

# Jerarquía de severidad molecular
effect_priority = {
  "stop_gained": 1,
  "frameshift_variant": 2,
  "splice_acceptor_variant": 3,
  "splice_donor_variant": 4,
  "splice_region_variant": 5,
  "missense_variant": 6,
  "synonymous_variant": 7,
  "intron_variant": 8,
  "upstream_gene_variant": 9,
  "downstream_gene_variant": 10
}

# Rank intra-categoría
df["effect_rank"] = (df["EFFECT"].map(effect_priority).fillna(999))

# Ordenar por atributos, mayor impacto y efecto
df = df.sort_values(by=[
  "CHROM",
  "POS",
  "REF",
  "ALT",
  "impact_rank",
  "effect_rank"
  ])

# Una annotación por variante
df_unique = df.drop_duplicates(subset=["CHROM", "POS", "REF", "ALT"])

# Exportar dataset limpio
df_unique.to_csv(
  "data/processed/final_variant_annotations.tsv",
  sep="\t", index=False
)

# Conteo por región
region_summary = (
  df_unique["REGION"].value_counts().reset_index()
)
region_summary.columns = ["REGION", "COUNT"]

# Porcentaje
region_summary["PERCENT"] = (
  region_summary["COUNT"] /
  region_summary["COUNT"].sum()
) * 100

# Exportar resumen
region_summary.to_csv(
  "results/region_summary.tsv",
  sep="\t", index=False
)

print("\nResumen por región:")
print(region_summary)

print(df_unique.head())

print("\nTotal unique variants:")
print(len(df_unique))