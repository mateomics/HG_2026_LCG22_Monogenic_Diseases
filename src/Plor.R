library(tidyverse)

# Leer efectos
effects <- read.table(
  "results/effects.txt",
  header = FALSE,
  col.names = c("Effect")
)

# Clasificar regiones
effects$Region <- case_when(
  effects$Effect %in%
    c(
      "missense_variant",
      "stop_gained",
      "synonymous_variant",
      "frameshift_variant",
      "stop_lost",
      "start_lost",
      "protein_altering_variant",
      "coding_sequence_variant",
      "5_prime_UTR_variant",
      "3_prime_UTR_variant"
    ) ~ "Coding",

  effects$Effect %in%
    c(
      "intron_variant",
      "splice_acceptor_variant",
      "splice_donor_variant",
      "splice_region_variant"
    ) ~ "Intronic/Splicing",

  effects$Effect %in%
    c(
      "intergenic_region",
      "upstream_gene_variant",
      "downstream_gene_variant"
    ) ~ "Intergenic/Regulatory",

  TRUE ~ "Other"
)

# Contar
region_counts <- effects %>%
  count(Region) %>%
  arrange(desc(n))

region_counts
