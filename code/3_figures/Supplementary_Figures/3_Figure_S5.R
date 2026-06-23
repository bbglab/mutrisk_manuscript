# supplementary figure 5:
library(Biostrings)
source("code/0_functions/analysis_variables.R")

# determine the whole-genome trinucleotide rates
sites_whole_genome = mutrisk::hg19_trinuc_counts

# trinucleotide sites in TP53
match = setNames(object = c("A", "C", "G", "T"), c("T", "G", "C", "A"))
cancer_boostdm = fread("processed_data/boostdm/boostdm_genie_cosmic/colon_boostDM_cancer.txt.gz") |>
  filter(gene_name == "TP53") |>
  mutate(trinuc = substr(mut_type, 1,3)) |>
  as.data.frame()

trinucs = DNAStringSet(cancer_boostdm$trinuc)
trinucs[substr(trinucs, 2,2) %in% c("A", "G")] = reverseComplement(trinucs[substr(trinucs, 2,2) %in% c("A", "G")])

trinucs = as.character(trinucs) |>  table()
trinuc_counts_TP53_all = data.frame(trinucleotide = names(trinucs),
                           TP_53_counts = as.numeric(trinucs))

# trinuc counts for driver mutations only
cancer_boostdm_driver = cancer_boostdm |> filter(boostDM_class == TRUE)
trinucs = DNAStringSet(cancer_boostdm_driver$trinuc)
trinucs[substr(trinucs, 2,2) %in% c("A", "G")] = reverseComplement(trinucs[substr(trinucs, 2,2) %in% c("A", "G")])

trinucs = as.character(trinucs) |>  table()
trinuc_counts_TP53_driver = data.frame(trinucleotide = names(trinucs),
                                    TP_53_driver_counts = as.numeric(trinucs))

trinuc_counts_all = left_join(sites_whole_genome, trinuc_counts_TP53_all) |>
  left_join(trinuc_counts_TP53_driver)
trinuc_counts_all |>
  pivot_longer(-trinucleotide) |>
  ggplot(aes(x = trinucleotide, y = value)) +
  geom_col() +
  facet_grid(name ~ . , scales = "free_y")

# get the trinucleotide mutation rates for POLD1 mutations
colon_muts = fread("processed_data/colon/colon_cell_muts.tsv.gz")
POLD1_muts = colon_muts |>
  filter(category == "POLD1")
POLD1_muts = get_mut_context(POLD1_muts)
POLD1_muts = POLD1_muts |> mutate(trinucleotide = trinuc)
POLD1_trinuc_counts = count(POLD1_muts, trinucleotide, name = "POLD1 mut count")

left_join(trinuc_counts_all, POLD1_trinuc_counts) |>
  pivot_longer(-trinucleotide) |>
  ggplot(aes(x = trinucleotide, y = value)) +
  geom_col() +
  facet_grid(name ~ . , scales = "free_y")


# make it slightly more specific: use the trinucloetide options:
triplets_TP53_all = cancer_boostdm |> select(-trinuc) |>
  left_join(triplet_match_substmodel) |>
  count(triplet, name = "TP53 full gene ")

triplets_TP53_driver = cancer_boostdm |> select(-trinuc) |>
  filter(boostDM_class == TRUE) |>
  left_join(triplet_match_substmodel) |>
  count(triplet, name = "TP53 driver mutations")

POLD1_counts = POLD1_muts |>
  count(triplet, name = "POLD1 overall\nmutation profile")
# 1. Define the label mapping
facet_labels <- c(
  `TP53 driver mutations` = "Mutations",
  `POLD1 overall\nmutation profile` = "Mutations",
  `trinucleotide presence\nwhole genome` = "Mutable sites"
)

figure_S5C = left_join(POLD1_counts, triplets_TP53_driver) |>
  left_join(triplet_match_substmodel |> select(-mut_type, -strand) |> distinct()) |>
  left_join(sites_whole_genome |> dplyr::rename(trinuc = trinucleotide, `trinucleotide presence\nwhole genome` = trinuc_counts)) |>
  pivot_longer(c(`TP53 driver mutations`, `POLD1 overall\nmutation profile`, `trinucleotide presence\nwhole genome`)) |>
  mutate(value = ifelse(is.na(value), 0, value)) |>
  mutate(triplet = factor(triplet, levels = TRIPLETS_96)) |>
  ggplot(aes(x = triplet, y= value, fill = type)) + 
  geom_col() +
  facet_grid(name ~ . , scales = "free_y", switch = "y", labeller = as_labeller(facet_labels)) + 
  cowplot::theme_cowplot() +
  scale_fill_manual(values = COLORS6) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        strip.text = element_text(size = 10),
        strip.placement = "outside",         
        strip.background = element_blank()) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = NULL)

saveRDS(figure_S5C, "manuscript/figure_panels/figure_s5/Figure_S5C.rds")

