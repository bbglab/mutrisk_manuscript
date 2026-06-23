# load all neccesary functions and tables:
source("code/0_functions/analysis_variables.R")

genie_data = fread("raw_data/GENIE_17/data_mutations_extended.txt")
genie_metadata = fread("raw_data/GENIE_17/GENIE_release_17.0/data_clinical_sample.txt")
genie_samples = genie_metadata[`Cancer Type` %in% c("Colorectal Cancer", "Non-Small Cell Lung Cancer",  "Melanoma", "Leukemia"),] |>
  select(`Sample Identifier`, `Age at Which Sequencing was Reported (Years)`, `Cancer Type`,
         `Sequence Assay ID`) |>
  `colnames<-`(c("Tumor_Sample_Barcode", "age", "cancer_type", "panel"))

genie_cancer_type = inner_join(genie_samples, genie_data)
genie_crc = genie_cancer_type[cancer_type == "Colorectal Cancer", ]
genie_crc$Variant_Classification |> table()

# check which cohorts do not test for our genes of interest:
cohort_gene_files = list.files("raw_data/GENIE_17/GENIE_release_17.0/", pattern = "data_gene_panel_",
                               full.names = TRUE)

select_genes = c("TP53", "APC", "KRAS")
# only select the cohorts which contain TP53, APC and KRAS:
genomic_info = fread("raw_data/GENIE_17/GENIE_release_17.0/genomic_information.txt")
cohort_all_genes = genomic_info |> select(Hugo_Symbol, SEQ_ASSAY_ID) |>
  filter(Hugo_Symbol %in% select_genes) |> distinct() |>
  dplyr::count(SEQ_ASSAY_ID) |>
  filter(n == 3 )
missing_cohorts = setdiff(unique(genomic_info$SEQ_ASSAY_ID), cohort_all_genes$SEQ_ASSAY_ID)
genie_crc = genie_crc |>
  filter(panel %in% cohort_all_genes$SEQ_ASSAY_ID)

# check panel size of the cohorts:
genomic_info = genomic_info |>
  mutate(length = End_Position - Start_Position)

panel_sizes = genomic_info |>
  filter(includeInPanel == TRUE) |>
  group_by(SEQ_ASSAY_ID) |>
  summarize(panel_length = sum(length)) |>
  arrange(panel_length)  |>
  dplyr::rename(panel = SEQ_ASSAY_ID) |>
  filter(panel_length > 0)

panel_muts = genie_crc |>
  count(Tumor_Sample_Barcode, panel)

panel_rate = inner_join(panel_muts, panel_sizes) |>
  mutate(rate = n*1e6 / panel_length,
         hypermut = rate > 10)
hypermut_rate = panel_rate |> group_by(panel) |>
  summarize(n_normal = sum(!hypermut),
            n_hypermut = sum(hypermut))

ggplot(hypermut_rate, aes(x = n_normal, y = n_hypermut)) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = panel)) +
  geom_abline(slope = 1)

# for simplicity, we will filter for centers which have a 'normal' proportion of hypermut vs non-hypermut
genie_crc = genie_crc |> filter(Center %in% c("MSK", "DFCI", "PROV"))
# also, filter out any hypermutated samples
genie_crc = left_join(genie_crc, panel_rate)  |>
  filter(hypermut == FALSE)

genie_crc = genie_crc |>
  filter(Hugo_Symbol %in% select_genes) |>
  mutate(variant_class = case_match(Variant_Classification,
                                    c("Frame_Shift_Del", "Frame_Shift_Ins") ~ "fs",
                                    c("Splice_Region", "Splice_Site") ~ "splice",
                                    "Translation_Start_Site" ~ "Nonsense_Mutation",
                                    c("3'Flank", "3'UTR", "5'UTR", "5'Flank", "Intron", "Silent", "Nonstop_Mutation",
                                      "In_Frame_Del", "In_Frame_Ins") ~ "other",
                                    .default = Variant_Classification)) |>
  filter(variant_class != "RNA")

# number of CRC samples
genie_crc$Tumor_Sample_Barcode |> unique() |> length()

# number of mutated samples per gene
genie_crc |> filter(Hugo_Symbol %in% select_genes) |>
  select(Hugo_Symbol, Tumor_Sample_Barcode, hypermut) |> distinct() |>
  count(Hugo_Symbol, hypermut)

# Make empty matrix for all sample/gene/muttype combinations
ids = unique(genie_crc$Tumor_Sample_Barcode)
mut_df = data.frame(Tumor_Sample_Barcode = rep(ids, each = length(select_genes)),
                    Hugo_Symbol = rep(select_genes, length(ids)),
                    Nonsense_Mutation = 0, fs = 0, Missense_Mutation = 0,  other = 0, splice = 0) |>
  mutate(idcol = paste0(Tumor_Sample_Barcode, "_", Hugo_Symbol))

crc_mtype = genie_crc |>
  count(Tumor_Sample_Barcode, variant_class, Hugo_Symbol) |>
  pivot_wider(values_from = n, names_from = variant_class, values_fill = 0)  |>
  dplyr::rename(missense = Missense_Mutation, nonsense = Nonsense_Mutation) |>
  mutate(idcol = paste0(Tumor_Sample_Barcode, "_", Hugo_Symbol))

idx = match(crc_mtype$idcol, mut_df$idcol)
all(mut_df[idx, "idcol"] == crc_mtype$idcol) # check if the order is exactly the same
mut_df[idx, ] = crc_mtype # replace all values with mutated sites

# At the moment, we are not taking into account the 'other' mutations. The impact of these is hard to estimate.
major_types = mut_df |>
  mutate(mutation_type = case_when(Nonsense_Mutation + Missense_Mutation + splice == 2 & fs == 0 ~ "double_snv",
                                   Nonsense_Mutation + Missense_Mutation + splice == 1 & fs == 0 ~ "single_snv",
                                   Nonsense_Mutation + Missense_Mutation + splice > 2 & fs == 0 ~ ">=3snv",
                                   Nonsense_Mutation + Missense_Mutation + splice == 1 & fs == 1 ~ "snv+fs",
                                   Nonsense_Mutation + Missense_Mutation + splice == 0 & fs == 1 ~ "single_fs",
                                   Nonsense_Mutation + Missense_Mutation + splice == 0 & fs == 2 ~ "double_fs",
                                   Nonsense_Mutation + Missense_Mutation + splice == 0 & fs >= 2 ~ "fs>2",
                                   Nonsense_Mutation + Missense_Mutation + splice > 1 & fs > 0 ~ ">=2snv_fs",
                                   Nonsense_Mutation + Missense_Mutation + splice == 1 & fs > 1 ~ "single_snv_>=2fs",
                                   other != 0 ~ "other",
                                   Nonsense_Mutation + Missense_Mutation + splice + other + fs == 0 ~ "WT"))

# Compute the cumulative percentages (top of each rectangle)
colon_apc_mut_counts = major_types |>
  group_by(mutation_type, Hugo_Symbol) |>
  count()

colon_apc_mut_counts |>
  filter(mutation_type == "WT")



types = major_types |>
  dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol, mutation_type)

# fraction of samples with KRAS mutations
APC_double = types |>
  mutate(m  = paste0(paste0(Hugo_Symbol, mutation_type))) |>
  filter(m == "APCdouble_snv") |>
  dplyr::count(Tumor_Sample_Barcode)
sum(APC_double$n == 1) / 8621
# 11.85 of all CRCs have two APC snv mutations

# fraction of samples with KRAS mutations
KRAS = types |>
  mutate(m  = paste0(paste0(Hugo_Symbol, mutation_type))) |>
  filter(m == "KRASsingle_snv") |>
  dplyr::count(Tumor_Sample_Barcode)
sum(KRAS$n == 1) / 8621
# 42.86% of all CRC have KRAS SNVs

KRAS_APC_double = types |>
  mutate(m  = paste0(paste0(Hugo_Symbol, mutation_type))) |>
  filter(m %in% c("APCdouble_snv", "KRASsingle_snv")) |>
  dplyr::count(Tumor_Sample_Barcode)
# about 5% of all CRC samples have the double APC snv and KRAS mutations
sum(KRAS_APC_double$n == 2) / 8621
# 5.52% of all CRC have KRAS + double APC SNV

KRAS_APC_single  = types |>
  mutate(m  = paste0(paste0(Hugo_Symbol, mutation_type))) |>
  filter(m %in% c("APCsingle_snv", "KRASsingle_snv")) |>
  dplyr::count(Tumor_Sample_Barcode)
# about 5% of all CRC samples have the single APC snv and KRAS mutations
sum(KRAS_APC_single$n == 2) / 8621
# 9.30% of all CRC have KRAS + single APC SNV

# make a dataframe with the output of mutation type and fraction mutated:
mut_types = list(APC_single_snv = "APCsingle_snv",
                 KRAS_single_snv = "KRASsingle_snv",
                 APC_double = "APCdouble_snv",
                 KRAS_APC_single =  c("APCsingle_snv", "KRASsingle_snv"),
                 KRAS_APC_double =  c("APCdouble_snv", "KRASsingle_snv"))

fraction_mut = data.frame(n_mutated = rep(NA, length(mut_types)), n_CRC = NA, percentages = NA)
rownames(fraction_mut) = names(mut_types)
for (i in names(mut_types)) {

  n_crc = types$Tumor_Sample_Barcode |> unique() |> length()

  n_mutated = types |>
    mutate(m  = paste0(paste0(Hugo_Symbol, mutation_type))) |>
    filter(m %in% mut_types[[i]]) |>
    dplyr::count(Tumor_Sample_Barcode)  |>
    filter(n == length(mut_types[[i]]))
  # about 5% of all CRC samples have the single APC snv and KRAS mutations
  fraction_mut[i, "n_mutated"] = nrow(n_mutated)
  fraction_mut[i, "n_CRC"] = n_crc
  fraction_mut[i, "percentages"] = nrow(n_mutated) / n_crc
  # 9.30% of all CRC have KRAS + single A
}
fraction_mut = fraction_mut |> rownames_to_column("mutated_combination")
fwrite(fraction_mut, "processed_data/GENIE_17/CRC_mutation_fractions.txt")

# plot treemap indicating the fraction of mutated reads:
FS6B = colon_apc_mut_counts |>
  group_by(Hugo_Symbol) |>
  mutate(percentage = round(n*100 / sum(n),2)) |>
  ggplot(aes(area = n, fill = mutation_type)) +
  treemapify::geom_treemap(alpha = 0.7) +  # use the treemapify package to make a treeplot
  ggsci::scale_fill_igv() +
  facet_wrap(Hugo_Symbol ~ . ) +
  cowplot::theme_cowplot() +
  treemapify::geom_treemap_text(aes(label = paste0(mutation_type, ":\n", n, "\n", percentage, "%")),
                                color = "black",
                                place = "centre",
                                size = 15) +
  theme(legend.position = "none") +
  ggtitle("fraction mutated types in CRC - GENIE")

output_dir_s6 = "manuscript/figure_panels/figure_s6"
dir.create(output_dir_s6, showWarnings = FALSE, recursive = TRUE)

saveRDS(FS6B, "manuscript/figure_panels/figure_s6/Figure_S6B.rds")

# ============================================================
# Panel A: Polyp/adenoma incidence (moved from Figure 4)
# ============================================================
# Load polyp incidence data
ad_incidence = fread("raw_data/polyp_incidence/Adenomas_by10_intestines.csv")
age_min_column = rep(c(40, 50, 60, 70, 80), 2)
age_max_column = rep(c(50, 60, 70, 80, 90), 2)

# 135 CNADs sequenced, of which 73 have APC mutation
# 54% of adenomas have APC mutation
ad_incidence = ad_incidence |>
  filter(!Sex %in% c("Males Total", "Females Total")) |>
  mutate(min_age = age_min_column,
         max_age = age_max_column,
         age = (min_age + max_age) / 2,
         fraction_adenoma_apc = (`All sizes` * 0.54) / 10)

# take the mean across samples
ad_incidence_mean = ad_incidence |>
  group_by(age) |>
  summarize(fraction_adenoma_apc = mean(fraction_adenoma_apc), .groups = "drop")

# Create Panel A: Polyp incidence
FS6B = ad_incidence_mean |>
  ggplot(aes(x = age, y = fraction_adenoma_apc * 100)) +  # Convert to per 100
  geom_line(color = "#4a4a4a", linewidth = 1.2) +
  geom_point(color = "#4a4a4a", size = 2) +
  scale_y_continuous(labels = scales::label_number(), limits = c(0, 70)) +
  scale_x_continuous(limits = c(20, 85)) +
  theme_cowplot() +
  theme(
    axis.title.y = element_text(color = "#4a4a4a", size = 10),
    axis.text.y = element_text(color = "#4a4a4a", size = 9),
    axis.title.x = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(y = "Polyp incidence per 100 people", x = "Age (years)")

# Add panel label
# FS6B = FS6B + theme(plot.tag.position = c(-0.08, 1.03)) + labs(tag = "B")

saveRDS(FS6B, file.path(output_dir_s6, "Figure_S6A.rds"))

