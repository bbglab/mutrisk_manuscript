# Supplementary Note III Figures
library(GenomicRanges)
library(rtracklayer)
library(ggh4x)
source("code/0_functions/analysis_variables.R")

# load metadata
metadata_files = c("processed_data/blood/blood_metadata.tsv", "processed_data/colon/colon_metadata.tsv",
                   "processed_data/lung/lung_metadata.tsv")
names(metadata_files) = str_split_i(metadata_files, "\\/", 2)
metadata = lapply(metadata_files, \(x) fread(x)[,c("sampleID", "category", "age", "donor")]) |>
  rbindlist(idcol = "tissue")

# Load gene_of_interest boostdm
boostdm_files = list.files("processed_data/boostdm/boostdm_genie_cosmic/", pattern = "[lung|colon|CH]_boostDM_cancer.txt.gz", full.names = TRUE)
names(boostdm_files) = c("blood", "colon", "lung")
boostdm = lapply(boostdm_files,  \(x) fread(x) |> mutate(driver = ifelse(boostDM_class == TRUE, "driver", "non-driver")))

# load the mutation rates
expected_rate_list = list()
ratio_list = list()

for (tissue in c("colon", "blood", "lung")) {
  expected_rate_list[[tissue]] = fread(paste0("processed_data/", tissue, "/", tissue, "_expected_rates.tsv.gz"))
  ratio_list[[tissue]] = fread(paste0("processed_data/", tissue, "/", tissue, "_mut_ratios.tsv.gz"))
}
expected_rates = rbindlist(expected_rate_list, idcol = "tissue", use.names = TRUE)
ratios = rbindlist(ratio_list, idcol = "tissue", use.names = TRUE)

# plot the number of mutations for TP53 as individual points
tissue_select = "colon"
gene_of_interest = "TP53"

get_n_muts = function(tissue_select, gene_of_interest) {
  dotplot_list = list()
  categories = expected_rates |>
    filter(tissue == tissue_select) |>
    pull(category) |> unique()

  for (category_select in categories) {

    print(category_select)
     if(!gene_of_interest %in% boostdm[[tissue_select]]$gene_name) {stop("gene of interest not in BoostDM cancer type database")}

    # make dotplots with the individual variation:
    dotplot_list[[category_select]] = merge_mutrisk_drivers(boostdm, ratios, expected_rates, gene_of_interest = gene_of_interest,
                                              tissue_select = tissue_select, category_select = category_select,filter_age = FALSE,
                                              individual = "all")[[1]] |>
      group_by(donor, driver) |>
      summarize(across(c(mle, cilow, cihigh), sum)) |>
      mutate(tissue = tissue_select, category = category_select)
  }
  mutated_list = rbindlist(dotplot_list) |>
    mutate(gene = gene_of_interest) |>
    left_join(metadata |> select(donor, age) |> distinct()) |>
    mutate(tissue_category = paste0(tissue, "_", category))
  return(mutated_list)
}



plot_n_muts = function(mutated_list, driver = TRUE) {

  if (driver == TRUE) {
    mutated_list = mutated_list |> filter(driver == "driver")
  }

  order_categories = color_df |>
    filter(tissue == unique(mutated_list$tissue))

  gene = unique(mutated_list$gene)
  mutated_list = mutated_list |>
    mutate(category = factor(category, levels = order_categories$category),
           tissue_category = factor(tissue_category, levels = order_categories$tissue_category))
  plt = ggplot(mutated_list, aes(x = age, y = mle, color = tissue_category)) +
    geom_pointrange(aes(ymin = cilow, ymax = cihigh)) +
    facet_nested(. ~ tissue + category, axes = "y", remove_labels = "y") +
    scale_color_manual(values = tissue_category_colors) +
    scale_y_continuous(labels = scales::label_comma(), limits = c(0,NA)) +
    theme_cowplot() +
    labs(x = "Age (years)",
    y = paste0("Number of cells with\n", gene,  " driver mutation"), title = gene) +
    theme(panel.grid = element_blank(), strip.background = element_blank(),
          legend.position = "none", ggh4x.facet.nestline = element_line(),
          plot.title = element_text(hjust = 0.5))
  plt
}

mutated_list = get_n_muts("colon", "SMAD4")

SMAD4 = get_n_muts("colon", "SMAD4") |> plot_n_muts(driver = TRUE)
BRAF = get_n_muts("colon", "BRAF") |> plot_n_muts(driver = TRUE)
PIK3CA = get_n_muts("colon", "PIK3CA") |> plot_n_muts(driver = TRUE)

colon_plot = SMAD4 / BRAF / PIK3CA

output_dir = "manuscript/Supplementary_notes/Supplementary_Note_III/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
ggsave("manuscript/Supplementary_notes/Supplementary_Note_III/SN_II_Fig1_colon.png",
       colon_plot, width = 10, height = 14)

EGFR = get_n_muts("lung", "EGFR") |> plot_n_muts(driver = TRUE) |> prep_plot("A")
KRAS = get_n_muts("lung", "KRAS") |> plot_n_muts(driver = TRUE) |> prep_plot("B")
TP53 = get_n_muts("lung", "BRAF") |> plot_n_muts(driver = TRUE) |> prep_plot("C")
lung_plot = EGFR / KRAS / TP53

ggsave("manuscript/Supplementary_notes/Supplementary_Note_III/SN_II_Fig2_lung.png",
       lung_plot, width = 9, height = 14)


PPM1D = get_n_muts("blood", "PPM1D") |> plot_n_muts(driver = TRUE)
ASXL1 = get_n_muts("blood", "ASXL1") |> plot_n_muts(driver = TRUE)
SRSF2 = get_n_muts("blood", "SRSF2") |> plot_n_muts(driver = TRUE)
blood_plot = PPM1D | ASXL1 | SRSF2
ggsave("manuscript/Supplementary_notes/Supplementary_Note_III/SN_II_Fig3_blood.png",
       blood_plot, width = 12, height = 5)

