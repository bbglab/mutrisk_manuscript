# Script to produce figure 3, and the barplots indicating mutation accumulatoin in figure 4
library(MutationalPatterns)
library(ggh4x)
library(gt)
source("code/0_functions/analysis_variables.R")
getwd()

# load data sources
metadata_files = c("processed_data/blood/blood_metadata.tsv", "processed_data/colon/colon_metadata.tsv",
                   "processed_data/lung/lung_metadata.tsv")

names(metadata_files) = str_split_i(metadata_files, "\\/", 2)
metadata = lapply(metadata_files, \(x) fread(x)[,c("sampleID", "category", "age", "donor")]) |>
  rbindlist(idcol = "tissue")

# Load gene_of_interest boostdm
boostdm_files = list.files("processed_data/boostdm/boostdm_genie_cosmic/", pattern = "[lung|colon|CH]_boostDM_cancer.txt.gz", full.names = TRUE)
names(boostdm_files) = c("blood", "colon", "lung")
boostdm = lapply(boostdm_files,  \(x) fread(x) |> mutate(driver = ifelse(boostDM_class == TRUE, "driver", "non-driver")))# change names for overview


# load the mutation rates
expected_rate_list = list()
ratio_list = list()

for (tissue in c("colon", "blood", "lung")) {
  expected_rate_list[[tissue]] = fread(paste0("processed_data/", tissue, "/", tissue, "_expected_rates.tsv.gz"))
  ratio_list[[tissue]] = fread(paste0("processed_data/", tissue, "/", tissue, "_mut_ratios.tsv.gz"))
}
expected_rates = rbindlist(expected_rate_list, idcol = "tissue", use.names = TRUE)
ratios = rbindlist(ratio_list, idcol = "tissue", use.names = TRUE)

# filters
gene_of_interest = "TP53"

# Make a barplot showing the probabilitites for TP53 (poster usage)
include_hotspots = setNames(c(175, 248, 273,282), nm = c("R175",  "R248","R273","R282"))

prob_barplot_lung = make_gene_barplot(boostdm, ratios, expected_rates,  gene_of_interest = "TP53", tissue_select = "lung", category_select = "non-smoker",
                                      include_hotspots = include_hotspots, individual = "PD34215", cell_probabilities = TRUE) + labs(title = NULL, subtitle = "TP53 - lung", y = NULL)
prob_barplot_blood = make_gene_barplot(boostdm, ratios, expected_rates,  gene_of_interest = "TP53", tissue_select = "blood",
                                       include_hotspots = include_hotspots, individual = "KX008", cell_probabilities = TRUE) + labs(title = NULL, subtitle = "TP53 - blood", y = NULL)
prob_barplot_colon = make_gene_barplot(boostdm, ratios, expected_rates, gene_of_interest = "TP53", tissue_select = "colon",
                                       include_hotspots = include_hotspots, individual = "O340", cell_probabilities = TRUE) + labs(title = NULL, subtitle = "TP53 - colon")

F1B = wrap_plots(prob_barplot_colon, prob_barplot_lung, prob_barplot_blood, ncol = 3, guides = "collect") &
  theme(plot.subtitle = element_text(hjust = 0.5, vjust = 3.5))
saveRDS(F1B, "manuscript/figure_panels/figure_1/figure_1B.rds")

# Manuscript numbers
# calculate number of CpG vs non-CpG rate
# rates_CpG = expected_rates |>
#   left_join(triplet_match_substmodel) |>
#   left_join(metadata) |>
#   filter(tissue == "colon" & category == "normal") |>
#   mutate(cpg = ifelse(substr(triplet, 3,3) == "C" & substr(triplet, 7,7) == "G", "CpG", "non-CpG")) |>
#   group_by(sampleID, age, cpg, trinuc) |>
#   summarize(sum_rate = sum(mle), .groups = "drop_last") |>
#   summarize(mean_rate = mean(sum_rate)) |>
#   pivot_wider(names_from = cpg, values_from = mean_rate) |>
# Error: unexpected symbol in:
#  mutate(fold_change = CpG / `non-CpG`) |>
#   mutate(fold_change = CpG / `non-CpG`) |>
#   arrange(fold_change) |>
#   as.data.table()
# max CpG / non-CpG rate: 30, lowest is 6


# print mutation rate differences CpG vs non-CpG
cpg_muts = expected_rates |> left_join(triplet_match_substmodel) |>
  mutate(cpg = ifelse(substr(triplet, 3,3) == "C" & substr(triplet, 7,7) == "G", "CpG", "non-CpG")) |>
  group_by(cpg, trinuc) |>
  summarize(sum_rate = sum(mle), .groups = "drop_last") |>
  summarize(mean_rate = mean(sum_rate),
            min = min(sum_rate),
            max = max(sum_rate))
print(cpg_muts$mean_rate[1]/cpg_muts$mean_rate[2]) # 11.16
print(cpg_muts$min[1]/cpg_muts$min[2]) # 20.3
print(cpg_muts$max[1]/cpg_muts$max[2]) # 9.5

# Make a barplot indicating the number of mutations across TP53 across the three tissues (colon, lung, blood)
# Manuscript numbers also provided by this script
barplot_colon = make_gene_barplot(boostdm, ratios, expected_rates, gene_of_interest = "TP53", tissue_select = "colon",
                                  tissue_name = "Colon", cell_probabilities = FALSE) + labs(y = NULL, x = NULL)
barplot_lung = make_gene_barplot(boostdm, ratios, expected_rates, gene_of_interest = "TP53", tissue_select = "lung",
                                 tissue_name = "Lung", category_select = "non-smoker",
                                 cell_probabilities = FALSE) + labs(x = NULL)
barplot_blood = make_gene_barplot(boostdm, ratios, expected_rates, gene_of_interest = "TP53", tissue_select = "blood",
                                  tissue_name = "Blood", cell_probabilities = FALSE) + labs(y = NULL)

# first version of figure 3A
F3A = wrap_plots(barplot_colon, barplot_lung, barplot_blood, ncol = 1, guides = "collect")
saveRDS(F3A, "manuscript/figure_panels/figure_3/figure_3A.rds")

# update figure 3A with driver numbers:
driver_numbers = sapply(boostdm, \(x) x |> filter(gene_name == "TP53") |> pull(driver) |> table()   ) |>
  t() |> as.data.frame() |>
  mutate(label = paste0("  | TP53 driver SNVs: ", driver, " | non-driver SNVs: ", `non-driver`))
# reformat driver numbers:
driver_numbers = driver_numbers |> select(-label) |>
  rownames_to_column("Tissue") |>
  mutate(across(c(driver, `non-driver`), ~  format(., big.mark = ","))) |>
  `colnames<-`(c("Tissue", "TP53 driver SNV sites", "TP53 non-driver SNV sites")) |>
  mutate(Tissue = c("Blood", "Colon", "Lung"))

cols <- c("#ff725c", "#3ca951" , "#4269d0")
gt_driver_number = driver_numbers %>%  gt() %>%
  text_transform(locations = cells_body(columns = 1),
    fn = function(x) paste0('<span style="background: ', cols, '; color: white; padding: 4px 12px;
                            border-radius: 20px; font-weight: 500; display: inline-block;">',      x, '</span>')) %>%
  cols_width(1 ~ px(100), everything() ~ px(155)) |>
  cols_align("center", 1) |>  tab_options(column_labels.font.weight = "bold", data_row.padding = px(10))

gt_driver_number
gtsave(gt_driver_number, "manuscript/figure_panels/figure_3/figure_3B_table.png")

barplot_colon$labels$subtitle = paste0(barplot_colon$labels$subtitle, driver_numbers["colon", 3])
barplot_lung$labels$subtitle = paste0(barplot_lung$labels$subtitle, driver_numbers["lung", 3])
barplot_blood$labels$subtitle = paste0(barplot_blood$labels$subtitle, driver_numbers["blood", 3])
F3A = wrap_plots(barplot_colon, barplot_lung, barplot_blood, ncol = 1, guides = "collect")
saveRDS(F3A, "manuscript/figure_panels/figure_3/figure_3A.rds")

# Supplementary table for figure 3B: count unique TP53 genomic sites per tissue by driver/non-driver status
# Site = chr:pos:ref:alt (unique row in boostdm)
tp53_site_counts = lapply(names(boostdm), function(t) {
  tp53_data = boostdm[[t]][gene_name == "TP53"]
  data.table(
    tissue = t,
    driver_sites = tp53_data[driver == "driver", .N],
    nondriver_sites = tp53_data[driver == "non-driver", .N]
  )
}) |> rbindlist()

saveRDS(tp53_site_counts, "manuscript/figure_panels/figure_3/figure_3B_table.rds")

# Supplementary Figure 2 - TP53 for all tissues
tissue_categories = ratios |> select(tissue, category) |> distinct()
tissue_categories = tissue_categories[-5]
plot_list = list()
for (i in 1:nrow(tissue_categories)) {
  tissue_name = paste(as.character(tissue_categories[i]), collapse = "\n")
  plot_list[[i]] = make_gene_barplot(boostdm, ratios, expected_rates, gene_of_interest = "TP53",
                                     tissue_select = tissue_categories$tissue[i],
                                     category_select = tissue_categories$category[i],
                                     include_hotspots = include_hotspots,
                                    tissue_name = tissue_name, cell_probabilities = TRUE) +
    labs(title = NULL) + theme_classic()
}

plot_list = c(plot_list)
figure_S3 = wrap_plots(plot_list, nrow = 4) + plot_layout(guides = "collect") +
  plot_annotation(title = 'TP53: Expected number of mutated cells')
ggsave("manuscript/Supplementary_Figures/Figure_S3/Figure_S3.png", figure_S3, width = 14, height = 12)
ggsave("manuscript/Supplementary_Figures/Figure_S3/Figure_S3.svg", figure_S3, width = 14, height = 12)
ggsave("manuscript/Supplementary_Figures/Figure_S3/Figure_S3.pdf", figure_S3, width = 14, height = 12)

# APC colon barplot
APC_colon_normal = make_gene_barplot(boostdm, ratios, expected_rates, gene_of_interest = "APC",
                                     tissue_select = "colon", category_select = "normal", cell_probabilities = FALSE) +
  ggh4x::facet_grid2(driver ~ ., strip = strip_themed(background_y = elem_list_rect(fill = c("#C03830", "#707071")),
                                                      text_y = elem_list_text(colour = c("white"), face = "bold")), axes = "all",
                     remove_labels = "x")

colon_normal = make_gene_barplot(boostdm, ratios,expected_rates, gene_of_interest = "TP53",
                                 tissue_select = "colon", category_select = "normal", cell_probabilities = FALSE) +
  ggh4x::facet_grid2(driver ~ ., strip = strip_themed(background_y = elem_list_rect(fill = c("#C03830", "#707071")),
                                                      text_y = elem_list_text(colour = c("white"), face = "bold")), axes = "all",
                     remove_labels = "x")
ggsave("plots/colon/TP53_driver_non-driver.png", colon_normal, width = 10, height = 4.5, bg = "white")

# APC colon barplot
APC_colon_normal = make_gene_barplot(boostdm, ratios, expected_rates, gene_of_interest = "APC",
                                 tissue_select = "colon", category_select = "normal", cell_probabilities = FALSE) +
  ggh4x::facet_grid2(driver ~ ., strip = strip_themed(background_y = elem_list_rect(fill = c("#C03830", "#707071")),
                                                      text_y = elem_list_text(colour = c("white"), face = "bold")), axes = "all",
                     remove_labels = "x")

# Add dots to the colon barplot to make the individual bars more visible
df_dots = APC_colon_normal@data |>
  group_by(position, driver) |>
  summarize(mle = sum(mle))

APC_colon_normal = APC_colon_normal + geom_point(data = df_dots, aes(x = position, y = mle), size = 1.5) +
  scale_y_continuous(limits = c(NA, 2250), expand = expansion(mult = c(0, 0.1)))
ggsave("plots/colon/APC_driver_non-driver.png", APC_colon_normal, width = 12, height = 5, bg = "white")

# Figure 4A
F4A1 = make_gene_barplot(boostdm, ratios, expected_rates,  gene_of_interest = "APC", tissue_select = "colon", cell_probabilities = FALSE) +
  scale_y_continuous(breaks = extended_breaks(4), expand = expansion(mult = c(0, 0.1))) +
  scale_x_continuous(expand = c(0,0)) +
  theme(legend.position = "none",
        legend.text = element_text(size = rel(0.8)),
        legend.title = element_text(size = rel(0.8)),
        legend.key.size = unit(0.8, "lines"), legend.background = element_blank())

F4A2 = make_gene_barplot(boostdm, ratios, expected_rates, gene_of_interest = "KRAS", tissue_select = "colon", cell_probabilities = FALSE) +
  scale_y_continuous(breaks = extended_breaks(4), expand = expansion(mult = c(0, 0.1))) +
  theme(legend.position = "inside", legend.position.inside = c(0.9, 1),
        legend.text = element_text(size = rel(0.8)),
        legend.title = element_text(size = rel(0.8)),
        legend.key.size = unit(0.8, "lines"), legend.background = element_blank()) +
  scale_x_continuous(expand = c(0,0))

saveRDS(list(F4A1, F4A2), "manuscript/figure_panels/figure_4/figures_AB.rds")

# plot the number of mutations for TP53 as individual points
dotplot_list = list()
for (i in 1:nrow(color_df)) {

  category_select = color_df$category[i]
  tissue_select = color_df$tissue[i]

  # make dotplots with the individual variation:
  dotplot_list[[i]] = merge_mutrisk_drivers(boostdm, ratios, expected_rates, gene_of_interest = "TP53",
                                            tissue_select = tissue_select, category_select = category_select,
                                            filter_age = FALSE, individual = "all")[[1]] |>
    group_by(donor, driver) |>
    summarize(across(c(mle, cilow, cihigh), sum)) |>
    mutate(tissue = tissue_select, category = category_select)
}

dotplot_df = rbindlist(dotplot_list) |>
  mutate(tissue = factor(tissue, levels = c("colon", "lung", "blood")),
                         tissue_category = paste0(tissue, "_", category))

# Manuscript numbers: colon drivers TP53:
dotplot_df |>
  filter(tissue_category == "colon_normal" & driver == "driver") |>
  left_join(metadata) |> filter(age > 35) |>
  group_by(donor) |>
  summarize(across(c(mle, cilow, cihigh), mean), groups = "drop") |>
  summarize(min = min(mle), max = max(mle), mean = mean(mle))

# Manuscript numbers: colon drivers Lung TP53 |>
dotplot_df  |>
  filter(tissue == "lung" & driver == "driver") |>
  left_join(metadata) |>
  filter(age > 35) |>
  group_by(donor) |>
  summarize(across(c(mle, cilow, cihigh), mean), groups = "drop") |>
  summarize(min = min(mle), max = max(mle), mean = mean(mle))

dotplot_df  |>
  filter(tissue == "lung" & driver == "driver") |>
  filter(category == "non-smoker") |>
  left_join(metadata) |>
  filter(age > 35) |>
  group_by(donor) |>
  summarize(across(c(mle, cilow, cihigh), mean), groups = "drop") |>
  summarize(min = min(mle), max = max(mle), mean = mean(mle))

dotplot_df |>
  filter(tissue == "lung" & driver == "driver") |>
  left_join(metadata) |> filter(age > 35) |>
  group_by(tissue_category, donor) |>
  summarize(across(c(mle, cilow, cihigh), mean), groups = "drop") |>
  summarize(min = min(mle), max = max(mle), mean = mean(mle))

# Manuscript numbers: colon drivers Lung TP53 - median
dotplot_df |>
  filter(tissue == "lung" & driver == "driver") |>
  left_join(metadata) |> filter(age > 35) |>
  group_by(category, donor) |>
  summarize(across(c(mle, cilow, cihigh), mean), groups = "drop") |>
  summarize(median(mle))

df_total_muts = dotplot_df |>
  left_join(metadata |> select(-sampleID) |> distinct()) |>
  mutate(category = factor(category, levels = unique(color_df$category)))

tissue_plots = tissue_plots_raw =  list()

tissue_order = c("colon", "lung", "blood")
for (i in 1:3) {

  select_tissue = tissue_order[[i]]
  df_tissue = df_total_muts |>
    filter(tissue %in% select_tissue) |>
    filter(driver == "driver")

  plt = ggplot(df_tissue, aes(x = age, y = mle, color = tissue_category)) +
    geom_pointrange(aes(ymin = cilow, ymax = cihigh)) +
    facet_nested(. ~ tissue + category, axes = "y", remove_labels = "y") +
    scale_color_manual(values = tissue_category_colors) +
    scale_y_continuous(labels = scales::label_comma(), limits = c(0,NA)) +
    theme_cowplot() +
    theme(panel.grid = element_blank(), strip.background = element_blank(),
          legend.position = "none", ggh4x.facet.nestline = element_line()) +
    labs(x = "Age (years)", y = "Number of cells with\nTP53 driver mutation")
  tissue_plots_raw[[select_tissue]] = plt
  tissue_plots[[select_tissue]] = prep_plot(plt, LETTERS[i+2])
}

tissue_plots_raw[[2]] = tissue_plots_raw[[2]] + labs(y = NULL)
tissue_plots_raw[[3]] = tissue_plots_raw[[3]] + labs(y = NULL)

saveRDS(tissue_plots_raw, "manuscript/figure_panels/figure_3/figure_3CDE.rds")

# Figures for poster (can be removed if needed)
# figure exploration for poster:
poster = wrap_plots(tissue_plots_raw, nrow = 1, widths = c(4, 3, 1.2))
F3C = wrap_plots(tissue_plots, nrow = 1, widths = c(4, 3, 1.2))
F3C

# alternative poster figure [two rows of figures ]
F3C_bottom = wrap_plots(tissue_plots[-1], widths = c(2.8, 1))
F3C = tissue_plots[[1]] / F3C_bottom
F3C