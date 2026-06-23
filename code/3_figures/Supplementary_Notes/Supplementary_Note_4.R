# Colon Figure 2 script
library(GenomicRanges)
library(rtracklayer)
source("code/0_functions/analysis_variables.R")

tissue = "colon"
ncells = tissue_ncells_ci$mid_estimate[1]
colon_colors = c(normal = "#3ca951", IBD = "#6cc5b0", POLE = "#145220", POLD1 = "#222e24")

# load colon metadata
metadata = fread(paste0("processed_data/", tissue, "/", tissue, "_metadata.tsv")) |>   distinct() |>
  mutate(category = factor(category, levels = c("normal", "IBD", "POLD1", "POLE"))) |>
  select(-sensitivity, -coverage) |> distinct()

# load the mutation rates
expected_rates = fread(paste0("processed_data/", tissue, "/", tissue, "_expected_rates.tsv.gz"))
ratios = fread(paste0("processed_data/", tissue, "/", tissue, "_mut_ratios.tsv.gz"))

# load the mutation data:
cancer_bDM = fread("processed_data/boostdm/boostdm_genie_cosmic/colon_boostDM_cancer.txt.gz")

# Get the actual values for the graphic to be correct
APC_1450_hotspot = cancer_bDM[gene_name == "APC" & aachange == "R1450*" , c("gene_name", "mut_type", "aachange", "position", "driver")]
expected_rate_APC_1450 = expected_rates |>
  inner_join(APC_1450_hotspot, by = "mut_type") |>
  left_join(metadata) |>
  left_join(ratios) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ratio * ncells)) |>
  group_by(donor, category) |>
  summarise(across(c(mle, cilow, cihigh, age), mean)) |>
  setDT()

# Values for in the manuscript: Individuals above 35 - APC 1450* mutation rate
expected_rate_APC_1450[category == "normal" & age > 35 , mle] |> mean()
setDT(expected_rate_APC_1450)[category == "normal" & age > 35, mle] |> min()
setDT(expected_rate_APC_1450)[category == "normal" & age > 35, mle] |> max()

#### TP53 mutations for colon:
# profile of the driver rates for colon:
TP53_muts = expected_rates |>
  group_by(category, mut_type) |>
  summarize(across(c(mle, cilow, cihigh), mean)) |>
  left_join(ratios[gene_name == "TP53"]) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ncells * ratio)) |>
  left_join(triplet_match_substmodel) |>
  left_join(cancer_bDM[gene_name == "TP53", ], by = c("mut_type", "gene_name"), relationship = "many-to-many") |>
  mutate(driver_status = ifelse(driver, "driver", "non-driver"))  |> setDT()

TP53_count = TP53_muts[, .(mle = sum(mle)) , by = c("position", "category", "driver_status") ]
normal_TP53 = TP53_muts[category == "normal",]
normal_TP53_count = TP53_count[category == "normal",]

label_df = data.frame(label = "TP53: expected number of crypt stem cells mutated",
                      driver_status = "driver")

expected_TP53 = ggplot(normal_TP53, aes(x = position, y = mle)) +
  geom_col(aes(fill = type)) +
  geom_point(data = normal_TP53_count) +
  ggpp::geom_text_npc(data = label_df,  aes(label = label),  npcx = 0.05, npcy = 0.96) +
  scale_fill_manual(values = COLORS6) +
  labs(y = "numer of mutated cells") +
  facet_grid(driver_status ~ .) +
  theme_cowplot() + panel_border() +
  scale_y_continuous(expand=expansion(mult=c(0,0.1)), labels = label_comma())
ggsave("plots/colon/TP53_expected.png", expected_TP53 ,  width = 15, height = 5.5, bg = "white")

#####
# get the percentage of mutated cells for APC
####
site_freqs_APC_drivers = cancer_bDM[gene_name == "APC" & driver == TRUE,  .N, by = c("mut_type", "gene_name")]

# percentage all drivers
all_drivers = cancer_bDM[driver == TRUE,  .N, by = c("mut_type", "gene_name")]
driver_rates  = get_gene_rate(exp_rates = expected_rates, metadata = metadata, site_freqs = all_drivers, ratios = ratios, ncells = 1)
mean_mutated = driver_rates |> group_by(category, gene_name) |>
  group_by(category, gene_name) |>
  filter(age > 18) |>
  summarize(across(c(mle, cilow, cihigh), mean)) |>
  arrange(mle) |>
  mutate(gene_name = factor(gene_name, unique(gene_name)))

ggplot(mean_mutated, aes(x = fct_reorder(gene_name, mle, .desc = TRUE), y = mle)) +
  geom_col() +
  facet_grid(category ~ ., scales = "free_y") +
  labs(x = NULL, y = "percent of cells carrying a mutation") +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.1))) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# adress comment ferran - inclusion exclusion:
probs = expected_rates |>
  left_join(distinct(metadata), by = c("category", "sampleID")) |>
  inner_join(all_drivers, by = "mut_type", relationship = "many-to-many") |>
  left_join(ratios, by = c("category", "gene_name")) |>
  group_by(category, mut_type, donor, age, ratio, gene_name, N) |>
  summarize(across(c(mle, cilow, cihigh),  ~ mean(. * ratio))) |>
  filter(!is.na(donor)) |>
  as.data.table()

Rcpp::sourceCpp("code/0_functions/inclusion_exclusion2.cpp")

fraction_mut_test = tibble(donor = unique(probs$donor), `inclusion exclusion` = NA, sum = NA)
for (select_donor in unique(probs$donor)) {
  print(select_donor)
  probs_donor = probs[donor %in% select_donor, ]
  seq = rep(probs_donor$mle, probs_donor$N)
  fraction_mut_test[fraction_mut_test$donor == select_donor, 2 ]  = inclusion_exclusion2(seq)
  fraction_mut_test[fraction_mut_test$donor == select_donor, 3 ] =  sum(seq)
}

fraction_mut_test = fraction_mut_test |>
  mutate(`fraction multiple muts` = (sum - `inclusion exclusion`) / sum) |>
  left_join(metadata |> select(donor, category) |> distinct())

barplot_ie_vs_sum = fraction_mut_test |> pivot_longer(c(`inclusion exclusion`, sum), values_to = "percent mutated", names_to = "method") |>
  ggplot(aes(x = donor, y = `percent mutated`, alpha = method,  fill = category)) +
  geom_col(position = "dodge", color = "black") +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = colon_colors) +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.1))) +
  theme_cowplot()  +
  theme(axis.text.x = element_blank(), legend.position = "inside", legend.position.inside = c(0.05, 0.75))  +
  labs()

percent_double_mut = ggplot(fraction_mut_test, aes(x = donor, y = `fraction multiple muts`, fill = category)) +
  geom_col() +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = colon_colors) +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.1))) +
  theme_cowplot()  +
  theme(axis.text.x = element_blank(), legend.position = "none")

percent_non_double_mut =  ggplot(fraction_mut_test, aes(x = donor, y = 1 - `fraction multiple muts`, fill = category)) +
  geom_col() +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = colon_colors) +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.1))) +
  theme_cowplot()  +
  theme(axis.text.x = element_blank(), legend.position = "none")

mg = 7
supplementary_note_II_figure = barplot_ie_vs_sum / (percent_double_mut | percent_non_double_mut )  +
  plot_annotation(tag_levels  = "A") &
  theme(plot.margin = margin(mg, mg, mg, mg, unit = "mm"))
supplementary_note_II_figure

# in the end only using figure "A" was sufficient
output_dir = "manuscript/Supplementary_notes/Supplementary_Note_4/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
ggsave("manuscript/Supplementary_notes/Supplementary_Note_4/Supplementary_Note_4_Fig1.png", barplot_ie_vs_sum, width = 12, height = 5, bg = "white")


# average across clones
site_freqs_APC_drivers = cancer_bDM[gene_name == "APC" & driver == TRUE,  .N, by = c("mut_type", "gene_name")]
mrate_APC_drivers = get_gene_rate(exp_rates = expected_rates, metadata = metadata,
                                  site_freqs = site_freqs_APC_drivers, ratios = ratios, ncells = ncells)
plot_APC_drivers_all = plot_mrate(mrate_APC_drivers, title = "number of cells with APC driver mutations", colors = colon_colors)
ggsave("plots/colon/ncells_APC_drivers.png", plot_APC_drivers_all, width = 5.5, height = 4.2, bg = "white")

# double check - the differences between the mutation rates are so small that this indicates rounding errors
#### APC double driver mutations
mrate_APC_double_drivers = get_double_gene_rate(exp_rates = expected_rates, metadata = metadata,
                                                site_freqs = site_freqs_APC_drivers, ratios = ratios, ncells = ncells)
plot_apc_double_drivers_all = plot_mrate(mrate_APC_double_drivers, title = "number of cells with double APC driver mutations - standard method", colors = colon_colors) +
  ggforce::facet_zoom(ylim = c(0, 10), zoom.size = 1)
ggsave("plots/colon/ncells_APC_double_drivers_all_points.png", plot_apc_double_drivers_all,  width = 10, height = 5, bg = "white")

plot_apc_double_drivers_normal = plot_mrate(mrate_APC_double_drivers |> filter(category == "normal"), title = "number of cells with APC driver mutations", colors = colon_colors) +
  labs(subtitle = "normal tissue - using all cells averaged for each donor")
ggsave("plots/colon/ncells_APC_double_drivers_normal_points.png", plot_apc_double_drivers_normal,  width = 6, height = 5, bg = "white")

# taking into account the diversity between different mutational programmes
mrate_APC_double_drivers_fraction = get_double_gene_rate_fraction(exp_rates = expected_rates, metadata = metadata,
                                                                site_freqs = site_freqs_APC_drivers, ratios = ratios, ncells = ncells)
plot_mrate(mrate_APC_double_drivers_fraction, title = "number of cells with double APC driver mutations - fraction method", colors = colon_colors) +
  ggforce::facet_zoom(ylim = c(0, 10), zoom.size = 1)
ggsave("plots/colon/ncells_APC_double_drivers_points_fractions.png", width = 10, height = 5, bg = "white")

# Manuscript: make general overview table of the different conditions. This can be a supplementary or main figure table
mrate_APC_double_drivers_fraction |>
  filter(category == "normal" & age > 35) |>
  pull(mle) |> unique() |> summary()

double_drivers = mrate_APC_double_drivers_fraction |>
  filter(age > 35) |>
  group_by(category) |>
  summarize(`mean exp. \ndouble APC driver` = paste0(round(mean(mle),1), " (", round(min(mle), 1), "-", round(max(mle), 1), ")"),
            `number of donors` = dplyr::n())

single_drivers = mrate_APC_drivers |>  filter(age > 35) |>
  group_by(category) |>
  summarize(`mean age\nin years\n(min-max)` = paste0(round(mean(age),1), " (", round(min(age), 1), "-", round(max(age), 1), ")"),
            `mean cells exp. \nAPC driver (min-max)` = paste0(
              format(round(mean(mle)), big.mark = ","), " (",
              format(round(min(mle)), big.mark = ","), "-",
              format(round(max(mle)), big.mark = ","), ")"))

metadata_counts = metadata |> group_by(category) |>
  summarize(ncrypts = dplyr::n())

table_2 = left_join(single_drivers, double_drivers) |>
  left_join(metadata_counts)
writexl::write_xlsx(table_2, "manuscript/Table_2.xlsx")

# Supplementary Figure 5 manuscript : get the confidence interval of the mean:
mrate_APC_drivers |>  filter(age > 35) |>
  group_by(category) |>
  filter(category == "normal") |>
  summarize(mean = mean(mle),
            low = mean(cilow),
            high = mean(cihigh)) |>
  as.data.frame()

rbindlist(list(`diversity clones` = mrate_APC_double_drivers_fraction |> filter(category == "normal", ),
               `mean across donor` = mrate_APC_double_drivers |> filter(category == "normal", )), idcol = "name") |>
  ggplot(aes(x = fct_reorder(donor, mle), y = mle, color = name)) +
  geom_point(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = cilow, ymax = cihigh), position = position_dodge(width = 0.9)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "inside", legend.position.inside = c(0.02, 0.9)) +
  labs(y = "predicted number of cells with double mutation", x = NULL, color = "computing method",
       subtitle = "Comparison calculation methods # cells with two APC driver mutations")
ggsave("manuscript/Supplementary_notes/Supplementary_Note_4/Supplementary_Note4_Fig2.png", width = 8, height = 6, bg = 'white')
ggsave("manuscript/Supplementary_notes/Supplementary_Note_4/Supplementary_Note4_Fig2.png", width = 8, height = 6, bg = 'white')
