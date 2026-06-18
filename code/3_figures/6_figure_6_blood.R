# Figure 5 scripts
library(geomtextpath)
library(ggpubr)
library(gridExtra)
source("code/0_functions/analysis_variables.R")
tissue = "blood"
blood_colors = "#ff725c"

# load blood metadata
ncells = 1e5
tissue = "blood"
metadata = fread(paste0("processed_data/", tissue, "/", tissue, "_metadata.tsv")) |>
  distinct()

# load the mutation rates
expected_rates = fread(paste0("processed_data/", tissue, "/", tissue, "_expected_rates.tsv.gz"))

# load the relative mutation ratio for each sample for the DNMT3A gene:
ratios = fread(paste0("processed_data/", tissue, "/", tissue, "_mut_ratios.tsv.gz"))
sig_donor_rates = fread(paste0("processed_data/", tissue, "/", tissue, "_sig_donor_rates.tsv.gz"))

# load the GENIE data
genie_blood = fread("processed_data/GENIE_17/GENIE_17_processed.txt.gz") |>
  filter(CANCER_TYPE %in% c("Leukemia","B-Lymphoblastic Leukemia/Lymphoma","Myeloproliferative Neoplasms",
                            "Myelodysplastic Syndromes","Mature T and NK Neoplasms",
                            "Myelodysplastic/Myeloproliferative Neoplasms"))

# load boostdm_ch-genie-cosmic intersections
CH_bDM = fread("processed_data/boostdm/boostdm_genie_cosmic/CH_boostDM_cancer.txt.gz")

##### PLOTTING #####
# mutation rates
mutation_rates = expected_rates |>
  mutate(tissue_category = paste0(tissue, "_", category)) |>
  left_join(metadata)

# Mean mutation rate for each trinucleotide
# Check if the blood rates actually make sense - this seems too low - this must be because of the cord blood donors
mean_rates = mutation_rates |>
  filter(age > 0) |>
  group_by(tissue_category, mut_type) |>
  summarize(mle = mean(mle),
            mean_age = mean(age)) |>
  mutate(tissue_category_age = paste0(tissue, "\n(", format(mean_age, digits = 3, nsmall = 1 ), ") years"))

# make specific plots for specific positions:
# Blood make specific mirror genes for specific sites:
boostdm_ch = fread("processed_data/boostdm/boostdm_genie_cosmic/CH_boostDM_cancer.txt.gz")

# UKBiobank DNMT3A mutations:
UKB_DNMT3A_muts = fread("raw_data/UKBiobank/UkBiobank_DNMT3A_mut_age.csv")
UKB_DNMT3A_counts = UKB_DNMT3A_muts[, .N, by = c("aa_change", "REF", "ALT")]  |>
  mutate(position = parse_number(aa_change),
         type = paste0(REF, ">", ALT),
         type = case_match(type, .default = type,
                           "A>T" ~ "T>A", "A>G" ~ "T>C", "A>C" ~ "T>G",
                           "G>T" ~ "C>A", "G>A" ~ "C>T", "G>C" ~ "C>G")) |>
  mutate(mrate = N, tissue_category = "UKBiobank CH") |>
  select(position, type, tissue_category, mrate)

# plot not used in the main script (can be removed)
# # first attempt ukbiobank plots:
# DNMT3A_plot = UKB_DNMT3A_counts |>
#   ggplot(aes(x = position, y = N, fill = type)) +
#   geom_col() +
#   scale_fill_manual(values = COLORS6) +
#   theme_cowplot() +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
#   labs(x  = NULL , y = "Number of CH mutations\nobserved in the UKBiobank cohort", title = "DNMT3A", fill = NULL)


genes_ch = "DNMT3A"
mean_rates_blood = mean_rates |> filter(tissue_category == "blood_normal")

# only plot the driver genes
DNMT3A = boostdm_ch |> filter(gene_name == "DNMT3A")

mutations_blood_DNMT3A = left_join(DNMT3A, mean_rates_blood, relationship = "many-to-many", by = "mut_type") |>
  left_join(triplet_match_substmodel, by = "mut_type") |>
  group_by(position, tissue_category, type) |>
  summarize(mrate = sum(mle) * ncells) |>
  select(position, type, tissue_category, mrate)

library(ggh4x)
df_mirror = bind_rows(mutations_blood_DNMT3A, UKB_DNMT3A_counts) |>
  mutate(
    tissue_category = ifelse(tissue_category == "blood_normal", "Expected mutrate\nblood", tissue_category),
    tissue_category = factor(tissue_category, levels = c("UKBiobank CH", "Expected mutrate\nblood")),
    mrate = ifelse(tissue_category == "Expected mutrate\nblood", 0-mrate, mrate)) |>
  ungroup()

# way to make the plot extend both upper and lower axes
df_point = df_mirror |>
  group_by(tissue_category, position) |>
  summarize(mrate = sum(mrate), .groups = "drop_last") |>
  summarize(mrate = max(abs(mrate)) * 1.1) |>
  mutate(position = 500,
         mrate = ifelse(tissue_category == "Expected mutrate\nblood", 0-mrate, mrate)) |>
  ungroup()

F5A = ggplot(df_point, aes(x = position, y = mrate)) +
  geom_point(color = "white") +
  geom_col(data = df_mirror, aes(fill = type)) +
  geom_text(data = data.frame(tissue_category = factor("UKBiobank CH"), position = 50, mrate = 1500, label = "DNMT3A"),
            aes(label = label)) +
  facet_grid2(tissue_category ~ . , scales = "free") +
  scale_fill_manual(values = COLORS6) +
  theme_cowplot() +
  panel_border() +
  theme(legend.position = "none", panel.spacing.y = unit(0, "mm")) +
  labs(y = "Number expected/\nobserved muts",  x = "AA position") +
  scale_y_continuous(expand=expansion(mult=c(0,0)), breaks = scales::breaks_extended(n = 3), labels = abs)
F5A

saveRDS(F5A, "manuscript/figure_panels/figure_5/figure_5A.rds")


#### Figure 5B-C-D

# calculate expected number of mutated cells
calc_exp_muts = function(expected_rates, mut_positions, metadata, ratios, ncells) {
  expected_rates |>
    inner_join(mut_positions, by = "mut_type", relationship = "many-to-many") |>
    left_join(metadata, by = c("sampleID", "coverage", "category")) |>
    left_join(ratios, by = c("category", "gene_name")) |>
    mutate(across(c(mle, cilow, cihigh), ~ . * ratio * ncells * N)) |>
    group_by(donor, category, sampleID) |>
    summarise(across(c(mle, cilow, cihigh), sum),
              age = mean(age), .groups = "drop_last") |>
    summarise(across(c(mle, cilow, cihigh, age), mean), .groups = "drop")
}

plot_figures = function(driver_sites, y_label, metadata = metadata) {
  driver_rates = calc_exp_muts(expected_rates, driver_sites, metadata= metadata, ratios = ratios, ncells = ncells)

  high_estimate = ggplot(driver_rates,
                         aes(x = age, y = mle*13)) +
    geom_pointrange(aes(ymin = cilow*13, ymax = cihigh*13), color = blood_colors)  +
    labs(y = y_label, x = "Age (years)", subtitle = "1.3 million HSCs") +
    theme_cowplot()


  mid_estimate = ggplot(driver_rates,
                        aes(x = age, y = mle)) +
    geom_pointrange(aes(ymin = cilow, ymax = cihigh), color = blood_colors)  +
    labs(y = y_label, x = "Age (years)", subtitle = "100,000 HSCs") +
    theme_cowplot()

  return(list(high_estimate = high_estimate, mid_estimate = mid_estimate))
}


# DNMT3A driver hotspot
DNMT3A_R882H_hotspot = CH_bDM[aachange == "R882H" & gene_name == "DNMT3A", .N, c("gene_name", "mut_type", "aachange", "position", "boostDM_class")]
DNMT3A_R882H_hotspot_plots = plot_figures(DNMT3A_R882H_hotspot, "Number of cells with\nDNMT3A R882H mutation", metadata = metadata )

figure_S8 = mapply(DNMT3A_R882H_hotspot_plots[2:1], c("A", "B"), FUN = prep_plot) |>
  wrap_plots()
ggsave("manuscript/Supplementary_Figures/Figure_S8/Figure_S8.png", figure_S8, width = 10, height = 5)
ggsave("manuscript/Supplementary_Figures/Figure_S8/Figure_S8.pdf", figure_S8, width = 10, height = 5)


# calculate the time for a clone to be present at a VAF of 0.02:
# assuming 100K HSCs
extension_age = log(26000) / log(1.148)
extension_age = log(2000) / log(1.148)
metadata_extension = metadata |>
  mutate(age = age + extension_age)

driver_rates = calc_exp_muts(expected_rates, DNMT3A_R882H_hotspot, metadata= metadata, ratios = ratios, ncells = ncells)
driver_rates_extension = calc_exp_muts(expected_rates, DNMT3A_R882H_hotspot, metadata= metadata_extension, ratios = ratios, ncells = ncells)


age_shift_figure = list(`mutation induction rate` = driver_rates,
                        `VAF detection limit CH` = driver_rates_extension) |>
  rbindlist(idcol = "type") |>
  ggplot(aes(x = age, y = mle)) +
  geom_pointrange(aes(ymin = cilow, ymax = cihigh, color = type))  +
  labs(x = "Age (years)", subtitle = "100,000 HSCs", y = "Fraction of individuals with ") +
  scale_color_manual(values = c(blood_colors, "#8c463a")) +
  theme_cowplot()
# add in the additional aging shift, make this also part of supplementary figure S7
age_shift_figure

figure_S9A = age_shift_figure

# figure of DNMT3A driver mutations across age
DNMT3A_drivers = CH_bDM[gene_name == "DNMT3A" & boostDM_class == TRUE , .N, c("gene_name", "mut_type", "aachange", "position", "boostDM_class")]
DNMT3A_driver_plot = plot_figures(DNMT3A_drivers, "Number of cells with\n DNMT3A any mutation", metadata)
ggsave("plots/blood/masha_exploration/DNMT3A_driver_plot.png", DNMT3A_driver_plot, width = 5, height = 4.5, bg = "white")


DNMT3A_all_muts = CH_bDM[gene_name == "DNMT3A", .N, c("gene_name", "mut_type", "aachange", "position", "driver")]
DNMT3A_all_plot = plot_figures(DNMT3A_all_muts, "Number of cells with\n DNMT3A any mutation", metadata)
ggsave("plots/blood/masha_exploration/DNMT3A_all_plot.png", DNMT3A_all_plot, width = 5, height = 4.5, bg = "white")


# comparison for Masha:
TET2_drivers = CH_bDM[gene_name == "TET2" & boostDM_class == TRUE , .N, c("gene_name", "mut_type", "aachange", "position", "boostDM_class")]
TET2_driver_plot = plot_figures(TET2_drivers, "Number of cells with\nTET2 driver mutation", metadata)
ggsave("plots/blood/masha_exploration/TET2_driver_plot.png", TET2_driver_plot, width = 5, height = 4.5, bg = "white")

TP53_drivers = CH_bDM[gene_name == "TP53" & boostDM_class == TRUE , .N, c("gene_name", "mut_type", "aachange", "position", "boostDM_class")]
TP53_driver_plot = plot_figures(TP53_drivers, "Number of cells with\nTP53 driver mutation", metadata)
ggsave("plots/blood/masha_exploration/TP53_driver_plot.png", TP53_driver_plot, width = 5, height = 4.5, bg = "white")

# TET2 takes very long and is not used
#TET2_drivers = CH_bDM[gene_name == "TET2", .N, c("gene_name", "mut_type", "aachange", "position", "boostDM_class")]
#TET2_all_plot = plot_figures(TET2_drivers, "Number of cells with\nTET2 any mutation", metadata)
#ggsave("plots/blood/masha_exploration/TET2_all_plot.png", TET2_all_plot, width = 5, height = 4.5, bg = "white")

TP53_drivers = CH_bDM[gene_name == "TP53", .N, c("gene_name", "mut_type", "aachange", "position", "boostDM_class")]
TP53_all_plot = plot_figures(TP53_drivers, "Number of cells with\nTP53 any mutation", metadata)
ggsave("plots/blood/masha_exploration/TP53_all_plot.png", TP53_all_plot, width = 5, height = 4.5, bg = "white")

# Manuscript Numbers
calc_exp_muts(expected_rates, DNMT3A_drivers, metadata= metadata, ratios = ratios, ncells = ncells) |>
  filter(age > 0 ) |>
  summarize(across(c(mle, cilow, cihigh, age), mean))

calc_exp_muts(expected_rates, DNMT3A_drivers, metadata= metadata, ratios = ratios, ncells = ncells) |>
  filter(age > 0 ) |>
  summarize(across(c(mle, cilow, cihigh, age), \(x) mean(x)*13))

# make the general figure:
plots = c(DNMT3A_driver_plot, TET2_driver_plot, TP53_driver_plot)
plots[[2]] = plots[[2]] + labs(title = "DNMT3A") + theme(plot.title = element_text(hjust = 0.5))
plots[[4]] = plots[[4]] + labs(title = "TET2") + theme(plot.title = element_text(hjust = 0.5))
plots[[6]] = plots[[6]] + labs(title = "TP53") + theme(plot.title = element_text(hjust = 0.5))
F5B = wrap_plots(plots[c(2,4,6)], byrow = FALSE) |> prep_plot(label = "B")
F5C = wrap_plots(plots[c(1,3,5)], byrow = FALSE) |> prep_plot(label = "C")

# get all the driver mutations for the watson figure:
watson_variants = c("R882C", "R729W", "R326C", "R320*", "R882H", "R736H",
                    "Y735C", "R736C", "W860R", "R771*", "R598*", "P904L")

DNMT3A_watson_drivers = CH_bDM[gene_name == "DNMT3A" & boostDM_class == TRUE  &
                                 aachange %in% watson_variants, .N, c("gene_name", "mut_type", "aachange", "position", "boostDM_class")]

mutation_list = list(
  DNMT3A_R882H = calc_exp_muts(expected_rates, DNMT3A_R882H_hotspot, metadata, ratios, ncells),
  DNMT3A_drivers = calc_exp_muts(expected_rates, DNMT3A_drivers, metadata, ratios, ncells),
  DNMT3A_watson_drivers = calc_exp_muts(expected_rates, DNMT3A_watson_drivers, metadata, ratios, ncells)) |>
  rbindlist(idcol = "name") |>
  mutate(name = factor(name, levels = c("DNMT3A_drivers", "DNMT3A_R882H", "DNMT3A_watson_drivers")),
         prob_mle = get_prob_mutated_N(risk = mle/ncells ,ncells = ncells, N =  1),
         prob_cilow = get_prob_mutated_N(risk = (mle)/ncells ,ncells = ncells/4, N =  1),
         prob_cihigh = get_prob_mutated_N(risk = (mle)/ncells ,ncells = ncells*4, N =  1),
         prob_mut_5 = get_prob_mutated_N(risk = mle/ncells ,ncells = ncells, N =  5),
         prob_mut_10 = get_prob_mutated_N(risk = mle/ncells ,ncells = ncells, N = 13))
# also needed would be the list of all the variants in the Watson analysis reported to be mutated

# Supplementary Figure: DNMT3A, TET2 and TP53 mutations in CH
# UKbiobank individual counts
UKB_age_frequencies = fread("raw_data/UKBiobank/UKB_age_frequencies_DNMT3A.tsv") |>
  select(Age, Individuals)


# start for loop:
UKB_plot_list = list()
CH_genes = c("DNMT3A", "TET2", "TP53")
for (i in 1:3) {

  gene = CH_genes[i]
  UKB_gene_muts = fread(paste0("raw_data/UKBiobank/UkBiobank_", gene, "_mut_age.csv"))
  UKB_gene_drivers = CH_bDM[gene_name %in% gene & boostDM_class == TRUE]

  age_samples = UKB_gene_muts |>
    filter(aa_change %in% UKB_gene_drivers$aachange) |>
    count(Age)

  UKB_age_incidence = left_join(UKB_age_frequencies, age_samples) |>
    filter(Individuals > 2000) |>
    mutate(relative_incidence = n / Individuals)

  print(gene)
  (UKB_age_incidence |>
    filter(Age == 55)  |> pull(relative_incidence) * 100 )|> print()

  # consider grouping the data as in colon in bins of 10
  plt = ggplot(UKB_age_incidence, aes(x = Age, y = relative_incidence)) +
    geom_point() +
    scale_y_continuous(limits = c(0, NA), labels = label_percent()) +
    theme_cowplot() +
    labs(x = "Age (years)", y = paste0("Incidence of CH in UKB with\n", gene, " driver mutation"), subtitle = gene)

  UKB_plot_list[[gene]] = plt
}
F5D = wrap_plots(UKB_plot_list) |> prep_plot(label = "D")

# save figures:
saveRDS(F5A, "manuscript/figure_panels/figure_6/figure_6A.rds")
saveRDS(F5B, "manuscript/figure_panels/figure_6/figure_6B.rds")
saveRDS(F5C, "manuscript/figure_panels/figure_6/figure_6C.rds")
saveRDS(F5D, "manuscript/figure_panels/figure_6/figure_6D.rds")

# for supplementary figure 9b, add the ukbiobank data
DNMT3A_age = fread("raw_data/UKBiobank/UKB_age_frequencies_DNMT3A.tsv")
DNMT3A_age_fraction = DNMT3A_age |>
  mutate(fraction_DNMT3A_R882H = `R/H` / Individuals) |>
  filter(Individuals >= 2000)

# Manuscript Number
DNMT3A_age_fraction |>
  filter(Age == 55) |>
  pull(`fraction_DNMT3A_R882H`) * 100


figure_S9B = DNMT3A_age_fraction |>
  ggplot(aes(x = Age, y = fraction_DNMT3A_R882H)) +
  geom_pointpath() +
  theme_cowplot() +
  scale_y_continuous(labels = label_percent()) +
  labs(x = "Age (years)", y = "Fraction of UKB individuals\nwith DNMT3A R882H CH")

# Save
figure_S9A = prep_plot(figure_S9A, "A")
figure_S9B = prep_plot(figure_S9B, "B")

figure_S9 = figure_S9A + figure_S9B + plot_layout(widths = c(1.5, 1))
ggsave("manuscript/Supplementary_Figures/Figure_S9/Figure_S9.png", figure_S9, width = 12, height = 5)
ggsave("manuscript/Supplementary_Figures/Figure_S9/Figure_S9.pdf", figure_S9, width = 12, height = 5)

# numbers for supplementary figures:
df_DNMT3A_R882H = calc_exp_muts(expected_rates, DNMT3A_R882H_hotspot, metadata, ratios, ncells)
model_HSCs = lm(mle ~ age, df_DNMT3A_R882H)
summary(model_HSCs)
prediction_DNMT3A_HSCs = predict(model_HSCs, data.frame(age = 60))
coef = coefficients(model_HSCs)
coef[1] + (coef[2] * 60)

# check the effect of setting the intercept at 0
model_HSCs_0 = lm(mle ~ 0 + age , df_DNMT3A_R882H)
predict(model_HSCs_0, data.frame(age = 60))
coef = coefficients(model_HSCs_0)
coef[1] * 60

# model the rate for DNMT3A R882H mutations
DNMT3A_age_fraction
model_CH = lm(fraction_DNMT3A_R882H ~ Age, DNMT3A_age_fraction)
CH_DNMT3A = predict(model_CH, data.frame(Age = 60))
coef = coefficients(model_CH)
coef
coef[1] * 60

# results for manuscript supplementary figure 9:
1/prediction_DNMT3A_HSCs
1/CH_DNMT3A
