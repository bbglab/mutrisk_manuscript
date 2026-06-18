# Script focusing on the UKbiobank analyses
library(GenomicRanges)
source("code/0_functions/analysis_variables.R")

# load the bowel cancer data:
tissue = "colon"

## UKBiobank analyses - UKBbiobank data needs to be obtained
crc_freq = fread("raw_data/UKBiobank/colorectal_cancer_frequency_UKB.csv")
colors = tissue_colors[[tissue]]

# calculate the incidence rates
ukbiobank_crc = data.frame(age = 0:max(crc_freq$current_age)) |>
  mutate(n_alive = sapply(age, \(x) sum(crc_freq$current_age >= x)),
         n_tumor = sapply(age, \(x) sum(crc_freq$var_Colorectal_age == x, na.rm = TRUE)),
         n_no_tumor = n_alive - n_tumor,
         risk = n_tumor / n_alive,
         CRC_cumulative_risk = cumsum(risk)) |>
  filter(n_alive > 5000) # filter for the years when at least 5000 patients are in the cohort

# tumor incidence
n_cohort = ggplot(ukbiobank_crc, aes(x = age)) +
  geom_col(aes(y = n_alive)) +
  theme_cowplot() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(y = "number of individuals in cohort", title = "Number of individuals in the cohort") +
  theme(plot.title = element_text(hjust = 0.5))

# tumor incidence
CRC_incidence = ggplot(ukbiobank_crc, aes(x = age)) +
  geom_col(aes(y = n_tumor))  +
  theme_cowplot() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(y = "CRC incidence", title = "CRC Incidence") +
  theme(plot.title = element_text(hjust = 0.5))

# tumor incidence
yearly_incidence_rates = ggplot(ukbiobank_crc, aes(x = age)) +
  geom_col(aes(y = risk)) +
  theme_cowplot() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels = label_percent()) +
  labs(y = "Risk percentage", title = "Yearly incidence rates") +
  theme(plot.title = element_text(hjust = 0.5))

# Cumulative incidence
cumulative_incidence = ggplot(ukbiobank_crc, aes(x = age)) +
  geom_col(aes(y = CRC_cumulative_risk)) +
  theme_cowplot() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels = label_percent()) +
  labs(y = "Risk percentage", title = "Cumulative Incidence") +
  theme(plot.title = element_text(hjust = 0.5))

# combine plots in an overview plot of the UKbiobank incidence rate
n_cohort + CRC_incidence + yearly_incidence_rates + cumulative_incidence

# group the UKBiobank cohort in groups of 10 years
ukbiobank_crc = ukbiobank_crc |>
  mutate(age_group = cut(age, seq(0, 80, 10), labels = seq(5, 75, 10))) |>
  group_by(age_group) |>
  summarize(risk = mean(risk),
            CRC_cumulative_risk = mean(CRC_cumulative_risk)) |>
  filter(!is.na(age_group)) |>
  dplyr::rename(age = age_group)
ukbiobank_crc$age = as.numeric(as.character(ukbiobank_crc$age))

# Compare the UK Biobank data to the globocan incidence
globocan_data = fread("raw_data/globocan_incidence/dataset-cumulative-risk-by-age-in-inc-both-sexes-age-0-84-in-2017-colorectum.csv")
globocan_data = globocan_data |> distinct() |>
  dplyr::rename(country = `Country label`, CRC_cumulative_risk = `Cumulative risk`) |>
  select(country, CRC_cumulative_risk) |>
  mutate(CRC_cumulative_risk = CRC_cumulative_risk / 100)

ukbiobank_crc_cumulative = ukbiobank_crc |>
  filter(age == 79) |>
  mutate(country = "UK - UKBiobank*\n UKB 79 years") |>
  select(country, CRC_cumulative_risk)

rbind(ukbiobank_crc_cumulative, globocan_data) |>
  ggplot(aes(x = country, y = CRC_cumulative_risk)) +
  geom_col() +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(x = NULL, y = "CRC cumulative risk", title = "Cumulative risk across countries",
       subtitle = "CRC cumulative risk at 84 years") +
  geom_text(aes(label = paste0(format(CRC_cumulative_risk*100, digits = 3), "%"), y = CRC_cumulative_risk),
            vjust = -0.2, position = position_dodge(0.9)) +
  theme_cowplot()

# load metadata
metadata = fread(paste0("processed_data/", tissue, "/", tissue, "_metadata.tsv")) |>   distinct() |>
  mutate(category = factor(category, levels = c("normal", "IBD", "POLD1", "POLE")))

# load the mutation rates
ncells = tissue_ncells_ci$mid_estimate[1]
expected_rates = fread(paste0("processed_data/", tissue, "/", tissue, "_expected_rates.tsv.gz"))
expected_rates_normal = expected_rates |>   filter(category == "normal")
ratios = fread(paste0("processed_data/", tissue, "/", tissue, "_mut_ratios.tsv.gz")) |>
  filter(gene_name %in% c("APC", "KRAS", "TP53", "BRAF"))

# load the mutation data:
colon_bDM = fread("processed_data/boostdm/boostdm_genie_cosmic/colon_boostDM_cancer.txt.gz")

# Calculate the driver mutation rates
KRAS_single_snv_muts = colon_bDM[gene_name == "KRAS" & driver == TRUE, .N,
                                  c("gene_name", "mut_type", "aachange", "position")]

KRAS_single_snv = expected_rates_normal |>
  inner_join(KRAS_single_snv_muts, by = "mut_type", relationship = "many-to-many") |>
  left_join(metadata) |>
  left_join(ratios) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ratio)) |>
  group_by(donor, category) |>
  summarise(across(c(mle, cilow, cihigh, age), mean), groups = "drop")

# Manuscript nubmers: compare TP53 R175H: normal min: 62.26459 max: 727.7099 mean: 442.3734
TP53_R = colon_bDM[gene_name == "TP53" & aachange == "R175H", .N, c("gene_name", "mut_type", "aachange", "position")]
expected_rates_normal |>
inner_join(TP53_R, by = "mut_type", relationship = "many-to-many") |>
  left_join(metadata) |>
  left_join(ratios) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ratio * ncells)) |>
  group_by( category, donor) |>
  summarise(mean = mean(mle), min = min(mle), max = max(mle), .groups = "drop_last")  |>
  summarise( min = min(mean), max = max(mean), mean = mean(mean)) |> as.data.frame()

# Manuscript numbers: BRAF V600E mutation:  normal min: 0.3791835 max: 6.559309 mean: 2.256367
BRAF_V600E = colon_bDM[gene_name == "BRAF" & aachange == "V600E", .N, c("gene_name", "mut_type", "aachange", "position")]
expected_rates_normal |>
  inner_join(BRAF_V600E, by = "mut_type", relationship = "many-to-many") |>
  left_join(metadata) |>
  left_join(ratios) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ratio * ncells)) |>
  group_by(category, donor) |>
  summarise(mean = mean(mle), min = min(mle), max = max(mle), .groups = "drop_last")  |>
  summarise( min = min(mean), max = max(mean), mean = mean(mean)) |> as.data.frame()

# Manuscript numbers: Colon TP53 mutations - min: 5847.626 max: 19059 mean; 12782.92
TP53_R = colon_bDM[gene_name == "TP53" & driver == "TRUE", .N, c("gene_name", "mut_type", "aachange", "position")]
TP53_R$N |> sum()
expected_rates_normal |>
  inner_join(TP53_R, by = "mut_type", relationship = "many-to-many") |>
  left_join(metadata) |>
  filter(age > 35) |>
  left_join(ratios) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * N *  ratio * ncells)) |>
  group_by( category, donor, sampleID) |>
  summarise(mle = sum(mle), .groups = "drop_last")  |>
  summarise(mean = mean(mle), min = min(mle), max = max(mle), .groups = "drop_last")  |>
  summarise( min = min(mean), max = max(mean), mean = mean(mean)) |> as.data.frame()

# Expected number of cells with double mutations:
apc_counts_boostdm = colon_bDM[gene_name == "APC" & driver == TRUE, .N, by = c("gene_name", "mut_type",  "driver")]
apc_single_snv = expected_rates_normal |>
  left_join(ratios |> filter(gene_name == "APC")) |>
  left_join(metadata) |>
  inner_join(apc_counts_boostdm |> filter(driver), by = "mut_type", relationship = "many-to-many") |>
  mutate(across(c(mle, cilow, cihigh), ~ . * N * ratio)) |>
  group_by(category, donor, age,  mut_type) |>
  summarize(across(c(mle, cilow, cihigh), mean), .groups = "drop_last") |>
  summarize(across(c(mle, cilow, cihigh), sum), .groups = "drop")

apc_muts_estimated = apc_single_snv |>
  mutate(mle = mle * ncells) |>
  filter(age > 35) |>
  arrange(mle) |>
  as.data.table()
apc_muts_estimated$mle |> max() # 32,098
apc_muts_estimated$mle |> min() # 9332.253 - 9,332 in manuscript
apc_muts_estimated$mle |> mean() # 20934.64 - 20935 in manuscript

# Expected number of cells with double mutations:
apc_muts_estimated = apc_single_snv |>
  mutate(mle = mle * ncells) |>
  arrange(mle) |>  as.data.table()


# for the individuals above 35 years, or of the age 60
apc_muts_estimated_35 = apc_muts_estimated |> filter(age > 35)
apc_muts_estimated$mle |> max()
apc_muts_estimated$mle |> min()
apc_muts_estimated$mle |> mean()

apc_muts_estimated_60 = apc_muts_estimated |> filter(age > 35)
apc_muts_estimated_60$mle |> max()
apc_muts_estimated_60$mle |> min()
apc_muts_estimated_60$mle |> mean()

double_apc = apc_single_snv |>
  mutate(across(c(mle, cilow, cihigh), ~ ((.^2) /4)))
double_apc_ncells = double_apc |>
  mutate(ncells_mut = mle * ncells) |>
  arrange(age) |>   as.data.table()

double_apc_ncells_35 = double_apc_ncells  |> filter(age > 35)
double_apc_ncells$ncells_mut |> min()
double_apc_ncells$ncells_mut |> mean()
double_apc_ncells$ncells_mut |> max()

double_apc_ncells_60 = double_apc_ncells  |> filter(age == 60)
double_apc_ncells_60$ncells_mut |> min() # 0.833 Manuscript: 0.8
double_apc_ncells_60$ncells_mut |> mean() # 2.074525 Manuscript: 2.07
double_apc_ncells_60$ncells_mut |> max() # 3.956599 Manuscript: 3.96

KRAS_single_snv = expected_rates_normal |>
  inner_join(KRAS_single_snv_muts, by = "mut_type", relationship = "many-to-many") |>
  left_join(metadata) |>
  left_join(ratios) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ratio)) |>
  group_by(donor, category) |>
  summarise(across(c(mle, cilow, cihigh, age), mean), groups = "drop")

# Combination of single APC + single KRAS driver mutations
apc_kras = apc_single_snv |>
  left_join(KRAS_single_snv |> select(donor, category, mle_kras = mle, cilow_kras = cilow, cihigh_kras = cihigh),
            by = c("donor", "category")) |>
  mutate(mle = mle * mle_kras,
         cilow = cilow * cilow_kras,
         cihigh = cihigh * cihigh_kras) |>
  select(-mle_kras, -cilow_kras, -cihigh_kras)

double_apc_kras = double_apc_ncells |>
  mutate(mle = mle * KRAS_single_snv$mle,
         cilow = cilow * KRAS_single_snv$cilow,
         cihigh = cihigh * KRAS_single_snv$cihigh)

# Manuscript number across 1-million individuals (average age of the cohort)
mean(double_apc_kras$mle * ncells) * 1e6

# manuscript: number of 60-year olds with 2x APC and KRAS mutation: 1.03 per million individuals
risk_double_apc_kras = double_apc_kras |>
  filter(age == 60) |>
  mutate(mle = mle * ncells) |>
  summarize(mle = mean(mle)) |>
  as.numeric()
risk_double_apc_kras * 1e6 # 1.033262 Manuscript: 1.03

# Plot the driver mutation rates
plot_driver_muts = function(driver_rates, y_axis = "INSERT TITLE") {
  driver_rates |>
    filter(category == "normal") |>
    mutate(across(c(mle, cilow, cihigh), ~ . * ncells)) |>
    ggplot(aes(x = age, y = mle, color = category)) +
    geom_pointrange(aes(ymin = cilow, ymax = cihigh)) +
    scale_y_continuous(labels = scales::label_comma(), limits = c(0, NA)) +
    scale_color_manual(values = colors) +
    theme_cowplot() +
    labs(y = y_axis, x = "Age (years)") +
    theme(legend.position = "none", axis.title = element_text(size = 12))
}

# Make plots not for the manuscript (only KRAS, or single APC + KRAS)
KRAS_single_snv_plot = plot_driver_muts(KRAS_single_snv, y_axis = "Number of cells with\nKRAS driver mutations")
apc_kras_plot = plot_driver_muts(apc_kras, y_axis = "Number of cells with\nAPC + KRAS driver mutations")

# Manuscript-specific plots
F4C = plot_driver_muts(apc_single_snv, y_axis = "Number of cells with\nan APC driver mutation")
F4F = plot_driver_muts(double_apc, y_axis = "Number of cells with\ndouble APC driver mutations")
  F4H = plot_driver_muts(double_apc_kras, y_axis = "Number of cells with\n doubleAPC + KRAS driver mutations") +
  scale_y_continuous(labels = scales::label_number_auto())

# make the figures for
figs = list(F4C = F4C, F4F = F4F, F4H = F4H)

for (fig_i in names(figs)) {
  ggsave(paste0("plots/colon/", fig_i, "_colon.png"), figs[[fig_i]], width = 5, height = 4)
}

saveRDS(figs, "manuscript/figure_panels/figure_4/figures_C_F_H.rds")

# TP53 mutations - check if this needs to be here or in Figure 3 script (the TP53 script)
TP53_single_driver = expected_rates |>
  left_join(ratios |> filter(gene_name == "TP53")) |>
  left_join(metadata) |>
  inner_join(apc_counts_boostdm |> filter(driver), by = "mut_type", relationship = "many-to-many") |>
  mutate(across(c(mle, cilow, cihigh), ~ . * N * ratio)) |>
  group_by(category, donor, age,  mut_type) |>
  summarize(across(c(mle, cilow, cihigh), mean), .groups = "drop_last") |>
  summarize(across(c(mle, cilow, cihigh), sum))

# make this plot for all tissues
TP53_single_driver_plot = TP53_single_driver |>
#  filter(category %in% c("normal", "POLD1")) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ncells),
         category = factor(category, levels = c("normal", "IBD", "POLD1", "POLE"))) |>
  ggplot(aes(x = age, y = mle, color = category)) +
  geom_pointrange(aes(ymin = cilow, ymax = cihigh)) +
  scale_y_continuous(labels = scales::label_comma()) +
  scale_color_manual(values = colors) +
  facet_grid(. ~ category ,axes = "all_y") +
  theme_cowplot() +
  labs(y = 'number of cells with TP53 mutations', x = "Age (years)")
TP53_single_driver_plot

# get the UKBiobank incidence
ukbiobank_crc_plot = ukbiobank_crc |>
  select(age, CRC_cumulative_risk) |>
  pivot_longer(CRC_cumulative_risk, names_to = "name", values_to = "incidence") |>
  mutate(category = "CRC Cumulative Incidence")

# Fraction of CRC mutated for APC or KRAS SNV combinations
mutated_fractions = fread("processed_data/GENIE_17/CRC_mutation_fractions.txt")
writexl::write_xlsx(mutated_fractions, "manuscript/Supplementary_Tables/Supplementary_Table_6.xlsx")

# make combinations of the different ages and CRC-mutation combinations
risk = expand.grid(percentages = mutated_fractions$percentages, CRC_cumulative_risk = ukbiobank_crc$CRC_cumulative_risk)
names = expand.grid(mutation_combination = mutated_fractions$mutated_combination, age = ukbiobank_crc$age)
mutated_fraction_CRC = cbind(names, risk) |> as_tibble() |>
 mutate(CRC_mut_fraction = percentages * CRC_cumulative_risk)

mutated_fraction_CRC |> filter(mutation_combination == "APC_double")
mutated_fraction_CRC |> filter(age == 75)

label_65 = mutated_fraction_CRC |>
  group_by(mutation_combination) |>
  filter(age == 65) |>
  mutate(label = paste0(format(signif(CRC_mut_fraction*100, digits = 2)), "%"))

# single overview plot with all the CRC mutations in there:
ggplot(mutated_fraction_CRC, aes(x = age, y = CRC_mut_fraction)) +
  geom_point() +
  geom_line() +
  ggrepel::geom_text_repel(data = label_65, aes(label = label), nudge_x = -15) +
  ggh4x::facet_wrap2(mutation_combination ~ . , nrow = 1, axes = "y" ) +
  scale_y_continuous(breaks = scales::breaks_pretty(2), labels = label_percent()) +
  scale_x_continuous(limits = c(0, 80)) +
  theme_cowplot() +
  labs(x = "Age (years)", y = "Lifetime risk for CRC with\nspecifc mutation")

# make a new figure splitting all the individual panels into separate ones (so that we can identify the labels)
list_CRC_figures = list()
for (i in unique(mutated_fraction_CRC$mutation_combination)) {
  list_CRC_figures[[i]] = mutated_fraction_CRC |>
    filter(mutation_combination %in% i) |>
    ggplot(aes(x = age, y = CRC_mut_fraction)) +
    geom_point() +
    geom_line() +
    ggrepel::geom_text_repel(data = label_65 |> filter(mutation_combination %in% i), aes(label = label), nudge_x = -15) +
    scale_y_continuous(breaks = scales::breaks_pretty(2), labels = label_percent()) +
    scale_x_continuous(limits = c(0, 80)) +
    theme_cowplot() +
    theme(axis.title = element_text(size = 12)) +
    labs(x = "Age (years)", y = "Lifetime risk for CRC with\nspecifc mutation")
}

saveRDS(list_CRC_figures, "manuscript/figure_panels/figure_4/figure_4_CRC.rds")

# for presentation make figure
plot_CRC_APC = ggplot(mutated_fraction_CRC |> filter(mutation_combination == "APC_double"), aes(x = age, y = CRC_mut_fraction)) +
  geom_point() +
  geom_line() +
  ggrepel::geom_text_repel(data = label_65 |> filter(mutation_combination == "APC_double"), aes(label = label), nudge_x = -15) +
  ggh4x::facet_wrap2(mutation_combination ~ . , nrow = 1, axes = "y" ) +
  scale_y_continuous(breaks = scales::breaks_pretty(2), labels = label_percent()) +
  scale_x_continuous(limits = c(0, 80)) +
  theme_cowplot() +
  labs(x = "Age (years)", y = "Lifetime risk for CRC with\nspecifc mutation")
ggsave("plots/colon/plot_CRC_APC_plot.png", plot_CRC_APC, width = 5, height = 4)

plot_CRC_KRAS_APC = ggplot(mutated_fraction_CRC |> filter(mutation_combination == "KRAS_APC_double"), aes(x = age, y = CRC_mut_fraction)) +
  geom_point() +
  geom_line() +
  ggrepel::geom_text_repel(data = label_65 |> filter(mutation_combination == "KRAS_APC_double"), aes(label = label), nudge_x = -15) +
  ggh4x::facet_wrap2(mutation_combination ~ . , nrow = 1, axes = "y" ) +
  scale_y_continuous(breaks = scales::breaks_pretty(2), labels = label_percent()) +
  scale_x_continuous(limits = c(0, 80)) +
  theme_cowplot() +
  labs(x = "Age (years)", y = "Lifetime risk for CRC with\nspecifc mutation")
ggsave("plots/colon/plot_CRC_KRAS_APC_plot.png", plot_CRC_KRAS_APC, width = 5, height = 4)

#########################################################
# overlay the number of expected polyps in the epithelium
#########################################################
ad_incidence = fread("raw_data/polyp_incidence/Adenomas_by10_intestines.csv")
age_min_column = rep(c(40,50, 60, 70, 80), 2 )
age_max_column = rep(c(50, 60, 70, 80, 90), 2)

# first get the lesions by age, make horizontal bars to add in the text
# Mean Number of Adenomas per Ten Intestines Examined, by Age, Sex, and Adenoma Size
# Multiply by the fraction of adenoma's mutated for APC

# 135 CNADs sequenced, of which 73 have APC mutation
73 / 135  # 54% of adenomas - check where this number is coming from - must be from the paper but

ad_incidence = ad_incidence |>
  filter(!Sex %in% c("Males Total", "Females Total")) |>
  mutate(min_age = age_min_column,
         max_age = age_max_column,
         age = (min_age + max_age) / 2,
         fraction_adenoma_apc = (`All sizes` * 0.54) / 10)

# take the mean across samples
ad_incidence_mean = ad_incidence |>
  group_by(age) |>
  summarize(fraction_adenoma_apc = mean(fraction_adenoma_apc))

AD_plot = ggplot(ad_incidence_mean, aes(x = age, y = fraction_adenoma_apc)) +
  geom_point() +
  geom_line() +
  theme_cowplot() +
  scale_y_continuous(breaks = scales::breaks_extended(5), labels = label_percent(), limits = c(0, .75)) +
  scale_x_continuous(limits = c(20, 85)) +
  labs(y = "Lifetime risk for Adenoma\nwith APC mutation",
       x = "Age (years)")  +
  theme(axis.title = element_text(size = 12))
AD_plot
saveRDS(AD_plot, "manuscript/figure_panels/figure_4/figures_adenoma.rds")

# Supplementary rates: Determine signature-specific effects:
# calculate for the individual donors (colibactin) the number of expected mutations:
# read in the signature-sepecific activity across donors
colon_sig_rates = fread("processed_data/colon/colon_sig_donor_rates.tsv.gz")

apc_counts_boostdm = colon_bDM[gene_name == "APC" & driver == TRUE, .N, by = c("gene_name", "mut_type",  "driver")]

apc_sig_rates = left_join(apc_counts_boostdm, colon_sig_rates) |>
  mutate(mle = mle * N * ncells)
apc_sig_rates |>
  group_by(signature, donor) |>
  summarize(muts = sum(mle)) |>
  left_join(metadata |> select(donor, age, category) |> distinct()) |>
  filter(category == "normal") |>
  ggplot(aes(x = age, y = muts, color = signature)) +
  geom_point() +
  facet_wrap(signature ~. ) +
  theme_cowplot() +
  ggsci::scale_color_igv()

APC_muts_w_o_SBS88 = apc_sig_rates |>
  filter(signature != "SBS88") |>
  group_by( donor) |>
  summarize(`no SBS88` = sum(mle)) |>
  left_join(metadata |> select(donor, age, category) |> distinct()) |>
  filter(category == "normal")


APC_muts_w_o_SBS89 = apc_sig_rates |>
  filter(signature != "SBS89") |>
  group_by( donor) |>
  summarize(`no SBS89` = sum(mle)) |>
  left_join(metadata |> select(donor, age, category) |> distinct()) |>
  filter(category == "normal")

APC_muts = apc_sig_rates |>
  group_by( donor) |>
  summarize(`all signatures` = sum(mle)) |>
  left_join(metadata |> select(donor, age, category) |> distinct()) |>
  filter(category == "normal")

df_APC_noSBS88 = left_join(APC_muts_w_o_SBS88, APC_muts) |>
  pivot_longer(cols = c(`no SBS88`, `all signatures`), values_to = "Cells with APC mutation")

figure_no_SBS88 = df_APC_noSBS88 |>
  ggplot(aes(x = age, y = `Cells with APC mutation`, color = name)) +
  geom_line(aes(group = donor)) +
  geom_point(size = 3) +
  theme_cowplot() +
  scale_color_manual(values = c("black", "darkolivegreen3")) +
  scale_y_continuous(labels = label_comma(), limits = c(0, NA)) +
  labs(x = "Age (years)", color = NULL) +
  theme(legend.position = "inside", legend.position.inside = c(0.1, 0.8))


df_APC_noSBS89 = left_join(APC_muts_w_o_SBS89, APC_muts) |>
  pivot_longer(cols = c(`no SBS89`, `all signatures`), values_to = "Cells with APC mutation")

figure_no_SBS89 = df_APC_noSBS89 |>
  ggplot(aes(x = age, y = `Cells with APC mutation`, color = name)) +
  geom_line(aes(group = donor)) +
  geom_point(size = 3) +
  theme_cowplot() +
  scale_color_manual(values = c("black", "darkgoldenrod3")) +
  scale_y_continuous(labels = label_comma(), limits = c(0, NA)) +
  labs(x = "Age (years)", color = NULL) +
  theme(legend.position = "inside", legend.position.inside = c(0.1, 0.8))

figure_S6A = figure_no_SBS88
saveRDS(figure_S6A, "manuscript/Supplementary_Figures/Figure_S6/figure_S6A.rds")

# summary plot comparing the mutational load
supplementerary_figure_B = prep_plot(figure_no_SBS89, "B")
figure_S6A + supplementerary_figure_B

