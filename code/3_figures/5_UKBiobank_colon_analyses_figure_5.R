# Script focusing on the UKbiobank analyses
library(GenomicRanges)
library(patchwork)
library(ggpubr)  # For prep_plot function
source("code/0_functions/analysis_variables.R")
source("code/0_functions/plot_crc_cells_dual_axis.R")  # For COLOR_OUTCOME, COLOR_CELLS

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
colon_bDM = fread("processed_data/boostdm/boostdm_genie_cosmic/colon_boostDM_cancer.txt.gz") |> mutate(driver = boostDM_class)

# Calculate the driver mutation rates
KRAS_single_snv_muts = colon_bDM[gene_name == "KRAS" & driver == TRUE, .N,
                                  c("gene_name", "mut_type", "aachange", "position")]

KRAS_single_snv = expected_rates_normal |>
  inner_join(KRAS_single_snv_muts, by = "mut_type", relationship = "many-to-many") |>
  left_join(metadata) |>
  left_join(ratios) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * N * ratio)) |>  # Fixed: multiply by N * ratio
  group_by(donor, category) |>
  summarise(across(c(mle, cilow, cihigh, age), mean), groups = "drop")

# Manuscript nubmers: compare TP53 vs BRAF mutation
TP53_R = colon_bDM[gene_name == "TP53" & aachange == "R175H", .N, c("gene_name", "mut_type", "aachange", "position")]
expected_rates_normal |>
inner_join(TP53_R, by = "mut_type", relationship = "many-to-many") |>
  left_join(metadata) |>
  left_join(ratios) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ratio * ncells)) |>
  group_by( category, donor) |>
  summarise(mean = mean(mle), min = min(mle), max = max(mle), .groups = "drop_last")  |>
  summarise( min = min(mean), max = max(mean), mean = mean(mean))

BRAF_V600E = colon_bDM[gene_name == "BRAF" & aachange == "V600E", .N, c("gene_name", "mut_type", "aachange", "position")]
expected_rates_normal |>
  inner_join(BRAF_V600E, by = "mut_type", relationship = "many-to-many") |>
  left_join(metadata) |>
  left_join(ratios) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ratio * ncells)) |>
  group_by(category, donor) |>
  summarise(mean = mean(mle), min = min(mle), max = max(mle), .groups = "drop_last")  |>
  summarise( min = min(mean), max = max(mean), mean = mean(mean))

# Manuscript numbers  ompare TP53 vs APC mutation
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

# Combination of single APC + single KRAS driver mutations
apc_kras = apc_single_snv |>
  left_join(KRAS_single_snv |> select(donor, category, mle_kras = mle, cilow_kras = cilow, cihigh_kras = cihigh),
            by = c("donor", "category")) |>
  mutate(mle = mle * mle_kras,
         cilow = cilow * cilow_kras,
         cihigh = cihigh * cihigh_kras) |>
  select(-mle_kras, -cilow_kras, -cihigh_kras)

apc_muts_estimated = apc_single_snv |>
  mutate(mle = mle * ncells) |>
  filter(age > 35) |>
  arrange(mle) |>
  as.data.table()
apc_muts_estimated$mle |> max()
apc_muts_estimated$mle |> min()
apc_muts_estimated$mle |> mean()

# Expected number of cells with double mutations:
apc_muts_estimated = apc_single_snv |>
  mutate(mle = mle * ncells) |>
  arrange(mle) |>  as.data.table()


# for the individuals above 35 years, or of the age 60
apc_muts_estimated_35 = apc_muts_estimated |> filter(age > 35)
apc_muts_estimated$mle |> max()
apc_muts_estimated$mle |> min()
apc_muts_estimated$mle |> mean()

# Expected number of cells with double mutations: TP53
tp53_counts_boostdm = colon_bDM[gene_name == "TP53" & driver == TRUE, .N, by = c("gene_name", "mut_type",  "driver")]

tp53_single_snv = expected_rates_normal |>
  left_join(ratios |> filter(gene_name == "TP53")) |>
  left_join(metadata) |>
  inner_join(tp53_counts_boostdm |> filter(driver), by = "mut_type", relationship = "many-to-many") |>
  mutate(across(c(mle, cilow, cihigh), ~ . * N * ratio)) |>
  group_by(category, donor, age,  mut_type) |>
  summarize(across(c(mle, cilow, cihigh), mean), .groups = "drop_last") |>
  summarize(across(c(mle, cilow, cihigh), sum), .groups = "drop")

tp53_muts_estimated = tp53_single_snv |>
  mutate(mle = mle * ncells) |>
  filter(age > 35) |>
  arrange(mle) |>
  as.data.table()
tp53_muts_estimated$mle |> max()
tp53_muts_estimated$mle |> min()
tp53_muts_estimated$mle |> mean()


tp53_muts_estimated = tp53_single_snv |>
  mutate(mle = mle * ncells) |>
  filter(age == 60) |>
  arrange(mle) |>
  as.data.table()
tp53_muts_estimated$mle |> max()
tp53_muts_estimated$mle |> min()
tp53_muts_estimated$mle |> mean()

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
double_apc_ncells_60$ncells_mut |> min()
double_apc_ncells_60$ncells_mut |> mean()
double_apc_ncells_60$ncells_mut |> max()

cat("\n=== Double APC only (no KRAS) ===\n")
cat("Unique ages in data:", sort(unique(double_apc_ncells$age)), "\n")
cat("\nCells per person summary by age:\n")
print(double_apc_ncells |>
        group_by(age) |>
        summarize(mean = mean(ncells_mut),
                  min = min(ncells_mut),
                  max = max(ncells_mut),
                  n = n()) |>
        mutate(cells_per_million = mean * 1e6) |>
        arrange(age) |>
        as.data.frame())
cat("\n")

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
  mutate(
    mle = mle * KRAS_single_snv$mle,
    cilow = cilow * KRAS_single_snv$cilow,
    cihigh = cihigh * KRAS_single_snv$cihigh,
    ncells_mut = mle * ncells  # Recalculate cell counts with KRAS included
  )

# Manuscript number
cat("=== Checking double_apc_kras ncells_mut values ===\n")
cat("Sample data (age 60):\n")
print(double_apc_kras |>
        filter(age == 60 & category == "normal") |>
        select(age, donor, mle, ncells_mut) |>
        head(10))
cat("\nMean ncells_mut at age 60:", mean(double_apc_kras$ncells_mut[double_apc_kras$age == 60 & double_apc_kras$category == "normal"]), "\n")
cat("Expected individuals per 1M (mean):", mean(double_apc_kras$ncells_mut[double_apc_kras$age == 60 & double_apc_kras$category == "normal"]) * 1e6, "\n\n")

mean(double_apc_kras$mle * ncells) * 1e6

# manuscript: number of 60-year olds with 2x APC and KRAS mutation: 1.03 per million individuals
risk_double_apc_kras = double_apc_kras |>
  filter(age == 60) |>
  mutate(mle = mle * ncells) |>
  summarize(mle = mean(mle)) |>
  as.numeric()
risk_double_apc_kras * 1e6

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

# ============================================================
# NEW: Create composite 5-panel figure (A-E) like Python script
# ============================================================

# First, compute the CRC mutation fractions needed for Panel E
# Note: This file is generated by code/3_figures/Supplementary_Figures/3_Figure_S7.R
# Make sure to run that script first, or the file should already exist
mutated_fractions = fread("processed_data/GENIE_17/CRC_mutation_fractions.txt")
writexl::write_xlsx(mutated_fractions, "manuscript/Supplementary_Tables/Supplementary_Table_6.xlsx")

# make combinations of the different ages and CRC-mutation combinations
risk = expand.grid(percentages = mutated_fractions$percentages, CRC_cumulative_risk = ukbiobank_crc$CRC_cumulative_risk)
names = expand.grid(mutation_combination = mutated_fractions$mutated_combination, age = ukbiobank_crc$age)
mutated_fraction_CRC = cbind(names, risk) |> as_tibble() |>
 mutate(CRC_mut_fraction = percentages * CRC_cumulative_risk)

label_65 = mutated_fraction_CRC |>
  group_by(mutation_combination) |>
  filter(age == 65) |>
  mutate(label = paste0(format(signif(CRC_mut_fraction*100, digits = 2)), "%"))

# Extract mutation fractions for each panel
# APC mutations can be estimated from epidemiological data (gatekeeper gene)
# KRAS+APC double can also be estimated from sequencing data
frac_apc_single = mutated_fractions$percentages[mutated_fractions$mutated_combination == "APC_single_snv"][1]
frac_apc_double = mutated_fractions$percentages[mutated_fractions$mutated_combination == "APC_double"][1]
frac_apc_kras_double = mutated_fractions$percentages[mutated_fractions$mutated_combination == "KRAS_APC_double"][1]

cat(sprintf("Mutation fractions: APC_single=%.4f (%.1f%%), APC_double=%.4f (%.1f%%), APC_KRAS_double=%.4f (%.1f%%)\n",
            frac_apc_single, frac_apc_single*100, frac_apc_double, frac_apc_double*100,
            frac_apc_kras_double, frac_apc_kras_double*100))

# Function to create dual-axis plots (CRC cases + cell counts)
# mutation_fraction: Fraction of CRC with this specific mutation (from GENIE)
create_dual_axis_panel = function(cells_data, crc_data, left_label, right_label,
                                   mutation_fraction = 1, use_scientific = FALSE) {
  # Prepare CRC data - multiply by mutation fraction to get mutation-specific CRC cases
  crc_df = crc_data |>
    select(age, CRC_cumulative_risk) |>
    mutate(crc_cases_per_million = CRC_cumulative_risk * mutation_fraction * 1e6)

  # Prepare cell data
  cells_df = cells_data |>
    filter(category == "normal") |>
    mutate(
      cell_count = mle * ncells,
      cell_count_lo = cilow * ncells,
      cell_count_hi = cihigh * ncells
    )

  # Calculate the transformation factor
  crc_max = max(crc_df$crc_cases_per_million, na.rm = TRUE)
  cell_max = max(cells_df$cell_count, na.rm = TRUE)
  trans_factor = crc_max / cell_max

  # Create plot with two y-axes
  p = ggplot() +
    geom_line(data = crc_df, aes(x = age, y = crc_cases_per_million),
              color = "#4a4a4a", linewidth = 1.2) +
    geom_point(data = crc_df, aes(x = age, y = crc_cases_per_million),
               color = "#4a4a4a", size = 2) +
    geom_pointrange(data = cells_df,
                    aes(x = age, y = cell_count * trans_factor,
                        ymin = cell_count_lo * trans_factor,
                        ymax = cell_count_hi * trans_factor),
                    color = "#2e8b57", fill = "#2e8b57",
                    shape = 21, size = 0.5, stroke = 0.3, alpha = 0.8) +
    scale_y_continuous(
      name = left_label,
      labels = scales::label_comma(),
      limits = c(0, max(crc_max, max(cells_df$cell_count * trans_factor, na.rm = TRUE)) * 1.1),
      sec.axis = sec_axis(~ . / trans_factor, name = right_label,
                          labels = if (use_scientific) scales::label_scientific()
                                   else scales::label_comma())
    ) +
    scale_x_continuous(limits = c(0, 85)) +
    theme_cowplot() +
    theme(
      axis.title.y = element_text(color = "#4a4a4a", size = 10),
      axis.text.y = element_text(color = "#4a4a4a", size = 9),
      axis.title.y.right = element_text(color = "#2e8b57", size = 10),
      axis.text.y.right = element_text(color = "#2e8b57", size = 9),
      axis.ticks.y = element_line(color = "#4a4a4a"),
      axis.ticks.y.right = element_line(color = "#2e8b57"),
      axis.title.x = element_text(size = 10),
      plot.margin = margin(10, 40, 10, 10)
    ) +
    labs(x = "Age (years)")

  return(p)
}

# Panel B: Single APC + CRC
pB = prep_plot(
  create_dual_axis_panel(
    cells_data = apc_single_snv,
    crc_data = ukbiobank_crc,
    left_label = "CRC cases with one APC driver mutation per 1M",
    right_label = "Cells with one APC driver mutation",
    mutation_fraction = frac_apc_single
  ), "B")

# Panel C: Double APC + CRC
pC = prep_plot(
  create_dual_axis_panel(
    cells_data = double_apc_ncells,
    crc_data = ukbiobank_crc,
    left_label = "CRC cases with double APC driver mutation per 1M",
    right_label = "Cells with double APC driver mutations",
    mutation_fraction = frac_apc_double
  ), "C")

# Panel D: Double APC + KRAS + CRC
pD = prep_plot(
  create_dual_axis_panel(
    cells_data = double_apc_kras,
    crc_data = ukbiobank_crc,
    left_label = "CRC cases with double APC + KRAS driver mutations per 1M",
    right_label = "Cells with double APC + KRAS driver mutations",
    mutation_fraction = frac_apc_kras_double,  # Use KRAS+APC double fraction
    use_scientific = TRUE  # Use scientific notation for small values
  ), "D")

# Panel E: Expected individuals with double APC + KRAS mutations
# Y-axis: Number of individuals out of 1M with at least one double APC + KRAS cell
# Calculation: P(at least 1 cell) = 1 - exp(-lambda) where lambda = cells per person
# Then: individuals_per_1M = P(at least 1) * 1e6
double_apc_kras_for_individuals = double_apc_kras |>
  filter(category == "normal") |>
  mutate(
    lambda = ncells_mut,  # Expected cells per person (from ncells_mut)
    prob_at_least_one = 1 - exp(-lambda),  # Poisson probability
    individuals_per_million = prob_at_least_one * 1e6
  )

# Print detailed step-by-step calculation
cat("\n=== Panel E: Step-by-step calculation for individuals with double APC+KRAS ===\n")
cat("\nSample donor-level data at age 60:\n")
sample_60 = double_apc_kras_for_individuals |>
  filter(age == 60) |>
  select(age, donor, lambda, prob_at_least_one, individuals_per_million) |>
  head(5) |>
  as.data.frame()
print(sample_60)

cat("\nSummary by age (mean across donors):\n")
summary_data = double_apc_kras_for_individuals |>
  group_by(age) |>
  summarize(
    lambda_mean = mean(lambda),
    prob_mean = mean(prob_at_least_one),
    individuals_per_1M = mean(individuals_per_million)
  ) |>
  arrange(age) |>
  as.data.frame()
print(summary_data)

cat("\nKey observation:\n")
cat("  For very small lambda (<<1), P(at least 1) ≈ lambda\n")
cat("  So individuals_per_1M ≈ lambda * 1e6\n")
cat("  Current values at age 60: lambda =", mean(sample_60$lambda),
    ", individuals_per_1M =", mean(sample_60$individuals_per_million), "\n\n")

# Print data for Panel D (same underlying data)
cat("=== Panel D: Cells with Double APC + KRAS Driver Mutations ===\n")
cat("Age\t\tDonor\t\tCells per person (mle * ncells)\n")
panel_d_data = double_apc_kras |>
  filter(category == "normal") |>
  mutate(cell_count_per_person = mle * ncells) |>
  select(age, donor, cell_count_per_person)
print(panel_d_data |> as.data.frame() |> format(digits = 3), row.names = FALSE)
cat("\n")

# Create CRC data for panel E (same mutation fraction as panel D)
crc_double_apc_kras = ukbiobank_crc |>
  select(age, CRC_cumulative_risk) |>
  mutate(
    crc_cases_per_million = CRC_cumulative_risk * frac_apc_kras_double * 1e6
  )

# Create panel E plot
# Single y-axis: Individuals with double APC + KRAS per 1M
# Line = CRC cases with double APC + KRAS, Scatter = Expected carriers (individuals with ≥1 cell)
pE = ggplot() +
  # CRC cases (line)
  geom_line(data = crc_double_apc_kras,
            aes(x = age, y = crc_cases_per_million),
            color = COLOR_OUTCOME, linewidth = 1.2) +
  geom_point(data = crc_double_apc_kras,
             aes(x = age, y = crc_cases_per_million),
             color = COLOR_OUTCOME, size = 2) +
  # Expected individuals (scatter) - using Poisson probability
  geom_point(data = double_apc_kras_for_individuals,
             aes(x = age, y = individuals_per_million),
             color = COLOR_CELLS, size = 2.5, alpha = 0.8) +
  scale_y_continuous(
    name = "Individuals / CRC cases per 1M",
    labels = scales::label_comma()
  ) +
  scale_x_continuous(limits = c(0, 85), name = "Age (years)") +
  theme_cowplot() +
  theme(
    axis.title.y = element_text(color = COLOR_OUTCOME, size = 10),
    axis.text.y = element_text(color = COLOR_OUTCOME, size = 9),
    axis.ticks.y = element_line(color = COLOR_OUTCOME),
    plot.margin = margin(10, 30, 10, 10)
  )

pE = prep_plot(pE, "E")

# Create composite figure: B / (C D) / E (Panel A in Figure S7)
composite_figure = pB / (pC + pD) / pE +
  plot_layout(heights = c(0.8, 1, 1.1), widths = c(1, 1)) +
  plot_annotation(
    title = "APC/KRAS driver mutation burden and colorectal cancer incidence",
    theme = theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0),
      plot.title.padding = margin(0, 0, 10, 0)
    )
  )

# Save composite figure
output_dir = "manuscript/figure_panels/figure_4"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(output_dir, "figure_apc_progression_BCDE.png"),
       composite_figure, width = 10, height = 11, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "figure_apc_progression_BCDE.svg"),
       composite_figure, width = 10, height = 11, bg = "white")

# Save individual panels for flexibility
panel_plots = list(B = pB, C = pC, D = pD, E = pE)
saveRDS(panel_plots, file.path(output_dir, "figures_BCDE.rds"))

print(paste("Composite figure saved to:", file.path(output_dir, "figure_apc_progression_BCDE.png")))

# End of composite figure section
# ============================================================

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

# The CRC mutation fractions are now computed earlier (before composite figure)
# This block is kept for the individual panel plots that follow

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

figure_S5A = figure_no_SBS88
dir.create("manuscript/figure_panels/figure_s5", showWarnings = FALSE, recursive = TRUE)
saveRDS(figure_S5A, "manuscript/figure_panels/figure_s5/Figure_S5A.rds")

# summary plot comparing the mutational load
supplementerary_figure_B = prep_plot(figure_no_SBS89, "B")
figure_S5A + supplementerary_figure_B
