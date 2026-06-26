# APC progression figure — panels B, C, D, E/F
# Run from mutrisk_manuscript root

library(GenomicRanges)
library(patchwork)
library(ggpubr)
source("code/0_functions/analysis_variables.R")
source("code/0_functions/plot_crc_cells_dual_axis.R")

tissue <- "colon"
colors <- tissue_colors[[tissue]]

# ---- 1. UK Biobank CRC incidence ----
crc_freq <- fread("raw_data/UKBiobank/colorectal_cancer_frequency_UKB.csv")

ukbiobank_crc <- data.frame(age = 0:max(crc_freq$current_age)) |>
  mutate(
    n_alive  = sapply(age, \(x) sum(crc_freq$current_age >= x)),
    n_tumor  = sapply(age, \(x) sum(crc_freq$var_Colorectal_age == x, na.rm = TRUE)),
    risk     = n_tumor / n_alive,
    CRC_cumulative_risk = cumsum(risk)
  ) |>
  filter(n_alive > 5000)

ukbiobank_crc <- ukbiobank_crc |>
  mutate(age_group = cut(age, seq(0, 80, 10), labels = seq(5, 75, 10))) |>
  group_by(age_group) |>
  summarise(risk = mean(risk), CRC_cumulative_risk = mean(CRC_cumulative_risk)) |>
  filter(!is.na(age_group)) |>
  dplyr::rename(age = age_group)
ukbiobank_crc$age <- as.numeric(as.character(ukbiobank_crc$age))

# ---- 2. Load data ----
metadata <- fread(paste0("processed_data/", tissue, "/", tissue, "_metadata.tsv")) |>
  distinct() |>
  mutate(category = factor(category, levels = c("normal", "IBD", "POLD1", "POLE")),
         tissue = 'colon')

ncells <- tissue_ncells_ci$mid_estimate[1]
expected_rates <- fread(paste0("processed_data/", tissue, "/", tissue, "_expected_rates.tsv.gz")) |> mutate(tissue = 'colon')
expected_rates_normal <- expected_rates |> filter(category == "normal")
ratios <- fread(paste0("processed_data/", tissue, "/", tissue, "_mut_ratios.tsv.gz")) |>
  filter(gene_name %in% c("APC", "KRAS")) |> mutate(tissue = "colon")

colon_bDM <- fread("processed_data/boostdm/boostdm_genie_cosmic/colon_boostDM_cancer.txt.gz")

# ---- 3. Per-gene per-donor driver rates (following merge_mutrisk_drivers flow) ----
# Helper: mean across samples per donor per mut_type, multiply by ratio
compute_gene_rates <- function(gene, boostdm_table, expected_rates_data, ratios_table, metadata_table) {
  gr <- ratios_table |> filter(gene_name == gene, category == "normal") |> pull(ratio)
  dp <- boostdm_table[gene_name == gene & boostDM_class == TRUE, .(mut_type, position, driver)]

  # Step 1: join metadata, mean across samples per donor per mut_type
  er <- expected_rates_data |>
    left_join(metadata_table, by = c("tissue", "sampleID", "category")) |>
    group_by(donor, mut_type, tissue) |>
    summarise(across(c(mle, cilow, cihigh), mean), .groups = "drop") |>
    mutate(across(c(mle, cilow, cihigh), ~ . * gr))

  # Step 2: keep only driver mut_types, expand to positions, sum per donor
  er <- er |> filter(mut_type %in% unique(dp$mut_type))
  dp |>
    full_join(er, by = "mut_type", relationship = "many-to-many") |>
    group_by(donor) |>
    summarise(across(c(mle, cilow, cihigh), sum), .groups = "drop") |>
    left_join(metadata_table |> select(donor, age, category) |> distinct(), by = "donor")
}

apc_single_snv <- compute_gene_rates("APC", colon_bDM, expected_rates_normal, ratios, metadata)
KRAS_single_snv <- compute_gene_rates("KRAS", colon_bDM, expected_rates_normal, ratios, metadata)

boostdm_list <- list(
  colon = colon_bDM |>
    filter(boostDM_class == TRUE) |>
    mutate(driver = boostDM_class) |>
    select(-boostDM_class) |>
    as.data.table()
  )

# ---- 4. Double-hit rates ----
double_apc <- apc_single_snv |>
  mutate(across(c(mle, cilow, cihigh), ~ ((.^2) / 4)))


double_apc_ncells <- double_apc |>
  mutate(ncells_mut = mle * ncells)

double_apc_kras <- double_apc_ncells |>
  left_join(
    KRAS_single_snv |> rename_with(~ paste0("kras_", .), c(mle, cilow, cihigh)),
    by = c("category", "donor", "age")
  ) |>
  mutate(
    mle        = mle * kras_mle,
    cilow      = cilow * kras_cilow,
    cihigh     = cihigh * kras_cihigh,
    ncells_mut = mle * ncells
  )

# ---- 5. Mutation fractions (from GENIE) ----
mutated_fractions <- fread("processed_data/GENIE_17/CRC_mutation_fractions.txt")
frac_apc_single      <- mutated_fractions$percentages[mutated_fractions$mutated_combination == "APC_single_snv"][1]
frac_apc_double      <- mutated_fractions$percentages[mutated_fractions$mutated_combination == "APC_double"][1]
frac_apc_kras_double <- mutated_fractions$percentages[mutated_fractions$mutated_combination == "KRAS_APC_double"][1]

cat(sprintf("Mutation fractions: APC_single=%.1f%%, APC_double=%.1f%%, APC_KRAS_double=%.1f%%\n",
            frac_apc_single*100, frac_apc_double*100, frac_apc_kras_double*100))

# ---- 6. Dual-axis panel function ----
create_dual_axis_panel <- function(cells_data, crc_data, left_label, right_label,
                                    mutation_fraction = 1, use_scientific = FALSE) {
  crc_df <- crc_data |>
    select(age, CRC_cumulative_risk) |>
    mutate(crc_cases_per_million = CRC_cumulative_risk * mutation_fraction * 1e6)

  cells_df <- cells_data |>
    filter(category == "normal") |>
    mutate(
      cell_count    = mle * ncells,
      cell_count_lo = cilow * ncells,
      cell_count_hi = cihigh * ncells
    )

  crc_max   <- max(crc_df$crc_cases_per_million, na.rm = TRUE)
  cell_max  <- max(cells_df$cell_count, na.rm = TRUE)
  trans_factor <- crc_max / cell_max

  # --- return some numbers:
  # Cases for age 75
  cases_75 <- crc_df |> filter(age == 75) |> pull(crc_cases_per_million)
  cells_75 <- cells_df |> filter(age == 79) |> pull(cell_count)
  print(left_label)
  print(paste0("Age 75: ", signif(cases_75, 3), " CRC cases per 1M, "))

  print(right_label)
  print(paste0("Age 79 donor: ", signif(cells_75, 3), " cells with driver mutation"))

  ggplot() +
    geom_line(data = crc_df, aes(age, crc_cases_per_million),
              color = COLOR_OUTCOME, linewidth = 1.2) +
    geom_point(data = crc_df, aes(age, crc_cases_per_million),
               color = COLOR_OUTCOME, size = 2) +
    geom_pointrange(data = cells_df,
                    aes(age, cell_count * trans_factor,
                        ymin = cell_count_lo * trans_factor,
                        ymax = cell_count_hi * trans_factor),
                    color = COLOR_CELLS, fill = COLOR_CELLS,
                    shape = 21, size = 0.5, stroke = 0.3, alpha = 0.8) +
    scale_y_continuous(
      name = left_label, labels = scales::label_comma(),
      limits = c(0, max(crc_max, max(cells_df$cell_count * trans_factor, na.rm = TRUE)) * 1.1),
      sec.axis = sec_axis(~ . / trans_factor, name = right_label,
                          labels = if (use_scientific) scales::label_scientific()
                                   else scales::label_comma())
    ) +
    scale_x_continuous(limits = c(0, 85), name = "Age (years)") +
    theme_cowplot() +
    theme(
      axis.title.y           = element_text(color = COLOR_OUTCOME, size = 10),
      axis.text.y            = element_text(color = COLOR_OUTCOME, size = 9),
      axis.title.y.right     = element_text(color = COLOR_CELLS,   size = 10),
      axis.text.y.right      = element_text(color = COLOR_CELLS,   size = 9),
      axis.ticks.y           = element_line(color = COLOR_OUTCOME),
      axis.ticks.y.right     = element_line(color = COLOR_CELLS),
      axis.title.x           = element_text(size = 10),
      plot.margin            = margin(10, 10, 10, 10)
    )
}

# ---- 7. Panels A, B, C - raw ggplot objects (prep_plot + composite in compose_figures.R) ----
pA <- create_dual_axis_panel(apc_single_snv, ukbiobank_crc,
  "Colorectal Cancer cases with one APC\ndriver per 1M", "Cells with one APC driver mutation",
  mutation_fraction = frac_apc_single)

pB <- create_dual_axis_panel(double_apc_ncells, ukbiobank_crc,
  "Colorectal Cancer cases with\ndouble APC driver per 1M", "Cells with double APC driver mutations",
  mutation_fraction = frac_apc_double)

pC <- create_dual_axis_panel(double_apc_kras, ukbiobank_crc,
  "Colorectal Cancer cases with\ndouble APC+KRAS driver per 1M", "Cells with double APC+KRAS driver mutations",
  mutation_fraction = frac_apc_kras_double, use_scientific = TRUE)

# ---- 8. Panel D: Expected individuals with double APC+KRAS per 1M ----
double_apc_kras_ind <- double_apc_kras |>
  filter(category == "normal") |>
  mutate(
    cells_per_person = mle * ncells,  # lambda
    per_1M    = mle * ncells * 1e6,
    per_1M_lo = cilow * ncells * 1e6,
    per_1M_hi = cihigh * ncells * 1e6
  )

# CRC line (same fraction as panel D)
crc_apc_kras <- ukbiobank_crc |>
  mutate(crc_per_million = CRC_cumulative_risk * frac_apc_kras_double * 1e6)

# Annotation values
crc_ann <- crc_apc_kras |> filter(age %in% c(55, 65, 75)) |>
  mutate(label = paste0("age ", age, "\n", signif(crc_per_million, 3), " in 1M"))
cell_79 <- double_apc_kras_ind |> filter(age == 79)
cell_79_val  <- cell_79$per_1M

make_pD <- function() {
  y_max <- max(crc_apc_kras$crc_per_million, double_apc_kras_ind$per_1M, na.rm = TRUE)

  ggplot() +
    geom_line(data = crc_apc_kras, aes(age, crc_per_million),
              color = COLOR_OUTCOME, linewidth = 1.2) +
    geom_point(data = crc_apc_kras, aes(age, crc_per_million),
               color = COLOR_OUTCOME, size = 2) +
    geom_pointrange(data = double_apc_kras_ind,
                    aes(age, per_1M, ymin = per_1M_lo, ymax = per_1M_hi),
                    color = COLOR_CELLS, fill = COLOR_CELLS,
                    shape = 21, size = 0.5, stroke = 0.3, alpha = 0.8) +
    geom_text(data = crc_ann, aes(x = age - 10, y = crc_per_million, label = label),
              color = COLOR_OUTCOME, size = 3.5, hjust = 1) +
    geom_segment(data = crc_ann,
                 aes(x = age - 8, xend = age - 1, y = crc_per_million, yend = crc_per_million),
                 color = COLOR_OUTCOME, linewidth = 0.4,
                 arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
    annotate("text", x = 84, y = 600,
             label = paste0("age 79\n", signif(cell_79_val, 3), " in 1M"),
             color = COLOR_CELLS, size = 3.5, hjust = 1) +
    annotate("segment", x = 82, xend = 79, y = 500, yend = cell_79_val + 20,
             color = COLOR_CELLS, linewidth = 0.4,
             arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
    scale_y_continuous(
      name = "Individuals per 1M",
      labels = scales::label_comma(),
      limits = c(0, y_max * 1.1)
    ) +
    scale_x_continuous(limits = c(0, 85), name = "Age (years)") +
    theme_cowplot() +
    theme(
      legend.position        = c(0.02, 0.98),
      legend.justification   = c(0, 1),
      legend.direction       = "vertical",
      legend.background      = element_rect(fill = alpha("white", 0.7), color = NA),
      legend.text            = element_text(size = 6.5),
      legend.key.size        = unit(0.25, "cm"),
      axis.title.y           = element_text(size = 10),
      axis.text.y            = element_text(size = 9),
      axis.title.x           = element_text(size = 10),
      plot.margin            = margin(40, 40, 10, 10)
    )
}

pD <- make_pD()

# ---- 9. Save raw panels as RDS (composed in compose_figures.R) ----
output_dir <- "manuscript/figure_panels/figure_5"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(pA, file.path(output_dir, "figure_5A.rds"))
saveRDS(pB, file.path(output_dir, "figure_5B.rds"))
saveRDS(pC, file.path(output_dir, "figure_5C.rds"))
saveRDS(pD, file.path(output_dir, "figure_5D.rds"))

cat("Done — figure_5 panel RDS files saved\n")
