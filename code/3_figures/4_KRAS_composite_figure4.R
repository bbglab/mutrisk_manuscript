# KRAS composite figure — three-column layout
#   Left:   IntOGen vs MutRisk mirror plots (colon / lung / blood)
#   Middle: G12/G13/G14 residue bars by mutational signature
#   Right:  G12 age-vs-cells scatter plots
#
# Run from the mutrisk_manuscript project root:
#   Rscript code/3_figures/kras_composite_figure.R

library(ggh4x)
source("code/0_functions/analysis_variables.R")

# ---- 1. Load shared data ----

# IntOGen KRAS distribution (from 6_intogen_kras_distribution.R)
intogen = fread("processed_data/intogen/kras_distribution.tsv.gz")

# Metadata per tissue (needed globally for merge_mutrisk_drivers)
metadata_files = c("processed_data/blood/blood_metadata.tsv",
                   "processed_data/colon/colon_metadata.tsv",
                   "processed_data/lung/lung_metadata.tsv")
names(metadata_files) = str_split_i(metadata_files, "\\/", 2)
metadata = lapply(metadata_files, \(x) fread(x)[, .(sampleID, category, age, donor)]) |>
  rbindlist(idcol = "tissue")

# BoostDM driver tables
boostdm_files = list.files("processed_data/boostdm/boostdm_genie_cosmic/",
                           pattern = "lung|colon|CH", full.names = TRUE)
names(boostdm_files) = c("blood", "colon", "lung")
boostdm = lapply(boostdm_files, fread)

# Expected rates and ratios (combined across tissues)
expected_rate_list = list()
ratio_list = list()
for (tissue in c("colon", "blood", "lung")) {
  expected_rate_list[[tissue]] = fread(
    paste0("processed_data/", tissue, "/", tissue, "_expected_rates.tsv.gz"))
  ratio_list[[tissue]] = fread(
    paste0("processed_data/", tissue, "/", tissue, "_mut_ratios.tsv.gz"))
}
expected_rates = rbindlist(expected_rate_list, idcol = "tissue", use.names = TRUE)
ratios = rbindlist(ratio_list, idcol = "tissue", use.names = TRUE)

# ---- 2. Left column — IntOGen / MutRisk mirror plots ----

mirror_colon = make_kras_mirror_plot(intogen, boostdm, expected_rates, ratios,
                                      "colon", "normal")
mirror_lung  = make_kras_mirror_plot(intogen, boostdm, expected_rates, ratios,
                                      "lung",  c("non-smoker", "ex-smoker", "smoker"))
mirror_blood = make_kras_mirror_plot(intogen, boostdm, expected_rates, ratios,
                                      "blood", "normal")

left_column = wrap_plots(mirror_colon, mirror_lung, mirror_blood, ncol = 1,
                         guides = "collect") &
  theme(legend.position = "bottom")
saveRDS(left_column, "manuscript/figure_panels/figure_4/kras_left_column.rds")

# ---- 3. Middle column — G12/G13/G14 signature bars ----

# Load per-tissue signature rates (normalize column names)
sig_colon = fread("processed_data/colon/colon_sig_donor_rates.tsv.gz")
sig_lung  = fread("processed_data/lung/lung_sig_donor_rates.tsv.gz")
sig_blood = fread("processed_data/blood/blood_sig_donor_rates.tsv.gz")

# colon uses "mut_type, mle"; lung/blood use "name, value" — harmonize
setnames(sig_lung,  c("name", "value"), c("mut_type", "mle"))
setnames(sig_blood, c("name", "value"), c("mut_type", "mle"))

# Generate panel plots
bars_colon = make_kras_g12_bars(boostdm, sig_colon, metadata, ratios,
                                 expected_rates, "colon", "normal")
bars_lung_ns = make_kras_g12_bars(boostdm, sig_lung, metadata, ratios,
                                   expected_rates, "lung", "non-smoker")
bars_lung_ex = make_kras_g12_bars(boostdm, sig_lung, metadata, ratios,
                                   expected_rates, "lung", "ex-smoker")
bars_lung_sm = make_kras_g12_bars(boostdm, sig_lung, metadata, ratios,
                                   expected_rates, "lung", "smoker")
bars_blood   = make_kras_g12_bars(boostdm, sig_blood, metadata, ratios,
                                   expected_rates, "blood", "normal")

# ---- Sort AA changes by prominence (lung non-smoker as reference) ----
kras_map = rbindlist(boostdm, idcol = "boostdm_tissue", fill = TRUE) |>
  filter(gene_name == "KRAS", position %in% 12,
         consequence != "synonymous") |>
  select(mut_type, aachange) |>
  distinct()
site_counts = kras_map[, .N, by = "mut_type"]
ns_donors = metadata |>
  filter(tissue == "lung", category == "non-smoker") |>
  pull(donor) |>
  unique()
aa_order = sig_lung |>
  filter(donor %in% ns_donors) |>
  inner_join(site_counts, by = "mut_type", relationship = "many-to-many") |>
  filter(!is.na(signature)) |>
  mutate(mle = mle * N * tissue_ncells$ncells[tissue_ncells$tissue == "lung"]) |>
  left_join(distinct(kras_map, mut_type, aachange), by = "mut_type") |>
  group_by(aachange) |>
  summarize(total = sum(mle), .groups = "drop") |>
  arrange(desc(total)) |>
  pull(aachange)

# Apply global X limits
bars_colon   = bars_colon   + scale_x_discrete(limits = aa_order)
bars_lung_ns = bars_lung_ns + scale_x_discrete(limits = aa_order)
bars_lung_ex = bars_lung_ex + scale_x_discrete(limits = aa_order)
bars_lung_sm = bars_lung_sm + scale_x_discrete(limits = aa_order)
bars_blood   = bars_blood   + scale_x_discrete(limits = aa_order)

# ---- Signature colors matching Figure S2 (dynamic from data) ----
sig_contrib_files = list.files(
  "processed_data/", recursive = TRUE, pattern = "signature_contributions",
  full.names = TRUE)
all_sbs_cols = lapply(sig_contrib_files, function(f) {
  cn = names(fread(f, nrows = 0))
  grep("^SBS", cn, value = TRUE)
}) |> unlist() |> unique()
all_sigs_sorted = all_sbs_cols[order(parse_number(all_sbs_cols), decreasing = FALSE)]
g12_colors = setNames(ggsci::pal_igv()(length(all_sigs_sorted)), all_sigs_sorted)

# Helper: signatures contributing at least min_prop of total mle in a panel
get_active_sigs = function(plot_data, min_prop = 0.01) {
  sig_totals = plot_data |>
    group_by(signature) |>
    summarize(total = sum(mle), .groups = "drop")
  gt = sum(sig_totals$total)
  if (gt == 0) return(character(0))
  sig_totals |> filter(total / gt > min_prop) |> pull(signature) |> as.character()
}

# Only keep signatures that contribute > 1 % in at least one panel
g12_sigs_present = unique(c(
  get_active_sigs(bars_colon$data),
  get_active_sigs(bars_lung_ns$data),
  get_active_sigs(bars_lung_ex$data),
  get_active_sigs(bars_lung_sm$data),
  get_active_sigs(bars_blood$data)
))
active_colors = g12_colors[names(g12_colors) %in% g12_sigs_present]

# Enforce identical fill scale on all panels
shared_fill = scale_fill_manual(values = active_colors, drop = FALSE)
bars_colon   = bars_colon   + shared_fill
bars_lung_ns = bars_lung_ns + shared_fill
bars_lung_ex = bars_lung_ex + shared_fill
bars_lung_sm = bars_lung_sm + shared_fill
bars_blood   = bars_blood   + shared_fill

# Layout styling adjustments
bars_colon   = bars_colon   + theme(axis.title.x = element_blank(),
                                    axis.text.x = element_blank(),
                                    axis.ticks.x = element_blank(),
                                    axis.title.y = element_blank(),
                                    legend.position = "none")
bars_lung_ns = bars_lung_ns + theme(axis.title.x = element_blank(),
                                     axis.text.x = element_blank(),
                                     axis.ticks.x = element_blank(),
                                     axis.title.y = element_blank(),
                                     legend.position = "none")
bars_lung_ex = bars_lung_ex + theme(axis.title.x = element_blank(),
                                     axis.text.x = element_blank(),
                                     axis.ticks.x = element_blank(),
                                     axis.title.y = element_blank(),
                                     legend.position = "none")
bars_lung_sm = bars_lung_sm + theme(axis.title.x = element_blank(),
                                     axis.text.x = element_blank(),
                                     axis.ticks.x = element_blank(),
                                     axis.title.y = element_blank(),
                                     legend.position = "none")
bars_blood   = bars_blood   + theme(axis.title.y = element_blank(),
                                     axis.text.x = element_text(angle = 30, hjust = 1),
                                     legend.position = "none")

# ---- Add total G12 cell count to each panel subtitle ----
fmt_cells = function(x) {
  if (x >= 1) format(round(x), big.mark = ",", scientific = FALSE)
  else if (x >= 0.001) format(round(x, 3), scientific = FALSE)
  else format(x, scientific = TRUE, digits = 2)
}
get_g12_total = function(p) fmt_cells(sum(p$data$mle))
bars_colon   = bars_colon   + labs(subtitle = paste0("colon - normal (", get_g12_total(bars_colon), " cells)"))
bars_lung_ns = bars_lung_ns + labs(subtitle = paste0("lung - non-smoker (", get_g12_total(bars_lung_ns), " cells)"))
bars_lung_ex = bars_lung_ex + labs(subtitle = paste0("lung - ex-smoker (", get_g12_total(bars_lung_ex), " cells)"))
bars_lung_sm = bars_lung_sm + labs(subtitle = paste0("lung - smoker (", get_g12_total(bars_lung_sm), " cells)"))
bars_blood   = bars_blood   + labs(subtitle = paste0("blood - normal (", get_g12_total(bars_blood), " cells)"))

# ---- Combine panels with shared y-axis + standalone legend ----
# Stack the 5 panels vertically
middle_stack = wrap_plots(bars_colon, bars_lung_ns, bars_lung_ex,
                           bars_lung_sm, bars_blood, ncol = 1)

# Standalone legend plot using the same color scale
p_legend = ggplot(data.frame(sig = names(active_colors), y = 1),
                   aes(x = sig, y = y, fill = sig)) +
  geom_col(show.legend = TRUE) +
  scale_fill_manual(values = active_colors, name = "Signature") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_void()
legend_grob = cowplot::get_legend(p_legend)

# Shared y-axis label
y_grob = grid::textGrob("Number of cells with\nKRAS mutation",
                         rot = 90, gp = grid::gpar(fontsize = 11))

middle_column = wrap_plots(y_grob, middle_stack, legend_grob,
                            ncol = 3, widths = c(0.04, 1, 0.2))


# ---- 4. Bottom row — KRAS G12 age vs cells scatter ----

# only use boostdm_class = TRUE
boostdm = lapply(boostdm, function(dt) dt |> dplyr::filter(boostDM_class == TRUE))

tissue_plots_raw = make_kras_g12_scatters(boostdm, expected_rates, ratios, metadata)
saveRDS(tissue_plots_raw, "manuscript/figure_panels/figure_4/kras_g12_scatters.rds")

# ---- 5. Final composite ----

# Row 1: mirror plots (left, 2/3) + G12 SBS bars (right, 1/3)
# Mirror plots share X axis — remove x title from top two
mirror_colon2 = mirror_colon + theme(axis.title.x = element_blank())
mirror_lung2  = mirror_lung  + theme(axis.title.x = element_blank())

mirror_stack = mirror_colon2 + mirror_lung2 + mirror_blood +
  plot_layout(ncol = 1, heights = c(1, 1, 1), guides = "collect") &
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        plot.margin = margin(t = 2, r = 5, b = 2, l = 5)) &
  guides(fill = guide_legend(nrow = 1))

# G12 bars stacked
bars_stack = middle_column + plot_annotation(
  title = "KRAS G12",
  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.05, vjust = 0.5)))

row1_left  = mirror_stack |> prep_plot("A", all_margin = 8)
row1_right = bars_stack   |> prep_plot("B", all_margin = 4)

row1 = row1_left + row1_right + plot_layout(widths = c(2, 1))

# Row 2: age scatter plots side by side
scatter_combined = wrap_plots(tissue_plots_raw, nrow = 1, widths = c(2, 3, 2))

# Add spacing between tissue panels + larger y-axis label
scatter_combined = scatter_combined &
  theme(plot.margin = margin(r = 8, unit = "mm")) &
  theme(axis.title.y = element_text(size = 12))

row2 = scatter_combined |> prep_plot("C", all_margin = 3)

figure = row1 / row2 + plot_layout(heights = c(2.5, 1))

dir.create("manuscript/Figure_4", showWarnings = FALSE, recursive = TRUE)
ggsave("manuscript/Figure_4/kras_composite.png", figure, width = 18, height = 16, dpi = 300)
ggsave("manuscript/Figure_4/kras_composite.svg", figure, width = 18, height = 16)
ggsave("manuscript/Figure_4/kras_composite.pdf", figure, width = 18, height = 16)

cat("Done — saved kras_composite.{png,svg,pdf} to manuscript/Figure_4/\n")
