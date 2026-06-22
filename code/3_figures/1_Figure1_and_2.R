### Exome-wide analyses:
# Aims of the script: Generate plots for figure 1 and 2 of the manuscript
# 1. get the exome-wide mutation rate estimates across the different cohort
# 2. get the exome-wide mutation "saturation" plots across the genome

library(ggh4x)
source("code/0_functions/analysis_variables.R")
set.seed(1234)

# load metadata
md_files = list.files("processed_data/", recursive = TRUE, pattern = "_metadata",
                      full.names = TRUE)
names(md_files) = gsub("_.*", "", basename(md_files))
metadata = lapply(md_files, fread) |>
  rbindlist(idcol = "tissue", use.names = TRUE, fill = TRUE) |>
  dplyr::select(any_of(c("tissue", "sampleID", "category", "age", "donor"))) |>
  dplyr::distinct()

metadata |> select(tissue, age, donor) |> distinct() |>
    filter(donor %in% c("O340", "PD34215", "KX008"))

save_plots = function(plot_list, path, name, width = 13, height = 8) {

  for (plot_name in names(plot_list))
    ggsave(paste0(path, "/", name, "_", plot_name, ".png"),
           plot_list[[plot_name]], width = width, height = height, bg = "white")
}


plot_intersect_boxplot = function(intersects, analysis_name, groupby) {

  intersects = intersects |>
    mutate(color = paste0(tissue, "_", category))
  ggplot(intersects, aes(x = category, y = prob, fill = color)) +
    geom_boxplot(color = "grey40", outliers = FALSE) +
    geom_jitter(width = 0.2, alpha = 0.5, shape = 21, color = "white", stroke = 0.2, size = 3) +
    facet_grid(~ tissue, scales = "free_x", space = "free") +
    cowplot::theme_cowplot() +
    scale_fill_manual(values = tissue_category_colors) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(x = NULL, y = "Proportion of SNVs\nin at least one cell")
}

plot_intersect_boxplot_ci = function(intersects_ci, analysis_name, groupby) {

  intersects_ci = intersects_ci |>
    select(-ncells) |>
    pivot_wider(values_from = prob, names_from = "name")

  intersects = intersects_ci |>
    mutate(color = paste0(tissue, "_", category))
  ggplot(intersects, aes(x = category, y = mid_estimate, fill = color)) +
    geom_boxplot(color = "grey40", outliers = FALSE) +
    geom_pointrange(aes(ymin = low_estimate, ymax = high_estimate), alpha = 0.6, position = "jitter") +
    facet_grid(~ tissue, scales = "free_x", space = "free") +
    cowplot::theme_cowplot() +
    scale_fill_manual(values = tissue_category_colors) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(x = NULL, y = "Proportion of SNVs\nin at least one cell")
}


plot_prob_curve = function(probability_rates, analysis_name, groupby, nrow = 2) {
  rates = probability_rates |>
    mutate(color_id = paste0(tissue, "_", category))
  tissue_ncells_plot = rates |>
    select(tissue, category, color_id) |> distinct() |>
    left_join(tissue_ncells)

  ggplot(rates, aes(x = ncells, y = prob, group = groupID, color = color_id)) +
    geom_line() +
    geom_vline(data = tissue_ncells_plot, aes(xintercept = ncells), linetype = "dashed") +
    facet_nested_wrap(. ~ tissue + category, nrow = nrow, nest_line = element_line(linetype = 1),
                      axes = 'x') +
    cowplot::theme_cowplot() +
    cowplot::panel_border() +
    scale_color_manual(values = tissue_category_colors) +
    labs(x = "Number of cells in tissue",
         y = "Fraction of all possible mutations\npresent in at least one cell",
         title = analysis_name, subtitle = paste("grouped by:", groupby)) +
    theme(legend.position = "none",
          panel.spacing = unit(0.2,"lines"),
          strip.background = element_rect(color = "white", fill = "white"),
          panel.border = element_rect(color = "black")) +
    scale_x_log10(guide = "axis_logticks", labels = scales::label_log())
}

plot_prob_curve_ci = function(probability_rates, analysis_name, groupby, ncell_df, nrow = 2) {
  rates = probability_rates |>
    mutate(color_id = paste0(tissue, "_", category))
  tissue_ncells_plot = rates |>
    select(tissue, category, color_id) |> distinct() |>
    left_join(ncell_df) |>
    mutate(label = format_bignum(mid_estimate))

  blood_df = data.frame(tissue = factor("blood", levels = levels(rates$tissue)),
                        x_coord = as.numeric(ncell_df[c("high_estimate", "low_estimate")][3,]),
                        label = c("1.3 million", "25,000"))

  ggplot(rates) +
    geom_rect(data = tissue_ncells_plot,
              aes(xmin = low_estimate, xmax = high_estimate, ymin = -Inf, ymax = 1),
              color = "white", alpha = 0.5) +
    geom_line(aes(x = ncells, y = prob, group = groupID, color = color_id)) +
    geom_segment(data = tissue_ncells_plot, aes(x = mid_estimate, xend =  mid_estimate, y = 0, yend = 1),
                 linetype = "dashed", color = "black") +
    geomtextpath::geom_textline(data = tissue_ncells_plot, aes(x = mid_estimate, label = label, y = 0.5), hjust = 0.5, vjust = 1.2, angle = 90,
                                color = c("black", "black", "white"), size = 5) +
    geom_text(data = blood_df, aes(x = x_coord, y = 0.5, label = label), color = "black", angle = 90, vjust = c(1.1, -.1)) +
    facet_nested_wrap(. ~ tissue,nrow = nrow, nest_line = element_line(linetype = 1), axes = 'all', remove_labels = "y",
                      labeller = labeller(tissue = c(colon = "Colon", lung = "Lung", blood = "Blood"))) +
    cowplot::theme_cowplot() +
    scale_color_manual(values = tissue_category_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, NA)), breaks = extended_breaks(3)) +
    scale_x_log10(guide = "axis_logticks", labels = scales::label_log(), expand = expansion(mult = 0.1, 0.15)) +
    labs(x = "Number of cells in tissue",
         y = "Fraction of all possible mutations\npresent in at least one cell") +
    theme(legend.position = "none",
          panel.spacing = unit(0.2,"lines"),
          strip.background = element_blank()) +
    coord_cartesian(clip = "off")
}

# saturation vs age functions:
# preprocessing function for the saturation function
get_saturation_ci = function(intersects, groupby) {
  meta_donor = metadata |> select(contains(c(groupby, "age"))) |> distinct()
  intersects[[groupby]] = intersects$groupID

  intersects |>
    left_join(meta_donor) |>
    select(-ncells) |>
    pivot_wider(names_from = "name", values_from = "prob") |>
    mutate(tissue_category = paste0(tissue, "_", category))
}

get_saturation = function(intersects, groupby) {
  meta_donor = metadata |> select(contains(c(groupby, "age"))) |> distinct()
  intersects[[groupby]] = intersects$groupID
  intersects |>
    left_join(meta_donor) |>
    mutate(tissue_category = paste0(tissue, "_", category))
}

# function to code the age vs. saturation as a dotplot
plot_saturation_age = function(intersects_ci, analysis_name, groupby) {
  get_saturation(intersects_ci, groupby) |>
    select(-ncells) |>
    pivot_wider(values_from = "prob") |>
    ggplot(aes(x = age, y = mid_estimate, fill = tissue_category)) +
    geom_point(shape = 21, color = "white", stroke = 0.2, size = 3) +
    cowplot::theme_cowplot() +
    scale_fill_manual(values = tissue_category_colors) +
    theme(axis.text.x = element_text(vjust = 1, hjust = 1)) +
    labs(x = NULL,
         y = "Fraction of all possible mutations\npresent in at least one cell")
}

# function to code the age vs. saturation as a pointrange plot with confidence interval of the mutation rate
plot_saturation_age_ci = function(intersects_ci, analysis_name, groupby)  {
  get_saturation_ci(intersects_ci, groupby) |>
    ggplot(aes(x = age, fill = tissue_category)) +
    geom_pointrange(aes(ymin = low_estimate, y = mid_estimate, ymax = high_estimate), shape = 21, color = "black", stroke = 0.5) +
    facet_grid(~ tissue, scales = "free_x", space = "free") +
    cowplot::theme_cowplot() +
    scale_fill_manual(values = tissue_category_colors) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(x = NULL, y = "Proportion of SNVs\nin at least one cell",
         title = analysis_name, subtitle = paste("grouped by:", groupby))
}

analyze_probability = function(gene_counts, analysis_name, groupby = "donor",
                               filter_normal = FALSE) {

  prob_rates_tissue = list()
  prob_intersect_tissue = list()
  prob_intersects_ci_tissue = list()
  result_plot_list = list()

  for (tissue in c("colon", "lung", "blood")) {
    print(tissue)

    expected_rates = fread(paste0("processed_data/", tissue, "/", tissue, "_expected_rates.tsv.gz"))
    if (filter_normal == TRUE) {
      expected_rates = expected_rates |> filter(category %in% c("normal", "non-smoker"))
    }

    expected_rates = expected_rates |>
      left_join(metadata, by = c("sampleID", "category")) |>
      as.data.frame()

    expected_rates = expected_rates |>
      mutate(groupid = expected_rates[[groupby]]) |>
      group_by(category, mut_type, groupid, age) |>
      summarize(across(c(mle, cilow, cihigh), mean), .groups = "drop") |>
      setDT()

    ratios = fread(paste0("processed_data/", tissue, "/", tissue, "_mut_ratios.tsv.gz"))

    category_list = list()
    for (category_select in unique(expected_rates$category)) {

      print(category_select)

      # calculate the risk for specific categories:
      ratios_cat = ratios[category == category_select,]
      expected_rates_cat = expected_rates[category == category_select, ]
      expected_rates_list = split(expected_rates_cat, expected_rates_cat$groupid)

      m = matrix(NA, ncol = length(expected_rates_list), nrow = nrow(gene_counts))
      colnames(m) = names(expected_rates_list)

      for (groupid in names(expected_rates_list)) {
        # print(groupid)
        expected_rates_sample = expected_rates_list[[groupid]]
        genome_rates = get_adjusted_rates(expected_rates_sample = expected_rates_sample,
                                          gene_counts = gene_counts, ratios_cat = ratios_cat)
        m[,groupid] = genome_rates$adjusted_rates
      }

      m_sorted = apply(m, 2, sort)

      if (class(m_sorted)[1] != "matrix") {
        m_sorted = t(m_sorted) # prevent apply converting the matrix to a vector
      }

      # Compute means per bin if the dataframe is bigger or equal to 1e5
      if (nrow(m_sorted) >= 1e5) {
        binsize = 1000
        group <- rep(1:(nrow(m_sorted) %/% binsize), each = binsize)
        output_m <- rowsum(m_sorted, group) / binsize

      } else { output_m = m_sorted}
      colnames(output_m) = colnames(m)
      category_list[[category_select]] = as.data.table(output_m)
    }

    plot_data = lapply(category_list, pivot_longer, everything()) |>
      rbindlist(idcol = "category") |> arrange(value)

    age_df = metadata |> select(contains(groupby), "age") |>
      `colnames<-`(c("name", "age")) |> distinct()

    result_plot = plot_data |>
      left_join(age_df, by = 'name') |>
      arrange(category, name, value) |>
      group_by(category, name)

    binsize = unique(count(result_plot)$n)

    # get the cell counts
    ncells = tissue_ncells$ncells[tissue_ncells$tissue == tissue]
    ncells_ci = tissue_ncells_ci[tissue_ncells_ci$tissue == tissue,]
    ncells_ci_wide = tissue_ncells_ci_wide[tissue_ncells$tissue == tissue,]

    # This code was needed to be added back from f29e078 where it was removed
    result_plot_df = result_plot |>
      ungroup() |>
      mutate(n_mutated_cells = value * ncells,
             ncells = ncells,
             x = rep(1:binsize, n_groups(result_plot)) / 1e4,
             x = x/max(x)) # normalize to 1

    result_plot_list[[tissue]] = result_plot_df # save for later analysis

    ##### get the mutated probability
    gene_site_counts = gene_counts

    # plot the overall rates
    prob_rates = list()
    prob_intersects = list()
    prob_intersects_ci = list()
    for (category_select in unique(expected_rates$category)) {

      print(category_select)
      # calculate the risk for specific categories:
      ratios_cat = ratios[category == category_select,]
      expected_rates_cat = expected_rates[category == category_select, ]
      expected_rates_list = split(expected_rates_cat, expected_rates_cat$groupid)

      list_probs  = list()
      list_intersects = list()
      list_intersects_ci = list()
      for (name in names(expected_rates_list)) {
        # print(name)
        list_probs[[name]] = get_prob_mutated_range(expected_rates_sample = expected_rates_list[[name]],
                                                    gene_counts = gene_site_counts,
                                                    ratios_cat = ratios_cat,
                                                    step = 0.2)

        # get the rates for the tissue-specific number of cells
        list_intersects[[name]] = get_prob_mutated_range(expected_rates_sample = expected_rates_list[[name]],
                                                         gene_counts = gene_site_counts,
                                                         ratios_cat = ratios_cat, range = ncells)

        list_intersects_ci[[name]] = get_prob_mutated_range(expected_rates_sample = expected_rates_list[[name]],
                                                            gene_counts = gene_site_counts,
                                                            ratios_cat = ratios_cat, range = as.numeric(ncells_ci_wide[-1])) |>
          mutate(name = ncells_ci_wide[-1] |> names())
      }

      prob_rates[[category_select]] = rbindlist(list_probs, idcol = "groupID")
      prob_intersects[[category_select]] = rbindlist(list_intersects, idcol = "groupID")
      prob_intersects_ci[[category_select]] = rbindlist(list_intersects_ci, idcol = "groupID")

    }

    prob_rates_tissue[[tissue]] = rbindlist(prob_rates, idcol = "category")
    prob_intersect_tissue[[tissue]] = rbindlist(prob_intersects, idcol = "category")
    prob_intersects_ci_tissue[[tissue]] = rbindlist(prob_intersects_ci, idcol = "category")
  }

  # rbind, set tissue as factor - to save lines
  bind_dfs = function(df_list) {
    rbindlist(df_list, idcol = "tissue") |>
      mutate(tissue = factor(tissue, levels = c("colon", "lung", "blood")))
  }

  intersects = bind_dfs(prob_intersect_tissue)
  intersects_ci = bind_dfs(prob_intersects_ci_tissue)
  probability_rates = bind_dfs(prob_rates_tissue)
  plot_results = bind_dfs(result_plot_list)

  plot_list = list()
  plot_list[["plot_saturation"]] = plot_intersect_boxplot(intersects, analysis_name, groupby)
  plot_list[["plot_saturation_ci"]] = plot_intersect_boxplot_ci(intersects_ci, analysis_name, groupby)

  if (filter_normal == TRUE) {
    nrow = 1
    plot_list[["plot_saturation_curve_ci"]] = plot_prob_curve_ci(probability_rates = probability_rates, analysis_name, groupby, ncell_df = tissue_ncells_ci, nrow = nrow)

    metadata_donor = metadata |> select(donor, sampleID) |> distinct()
    metadata_donor$groupID = metadata_donor[[groupby]]
    probability_rates_select = left_join(probability_rates, metadata_donor) |>
      filter(donor %in% c("O340", "PD34215", "KX008"))
    plot_list[["plot_saturation_curve_ci_single_ncells"]] = plot_prob_curve_ci(probability_rates = probability_rates_select,
                                                                        analysis_name, groupby, ncell_df = tissue_ncells_ci, nrow = nrow)

    plot_list[["plot_saturation_curve_ci_single_ncells_line"]] = plot_prob_curve_ci(probability_rates = probability_rates_select,
                                                                                  analysis_name, groupby, ncell_df = tissue_ncells_ci |> filter(tissue == "none"), nrow = nrow)



  } else { nrow = 2}

  plot_list[["plot_saturation_curve"]] = plot_prob_curve(probability_rates = probability_rates, analysis_name, groupby,nrow = nrow)
  #plot_list[["plot_saturation_curve_ci_wide"]] = plot_prob_curve_ci(probability_rates = probability_rates, analysis_name, groupby, ncell_df = tissue_ncells_ci_wide, nrow = nrow)
  plot_list[["plot_probabilites_normal_colon"]] = plot_prob_curve(probability_rates |> filter(tissue == "colon" & category == "normal"), analysis_name, groupby, nrow = 1) +
    labs(subtitle = NULL)


  table(intersects_ci$tissue)/3


  plot_list[["plot_saturation_age"]] = plot_saturation_age(intersects_ci |> filter(name == "mid_estimate"), analysis_name, groupby)
  plot_list[["plot_saturation_age_ci"]] = plot_saturation_age_ci(intersects_ci, analysis_name, groupby)

  return(list(intersects = intersects, rates = probability_rates, result_plot_df = plot_results,  plot_list = plot_list))
}

# Exome: get the individual counts across the entire exome
genes = lapply(wintr::RefCDS_WGS, \(x) as.data.table(x[["L"]][,1:4]))
names(genes) = sapply(wintr::RefCDS_WGS, \(x) x[["gene_name"]])
gene_counts = rbindlist(genes, idcol = "gene_name")
gene_counts$mut_type = rep(mutrisk::triplet_match_substmodel$mut_type, nrow(gene_counts)/192)
colnames(gene_counts)[-c(1,6)] = c("syn", "mis", "nonsense", "splice")
gene_counts = gene_counts |>
  pivot_longer(-c(gene_name, mut_type), names_to = "consequence", values_to = "count") |>
  filter(count != 0) |>
  setDT()
gene_counts = sample_n(gene_counts, 1e5)

# analyze the frequency of mutations using the groupby setting to "donor". Used for the main figures
exome_analysis_normal = analyze_probability(gene_counts = gene_counts, analysis_name = "exome analysis", groupby = "donor", filter_normal = TRUE)
save_plots(exome_analysis_normal$plot_list, path = "plots/coverage_saturation/", name = "normal_exome", width = 7, height = 5)
save_plots(exome_analysis_normal$plot_list, path = "plots/coverage_saturation/", name = "normal_exome_wideplot", width = 10, height = 5)
exome_analysis = analyze_probability(gene_counts = gene_counts, analysis_name = "exome analysis", groupby = "donor")

# save for presentation:
ggsave("plots/coverage_saturation/presentation_line.png", exome_analysis_normal$plot_list$plot_saturation_curve_ci_single_ncells,
       width = 10, height = 3.7)
ggsave("plots/coverage_saturation/presentation_line1.png", exome_analysis_normal$plot_list$plot_saturation_curve_ci_single_ncells_line,
       width = 10, height = 3.7)
ggsave("plots/coverage_saturation/presentation_line2.png", exome_analysis_normal$plot_list$plot_saturation_curve_ci,
       width = 10, height = 3.7)

# perform similar analysis, now setting the groupby to "sampleID". Used for supplementary figures
#sampleID_exome_analysis = analyze_probability(gene_counts = gene_counts, analysis_name = "exome analysis", groupby = "sampleID")
sampleID_exome_analysis_normal = analyze_probability(gene_counts = gene_counts, analysis_name = "exome analysis", groupby = "sampleID", filter_normal = TRUE)
save_plots(sampleID_exome_analysis_normal$plot_list, path = "plots/coverage_saturation/", name = "sampleID_normal_exome", width = 7, height = 5)
save_plots(sampleID_exome_analysis_normal$plot_list, path = "plots/coverage_saturation/", name = "sampleID_normal_exome_wideplot", width = 10, height = 5)

figure_S4A = sampleID_exome_analysis_normal$plot_list$plot_saturation_curve_ci |>
  prep_plot("A", all_margin = 8)
figure_S4B = sampleID_exome_analysis_normal$plot_list$plot_saturation_age +
  facet_wrap(. ~ tissue, scales = "free_y") +
  cowplot::panel_border() +
  labs(x = "Age (years)") +
  theme(legend.position = "none")
figure_S4B = prep_plot(figure_S4B, "B", all_margin = 8)
figure_S4 = figure_S4A / figure_S4B + plot_layout(widths = c(1.65, 1))
ggsave("manuscript/Supplementary_Figures/Figure_S4/Figure_S4.png", figure_S4, width = 12, height = 8)

# Save as supplementary plot:
exome_analysis_normal$plot_list$plot_saturation_curve_ci
figure_2D = exome_analysis_normal$plot_list$plot_saturation_age +
  theme(legend.position = "inside", legend.position.inside = c(0.6, 0.5)) +
  labs(fill = NULL) +
  scale_fill_manual(labels = c("Blood normal", "Colon normal", "Lung non-smoker"), values = tissue_category_colors)

exome_analysis_normal$plot_list$plot_saturation_age +
  geom_smooth(aes(color = tissue_category), formula = y ~ x,se = FALSE,
              method = "glm", fullrange = TRUE,
              method.args = list(family = quasibinomial(link = "probit")), show.legend = FALSE) +
  theme(legend.position = "inside", legend.position.inside = c(0.6, 0.5)) +
  labs(fill = NULL) +
  scale_fill_manual(labels = c("Blood normal", "Colon normal", "Lung non-smoker"), values = tissue_category_colors)

# For the text: Get the proportion of mutated sites with at least one mutation:
exome_analysis_normal$intersects |>
  filter(tissue == "blood") |>
  left_join(select(metadata, donor, age) |> distinct() |> dplyr::rename(groupID = donor)) |>
  filter(age > 35) |>
  pull(age) |> mean()

exome_analysis_normal$intersects |>
  filter(tissue == "blood") |>
  left_join(select(metadata, donor, age) |> distinct() |> dplyr::rename(groupID = donor)) |>
  filter(age > 35) |>
  pull(prob) |> mean() * 100

#######
# load pan-cancer BoostDM driver predictions:
#######
pancancer_drivers = fread("processed_data/boostdm/boostdm_genie_cosmic/pancancer_boostDM_intersect.txt.gz")

## TP53
TP53_driver_counts = pancancer_drivers[driver == TRUE & gene_name == "TP53",
                                       .(count = .N), by = .(gene_name, mut_type)]

TP53_analysis = analyze_probability(gene_counts = TP53_driver_counts, analysis_name = "TP53 driver muts", groupby = "donor")
TP53_analysis$plot_list$intersects
save_plots(TP53_analysis$plot_list, path = "plots/coverage_saturation/", name = "TP53")

# analysis for KRAS G12V mutations:
KRAS_G12V_counts = pancancer_drivers[gene_name == "KRAS" & aachange == "G12V",
                                     .(count = .N), by = .(gene_name, mut_type)]
KRAS_G12V_rates = analyze_probability(gene_counts = KRAS_G12V_counts, analysis_name = "KRAS_G12V_muts", groupby = "donor")
KRAS_G12D_counts = pancancer_drivers[gene_name == "KRAS" & aachange == "G12D",
                                     .(count = .N), by = .(gene_name, mut_type)]
KRAS_G12D_rates = analyze_probability(gene_counts = KRAS_G12D_counts, analysis_name = "KRAS_G12D_muts", groupby = "donor")

# another TP53 driver mutations
TP53_R248Q_counts = pancancer_drivers[gene_name == "TP53" & aachange == "R248Q",
                                      .(count = .N), by = .(gene_name, mut_type)]
TP53_R248Q_rates = analyze_probability(gene_counts = TP53_R248Q_counts, analysis_name = "TP53_R248Q_muts", groupby = "donor")

TP53_R175H_counts = pancancer_drivers[gene_name == "TP53" & aachange == "R175H",
                                      .(count = .N), by = .(gene_name, mut_type)]
TP53_R175H_rates = analyze_probability(gene_counts = TP53_R175H_counts, analysis_name = "TP53_R175H_muts", groupby = "donor")

# APC driver mutation
APC_R1450_counts = pancancer_drivers[gene_name == "APC" & aachange == "R1450*",
                                     .(count = .N), by = .(gene_name, mut_type)]
APC_R1450_rates = analyze_probability(gene_counts = APC_R1450_counts, analysis_name = "APC_R1450_muts", groupby = "donor")

# BRAF driver mutation
BRAF_V600E_counts = pancancer_drivers[gene_name == "BRAF" & aachange == "V600E",
                                      .(count = .N), by = .(gene_name, mut_type)]
BRAF_V600E_rates = analyze_probability(gene_counts = BRAF_V600E_counts, analysis_name = "BRAF_V600E_muts", groupby = "donor")

BRAF_values = BRAF_V600E_rates$result_plot_df[tissue == "colon"] |>
  group_by(tissue, category) |>
  summarize(across(c(value, ncells), mean)) |>
  mutate(ncells_mut = value * ncells)

exome_list = exome_analysis$result_plot_df
driver_list = list(KRAS_G12V_rates$result_plot_df,
                   KRAS_G12D_rates$result_plot_df,
                   TP53_R248Q_rates$result_plot_df,
                   TP53_R175H_rates$result_plot_df,
                   APC_R1450_rates$result_plot_df,
                   BRAF_V600E_rates$result_plot_df)
driver_names = c("KRAS G12V", "KRAS G12D", "TP53 R248Q", "TP53 R175H", "APC R1450*", "BRAF V600E")
names(driver_list) = driver_names
drivers = rbindlist(driver_list, idcol = "driver_name")

# TODO: Change the 'x' column
plot_driver_incidence = function(mutation_list, drivers, name, plot_rows = 2, specific_individuals = FALSE) {

  mutation_list$driver_name = NA
  if (specific_individuals[1] == FALSE) {
    print("no specific individual assigned - taking the mean of the cohort")
    driver_summary = drivers[,.(mutated_rate = mean(n_mutated_cells)), by = c("driver_name", "tissue", "category", "x")]
    mut_summary = mutation_list[, .(mutated_rate = mean(n_mutated_cells)), by = c("tissue", "category", "x", "ncells")] |>
      mutate(category_tissue = paste0(tissue, "_", category))
  } else if (all(specific_individuals %in% unique(mutation_list$name))) {
    driver_summary = drivers[name %in% specific_individuals, c("driver_name", "tissue", "category", "n_mutated_cells", "x")] |>
      dplyr::rename(mutated_rate = n_mutated_cells)
    mut_summary = mutation_list[name %in% specific_individuals, c("tissue", "category", "n_mutated_cells", "x", "ncells")] |>
      dplyr::rename(mutated_rate = n_mutated_cells) |>
      mutate(category_tissue = paste0(tissue, "_", category))
  } else {
    print("specific_individuals must be part of the mutation_list and drivers dataframes")
  }

  # add the matching "x" which is the closest to the driver rate observed.
  for (i in 1:nrow(driver_summary)) {
    tissue_select = driver_summary[["tissue"]][i]
    category_select = driver_summary[["category"]][i]
    mutrate_driver = driver_summary[["mutated_rate"]][i]

    closest_x = mut_summary |>
      filter(tissue == tissue_select,
             category == category_select) |>
      mutate(diff = abs(mutated_rate - mutrate_driver)) |>
      filter(diff == min(diff)) |> pull(x)
    driver_summary[["x"]][i] = closest_x
  }

  pl = list()

  # filter drivers for their tissue-specficiity, and make label so that the number of mutations is annotated
  driver_summary = driver_summary |>
    filter(driver_name == "APC R1450*" & tissue == "colon"|
           driver_name %in% c("KRAS G12D", "TP53 R175H", "KRAS G12V", "BRAF V600E")) |>
    group_by(tissue) |>
    mutate(
      nmuts = format(signif(mutated_rate, 2), scientific = FALSE, trim = TRUE, big.mark = ",", drop0trailing = TRUE),
      label = paste0(driver_name, ": ", nmuts),
      max = max(mutated_rate)) |> ungroup()

  pl[["line_plot_nodriver"]] = ggplot(mut_summary, aes(x = x, y = mutated_rate)) +
    geom_line(aes(color = category_tissue), linewidth = 1) +
     facet_nested_wrap(. ~ factor(tissue, c("colon", "lung", "blood")) + category,
                      nrow = plot_rows, nest_line = element_line(linetype = 1),
                      axes = 'all', scales = "free_y") +
    cowplot::theme_cowplot() +
    scale_color_manual(values = tissue_category_colors) +
    scale_x_continuous(labels = label_percent()) +
    theme(legend.position = "none",
          panel.spacing = unit(0.2,"lines"),
          strip.background = element_rect(color = "white", fill = "white")) +
    labs( x = 'All possible single nucleotide variants (SNVs) in the exome sorted by mutation probability',
          y  = 'Number of cells with mutation')
  pl[["line_plot_nodriver"]]


  pl[["line_plot"]]  = pl[["line_plot_nodriver"]] +  ggrepel::geom_text_repel(data = driver_summary, aes(label = label, y =  mutated_rate),
                                                force = 4, size = 4,min.segment.length = 0,
                                                nudge_y = driver_summary$mutated_rate * 1.2 + driver_summary$max * 0.3 + 0.1,
                                                nudge_x = -driver_summary$x*0.5)






  # measurements of unevenness barplot decile:
  mut_deciles = mut_summary |>
    mutate(decile = ntile(x, 10)) %>%
    group_by(tissue, category, category_tissue, decile) %>%
    summarise(sum_mutated = sum(mutated_rate),
              mean_mutrate = mean(mutated_rate),
              median_mutrate = median(mutated_rate), .groups = "drop_last") |>
    mutate(percentage = sum_mutated / sum(sum_mutated),
           percentage_label = paste0(round(percentage * 100, 1), "%"))

  pl[["barplot_decile"]] = ggplot(mut_deciles, aes(x = decile,  y = mean_mutrate, fill = category_tissue)) +
    geom_col() +
    geom_text(aes(label = percentage_label, y = mean_mutrate), vjust = -0.2, position = position_dodge(0.9)) +
    ggh4x::facet_nested_wrap(. ~ tissue + category, nrow = plot_rows, nest_line = element_line(linetype = 1), scales = "free") +
    scale_fill_manual(values = tissue_category_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_x_continuous(breaks = 1:10, labels = scales::label_ordinal()) +
    labs(x = "decile sites ordered by mutagenicity", y = "mean number of mutations/site",
         title = paste("Mutation accumulation by decile:", name),
         subtitle = "percentage of total mutations contributed in %") +
    cowplot::theme_cowplot() +
    coord_cartesian(clip = "off") +
    theme(legend.position = "none",
          panel.spacing = unit(0.2,"lines"),
          strip.background = element_rect(color = "white", fill = "white"))

  # 1. Percentile barplot for the measurements of unevenness barplot percent for the probability
  mut_percent = mut_summary |>
    mutate(percent = ntile(x, 100)) %>%
    group_by(tissue, category, category_tissue, percent) %>%
    mutate(mutated_rate = mutated_rate / ncells ) |>  # transform the mean mutated fraction to the probability (by dividing by the number of cells)
    summarise(sum_mutated = sum(mutated_rate),
              mean_mutrate = mean(mutated_rate),
              median_mutrate = median(mutated_rate), .groups = "drop_last") |>
    mutate(percentage = sum_mutated / sum(sum_mutated),
           percentage_label = paste0(round(percentage * 100, 1), "%"))

  mut_percent_50 = mut_percent |>
    group_by(tissue, category, category_tissue) |>
    mutate(ymax = max(mean_mutrate)) |>
    filter(percent <= 50) |>
    group_by(tissue, category, category_tissue, ymax) |>
    summarize(bottom_50 = sum(percentage),
              y_position = max(mean_mutrate)) |>
    mutate(label = paste("probability of SNV in \nlowest 50%:\n", round(bottom_50 * 100, 1), "%"))

  pl[["barplot_percent_probability"]] = ggplot(mut_percent, aes(x = percent/100,  y = mean_mutrate, fill = category_tissue)) +
    geom_col() +
    ggrepel::geom_text_repel(data = mut_percent |> filter(percent > 97),
                             aes(label = percentage_label, y = mean_mutrate), vjust = -0.2, nudge_x = -0.12) +
    ggpubr::geom_bracket(data = mut_percent_50, aes(xmin = 0, xmax = 0.5,
                                                    y.position = ymax * 0.15, label = label)) +
    ggh4x::facet_nested_wrap(. ~ tissue, nrow = plot_rows, scales = "free", nest_line = element_line(linetype = 1)) +
    scale_fill_manual(values = tissue_category_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks = extended_breaks(4),
                       labels = function(x) x * 1e6) +
    scale_x_continuous(breaks = c(0, 0.5, 1), labels = label_percent()) +
    labs(x = "All possible single nucleotide variants (SNVs) in the exome sorted by mutation probability\n(in percentiles)",
         y = "Probability of mutation\nper cell (x10⁻⁶)") +
    cowplot::theme_cowplot() +
    coord_cartesian(clip = "off") +
    theme(legend.position = "none",
          panel.spacing = unit(0.2,"lines"),
          strip.background = element_rect(color = "white", fill = "white"),
          axis.line = element_line(lineend = "square"))
  pl[["barplot_percent_probability"]]

  # 2:Percentile barplot now for the number of mutated cells
  mut_ncells = mut_summary |>
    mutate(percent = ntile(x, 100)) %>%
    group_by(tissue, category, category_tissue, percent, ncells) %>%
    summarise(sum_mutated = sum(mutated_rate),
              mean_mutrate = mean(mutated_rate),
              median_mutrate = median(mutated_rate), .groups = "drop_last") |>
    mutate(percentage = sum_mutated / sum(sum_mutated),
           percentage_label = paste0(round(percentage * 100, 1), "%"),
           ncells = format_bignum(ncells))

  mut_ncells_label = mut_ncells |>
    group_by(tissue, category, category_tissue, ncells) |>
    mutate(ymax = max(mean_mutrate)) |>
    filter(percent <= 50) |>
    group_by(tissue, category, category_tissue, ymax, ncells) |>
    summarize(min = min(mean_mutrate),
              bottom_50 = mean(mean_mutrate),
              max = max(mean_mutrate),
              y_position = max(mean_mutrate)) |>
    mutate(foldchange = ymax / min,
           label = paste0("cells with SNV\naverage lowest 50%:", format(bottom_50, digits =  2), "\n(range = ", format(min, digits = 1), "-",format(max, digits = 1), ")"))

  pl[["barplot_percent_ncells"]] = ggplot(mut_ncells,
                                          aes(x = percent/100,  y = mean_mutrate, fill = category_tissue)) +
    geom_col() +
    ggrepel::geom_text_repel(data = mut_ncells |> filter(percent > 99),
              aes(label = prettyNum(mean_mutrate, digits = 1, big.mark = ","), y = mean_mutrate),
              vjust = -0.2, nudge_x = -0.12) +
    ggpubr::geom_bracket(data = mut_ncells_label, aes(xmin = 0, xmax = 0.5,
                                                      y.position = ymax * 0.15, label = label)) +
    ggh4x::facet_nested_wrap(. ~ tissue + ncells, nrow = plot_rows, nest_line = element_line(linetype = 1),
                             scales = "free") +
    scale_fill_manual(values = tissue_category_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks =  extended_breaks(4)) +
    scale_x_continuous(breaks = c(0,0.5, 1), labels = label_percent()) +
    labs(x = "exome SNV sites\nsorted by mutation probability",
         y = "Number of cells with SNV") +
    cowplot::theme_cowplot() +
    coord_cartesian(clip = "off") +
    theme(legend.position = "none",
          panel.spacing = unit(0.2,"lines"),
          strip.background = element_rect(color = "white", fill = "white"))
  pl[["barplot_percent_ncells"]]

  # Lorenz plot with gini statistics
  if (name == "normal_exome") {

    gini_table = mut_summary |>
      group_by(tissue, category, category_tissue) |>
      summarize(gini = ineq::Gini(mutated_rate)) |>
      mutate(label = paste(tissue, " Gini:", round(gini,2)))

    lc_values = mut_summary |>
      select(tissue, mutated_rate, x) |>
      pivot_wider(names_from = tissue, values_from = mutated_rate) |>
      arrange(x)

    df_lc = sapply(lc_values[,-1], \(x) ineq::Lc(x)$L) |>
      as.data.frame() |>
      mutate(x = ineq::Lc(lc_values[[2]])$p)
    df_lc_plot = df_lc |>
      pivot_longer(-x, names_to = "tissue") |>
      left_join(gini_table)

    tissue_basic_colors_plot = tissue_basic_colors
    names(tissue_basic_colors_plot) = gini_table[["label"]][match(names(tissue_basic_colors), gini_table$tissue)]

    pl[["lorenz_plot"]] = df_lc_plot |>
      ggplot(aes(x = x, color = label, y = value)) +
      geom_line() +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      coord_cartesian(clip = 'off') +
      scale_color_manual(values = tissue_basic_colors_plot) +
      scale_x_continuous(labels = label_percent(), expand = expansion(mult = c(0, 0)), limits = c(NA,1)) +
      scale_y_continuous(expand = expansion(mult = c(0.0, 0)), limits = c(NA,1)) +
      labs(x = "Exome SNVs sorted by mutation\nprobability",
           y = 'Cumulative share of probability', color = NULL) +
      cowplot::theme_cowplot() +
      theme(legend.position = "inside", legend.position.inside = c(0.05, 0.83))
    pl[["lorenz_plot"]]
  }

  return(pl)
}

# Select individuals for the plot:
individuals = metadata |> select(tissue, category, age, donor) |> distinct() |>
  filter(age > 50 & category %in% c("normal", "non-smoker")) |> arrange(age)

exome_list = exome_analysis$result_plot_df
exome_normal_list = exome_list |> filter(category %in% c("normal", "non-smoker"))
drivers_normal = drivers |> filter(category %in% c("normal", "non-smoker"))
plot_list_normal = plot_driver_incidence(mutation_list = exome_normal_list,
                                         drivers = drivers_normal, name = "normal_exome", plot_rows = 1)


plot_list_normal$barplot_percent_probability

drivers_normal = drivers |> filter(category %in% c("normal", "non-smoker"))
plot_list_normal_individuals = plot_driver_incidence(mutation_list = exome_normal_list,
                                         drivers = drivers_normal, name = "normal_exome", plot_rows = 1,
                                         specific_individuals = c("O340", "PD34215", "KX008"))

ggsave("plots/coverage_saturation/line_plot_nodriver.png", plot_list_normal$line_plot_nodriver, width = 12, height = 4)
ggsave("plots/coverage_saturation/line_plot.png", plot_list_normal$line_plot, width = 12, height = 4)

plot_list_normal = plot_driver_incidence(mutation_list = exome_normal_list,
                                                     drivers = drivers_normal, name = "normal_exome", plot_rows = 1,
                                                     specific_individuals = c("O340", "PD34215", "KX008"))

save_plots(plot_list_normal, "plots/coverage_saturation/", name = "normal_exome_wide", width = 10, height = 5)
save_plots(plot_list_normal, "plots/coverage_saturation/", name = "normal_exome", width = 6.5, height = 5)
ggsave("plots/coverage_saturation/fig1_lorenz_plot.png", plot_list_normal$lorenz_plot, width = 5, height = 4, bg = "white")

# prepare figures for saving and formatting as a single figure
# save Figure 1 panels to be used in main figure script
plot_list_normal = plot_driver_incidence(mutation_list = exome_normal_list,
                                         drivers = drivers_normal, name = "normal_exome", plot_rows = 1)
saveRDS(plot_list_normal_individuals$barplot_percent_probability, "manuscript/figure_panels/figure_1/figure_1C.rds")
saveRDS(plot_list_normal$lorenz_plot, "manuscript/figure_panels/figure_1/figure_1D.rds")

# Figure 2:
figure_2A = prep_plot(exome_analysis_normal$plot_list$plot_saturation_curve_ci, label = "A")
figure_2C = exome_analysis_normal$plot_list$plot_saturation_age +
  theme(legend.position = "inside", legend.position.inside = c(0.75, 0.5)) +
  labs(fill = NULL, title = NULL, subtitle = NULL, x = 'Age (years)') +
  scale_fill_manual(labels = c("blood", "colon", "lung"), values = tissue_category_colors)

# save the individual figures (raw - so margins can be adjusted):
saveRDS(exome_analysis_normal$plot_list$plot_saturation_curve_ci, "manuscript/figure_panels/figure_2/figure_2A.rds")
saveRDS(figure_2C, "manuscript/figure_panels/figure_2/figure_2C.rds")
saveRDS(plot_list_normal$line_plot, "manuscript/figure_panels/figure_2/figure_2D.rds")

# exploration of mutation distribution plots for TP53 driver mutations:
TP53_plots = plot_driver_incidence(mutation_list = TP53_analysis$result_plot_df,
                                   drivers = drivers[grepl("TP53", driver_name)],
                                   name = "TP53")
save_plots(TP53_plots, "plots/coverage_saturation/", "TP53_drivers")