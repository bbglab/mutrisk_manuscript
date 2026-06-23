# Age-vs-cells scatter plots for KRAS G12 driver mutations.
# Uses facet_nested with nestline matching Figure 3D style.
#
# Returns a named list of 3 ggplots:
#   colon (1 facet), lung (3 facets), blood (1 facet)
#
# Requires: analysis_variables.R sourced, ggh4x loaded

make_kras_g12_scatters = function(boostdm, expected_rates, ratios, metadata) {

  # Per-donor KRAS G12 cell estimates by tissue x category
  dotplot_list = list()
  for (i in 1:nrow(color_df)) {

    category_select = color_df$category[i]
    tissue_select   = color_df$tissue[i]

    dotplot_list[[i]] = merge_mutrisk_drivers(boostdm, ratios, expected_rates,
                                               gene_of_interest = "KRAS",
                                               tissue_select    = tissue_select,
                                               tissue_name      = tissue_select,
                                               category_select  = category_select,
                                               filter_age       = FALSE,
                                               individual       = "all")[[1]] |>
      filter(position == 12) |>
      group_by(donor) |>
      summarize(across(c(mle, cilow, cihigh), sum), .groups = "drop") |>
      mutate(tissue = tissue_select, category = category_select)
  }

  dotplot_df = rbindlist(dotplot_list) |>
    mutate(tissue = factor(tissue, levels = c("colon", "lung", "blood")),
           tissue_category = paste0(tissue, "_", category))

  # Keep only: colon_normal, all lung categories, blood_normal
  dotplot_df = dotplot_df |>
    filter(category == "normal" | tissue == "lung")

  df_total_muts = dotplot_df |>
    left_join(metadata |>
              mutate(tissue = as.character(tissue)) |>
              select(tissue, donor, age) |> distinct(),
              by = c("tissue", "donor")) |>
    mutate(category = factor(category, levels = unique(color_df$category)))

  # Build one plot per tissue (colon, lung, blood), faceted by category
  tissue_plots = list()
  tissue_order = c("colon", "lung", "blood")

  for (i in 1:3) {

    select_tissue = tissue_order[[i]]
    df_tissue = df_total_muts |>
      filter(tissue == select_tissue)

    plt = ggplot(df_tissue, aes(x = age, y = mle, color = tissue_category)) +
      geom_pointrange(aes(ymin = cilow, ymax = cihigh)) +
      facet_nested(. ~ tissue + category, axes = "y", remove_labels = "y") +
      scale_color_manual(values = tissue_category_colors) +
      scale_y_continuous(labels = scales::label_comma(), limits = c(0, NA)) +
      theme_cowplot() +
      theme(panel.grid = element_blank(),
            strip.background = element_blank(),
            legend.position = "none",
            ggh4x.facet.nestline = element_line()) +
      labs(x = "Age (years)",
           y = "Number of cells with\nKRAS G12 driver mutation")

    tissue_plots[[select_tissue]] = plt
  }

  # Remove y-axis label from lung and blood (share with colon)
  tissue_plots[["lung"]]  = tissue_plots[["lung"]]  + labs(y = NULL)
  tissue_plots[["blood"]] = tissue_plots[["blood"]] + labs(y = NULL)

  tissue_plots
}
