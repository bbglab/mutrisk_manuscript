make_gene_barplot = function(boostdm, ratios, expected_rates,  gene_of_interest,
                             tissue_select = "colon", tissue_name = NULL,
                             category_select = "normal",
                             cell_probabilities = FALSE, individual = FALSE, older_individuals = TRUE,
                             include_hotspots = FALSE) {

  if (is.null(tissue_name)) {tissue_name = tissue_select}

  mr_drivers = merge_mutrisk_drivers(boostdm, ratios, expected_rates, gene_of_interest, tissue_select, tissue_name, category_select, cell_probabilities,
                                     individual, filter_age = FALSE)
  diff_high_low =  max(mr_drivers$expected_gene_muts$mle) / min(mr_drivers$expected_gene_muts$mle)
  print(paste0("Gene: ", gene_of_interest, ", Tissue: ", tissue_select, ", Individual: ", individual))
  print(paste0("difference max/min probability: ", round(diff_high_low, 1)))

  expected_gene_muts = mr_drivers$expected_gene_muts
  label = mr_drivers$label

  y_label = "Number of cells with mutation"
  if (cell_probabilities == TRUE) {
    ncells_select = 1
    y_label = "Probability of mutation\n per cell(x10⁻⁶)"
  }

  if (max(expected_gene_muts$position, na.rm = TRUE) > Inf) { # for now set the level to Inf to allow for large genes
    expected_gene_muts = expected_gene_muts |>
      mutate(position = (position - 1) %/% 5 + 1,
             position = position * 5) |>
      group_by(position, tissue, mut_type, driver)  |>
      summarise(mle = sum(mle, na.rm = TRUE), .groups = "drop")
    x_label = "AA position (5AA bins)"
  } else { x_label = "AA position"}

  expected_gene_muts_label = left_join(expected_gene_muts, mutrisk:::triplet_match_substmodel, by = "mut_type")

  # way to make the plot extend both upper and lower axes
  pl = ggplot(expected_gene_muts_label,
              aes(x = position, y = mle)) +
    geom_col(aes(fill = type)) +
    scale_fill_manual(values = mutrisk::COLORS6) +
    theme_cowplot() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = x_label, y = y_label, title = gene_of_interest, subtitle = label, fill = NULL)

  if (cell_probabilities == TRUE) {
    pl = pl + scale_y_continuous(expand = expansion(mult = c(0, 0.2)), labels = function(x) x * 1e6)
  }

  if (include_hotspots[1] != FALSE) {

    if (class(include_hotspots) != "numeric") {
      stop("'include_hotspots' needs to be a numeric vector indicating hotspot positions in the gene\n
           optional: include names of the hotspots to indicate AA positions")
    }

    # add layer behind the main bars of the mutation plot
    hotspot_layer = geom_vline(xintercept = include_hotspots, color = "grey")
    pl$layers = c(hotspot_layer, pl$layers)

    df_hotspots = data.frame(label = names(include_hotspots), x = include_hotspots)
    pl =  pl + ggrepel::geom_text_repel(data = df_hotspots, aes(label = label, x = x, y = Inf), angle = 90, bg.color =  "white", size = 3,
                                    bg.r = 0.25, force = 0.2) +  coord_cartesian(clip = "off")
    pl

    }

  pl
}