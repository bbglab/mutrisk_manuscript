# Mirror plot: IntOGen observed % (top) vs MutRisk expected cells (bottom)
# for a single tissue.  Uses facet_grid2 with negated Y-axis for the expected half.
#
# Requires: analysis_variables.R sourced, ggh4x loaded

SUBS_CLASSES = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")

# Tissue colours used for facet strip backgrounds
MIRROR_TISSUE_COLORS = c(colon = "#3ca951", lung = "#4269d0", blood = "#ff725c")

## Subtitle for complete tumor type (e.g. COADREAD, LUAD, AML)
FULL_CANCER_TYPE = c(COADREAD = "Colorectal Adenocarcinoma", LUSC = "Lung Squamous Cell Carcinoma", AML = "Acute Myeloid Leukemia")

make_kras_mirror_plot = function(intogen_data, boostdm, expected_rates, ratios,
                                  tissue_select, category_select = "normal") {

  tissue_color = MIRROR_TISSUE_COLORS[[tissue_select]]

  # Cancer type label (e.g. COADREAD, LUAD, AML)
  cancer_type = intogen_data |>
    filter(tissue == tissue_select) |>
    pull(CANCER_TYPE) |> unique() |> na.omit()
  if (length(cancer_type) == 0) cancer_type = tissue_select

  # ---- IntOGen: pivot wide columns (C>A .. T>G) to long ----
  intogen_strip = "IntOGen\nSamples"
  intogen_tissue = intogen_data |>
    filter(tissue == tissue_select) |>
    pivot_longer(all_of(SUBS_CLASSES), names_to = "type", values_to = "mrate") |>
    mutate(tissue_category = intogen_strip) |>
    select(position, type, tissue_category, mrate)

  # ---- MutRisk: reuse make_gene_barplot data flow ----
  # (merge_mutrisk_drivers + triplet_match, no summarisation)
  mr = merge_mutrisk_drivers(boostdm, ratios, expected_rates,
                              gene_of_interest  = "KRAS",
                              tissue_select     = tissue_select,
                              tissue_name       = tissue_select,
                              category_select   = category_select,
                              cell_probabilities = FALSE,
                              filter_age        = TRUE)

  # merge_mutrisk_drivers already multiplies by ncells internally
  expected_strip = "Estimated\ncells"
  mutrisk_tissue = mr$expected_gene_muts
  mutrisk_tissue = left_join(mutrisk_tissue, mutrisk:::triplet_match_substmodel,
                              by = "mut_type")
  mutrisk_tissue$mrate = mutrisk_tissue$mle
  mutrisk_tissue = mutrisk_tissue |>
    select(position, type, mrate) |>
    mutate(tissue_category = expected_strip)

  # ---- Combine, negate expected half ----
  category_levels = c(intogen_strip, expected_strip)

  df_mirror = bind_rows(intogen_tissue, mutrisk_tissue) |>
    mutate(
      tissue_category = factor(tissue_category, levels = category_levels),
      mrate = ifelse(tissue_category == expected_strip, -mrate, mrate))

  # Invisible point for symmetric y-limits
  df_point = df_mirror |>
    group_by(tissue_category) |>
    summarize(mrate = max(abs(mrate)) * 1.1, .groups = "drop") |>
    mutate(
      position = max(df_mirror$position, na.rm = TRUE),
      mrate    = ifelse(tissue_category == expected_strip, -mrate, mrate))

  # ---- Hotspot annotations (G12, G13, G14) --------
  # Both strips: text left of bar, horizontal arrow from text to bar top.
  # Text is offset slightly below the bar tip to avoid clipping at the plot edge.
  hotspot_positions <- c(12, 13, 14)
  hotspot_labels    <- c("G12", "G13", "G14")
  x_offset          <- 10
  arrow_gap         <- 5
  y_offset          <- 0.95  # multiply bar_top to pull text downward from edge

  top_data <- df_mirror |>
    filter(tissue_category == intogen_strip, position %in% hotspot_positions) |>
    group_by(position) |>
    summarise(bar_top = sum(mrate), .groups = "drop") |>
    mutate(
      tissue_category = factor(intogen_strip, levels = category_levels),
      label           = hotspot_labels[match(position, hotspot_positions)],
      x_label         = position - x_offset,
      x_arrow_start   = position - arrow_gap,
      x_arrow_end     = position,
      y_label         = ifelse(
                          position == 14 & bar_top < 0.08, bar_top + 2, 
                          ifelse(position == 14, bar_top + (bar_top * 1.8), bar_top * y_offset)),
      y_arrow         = bar_top
    )

  bottom_data <- df_mirror |>
    filter(tissue_category == expected_strip, position %in% hotspot_positions) |>
    group_by(position) |>
    summarise(bar_top = sum(mrate), .groups = "drop") |>
    mutate(
      tissue_category = factor(expected_strip, levels = category_levels),
      label           = hotspot_labels[match(position, hotspot_positions)],
      # G12: pull label higher (closer to zero) and horizontally closer to bar
      x_label         = position - x_offset,
      x_arrow_start   = position - arrow_gap,
      x_arrow_end     = position,
      y_label         = ifelse(position == 12, bar_top * 1.5, bar_top),
      y_arrow         = bar_top
    )

  ggplot(df_point, aes(x = position, y = mrate)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.5) +
    geom_point(color = "white") +
    geom_col(data = df_mirror, aes(fill = type)) +
    # Top strip: text left, horizontal arrow rightward
    geom_text(data = top_data,
              aes(x = x_label, y = y_label, label = label,
                  tissue_category = tissue_category),
              color = "black", size = 3, fontface = "plain", hjust = 0,
              inherit.aes = FALSE) +
    geom_segment(data = top_data,
                 aes(x = x_arrow_start, y = y_label,
                     xend = x_arrow_end, yend = bar_top,
                     tissue_category = tissue_category),
                 arrow = arrow(length = unit(0.04, "inches"), type = "closed"),
                 color = "black", linewidth = 0.4,
                 inherit.aes = FALSE) +
    # Bottom strip: text left, horizontal arrow rightward
    geom_text(data = bottom_data,
              aes(x = x_label, y = y_label, label = label,
                  tissue_category = tissue_category),
              color = "black", size = 3, fontface = "plain", hjust = 0,
              inherit.aes = FALSE) +
    geom_segment(data = bottom_data,
                 aes(x = x_arrow_start, y = y_label,
                     xend = x_arrow_end, yend = bar_top,
                     tissue_category = tissue_category),
                 arrow = arrow(length = unit(0.04, "inches"), type = "closed"),
                 color = "black", linewidth = 0.4,
                 inherit.aes = FALSE) +
    facet_grid2(tissue_category ~ ., scales = "free",
                strip = strip_themed(
                  background_y = elem_list_rect(
                    fill = c(tissue_color, tissue_color)),
                  text_y = elem_list_text(
                    colour = c("white", "white"),
                    face = "plain"))) +
    scale_fill_manual(values = COLORS6) +
    theme_cowplot() +
    theme(legend.position = "none",
          panel.spacing.y = unit(0, "mm"),
          strip.background.y = element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_line(color = "grey85"),
          axis.line.y = element_line(color = "grey85"),
          axis.text.y = element_text(margin = margin(r = 4)),
          strip.text.y = element_text(margin = margin(r = 4, l = 4))) +
    labs(y = NULL, x = "AA position",
         title = "KRAS", subtitle = paste0(FULL_CANCER_TYPE[cancer_type], "/", tissue_select)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0)),
                       breaks = scales::breaks_extended(n = 3),
                       labels = abs)
}
