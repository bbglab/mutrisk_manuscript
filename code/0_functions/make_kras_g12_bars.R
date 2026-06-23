make_kras_g12_bars = function(boostdm, sig_donor_rates, metadata, ratios,
                               expected_rates,
                               tissue_select, category_select = "normal",
                               positions = 12) {

  ncells = tissue_ncells$ncells[tissue_ncells$tissue == tissue_select]

  # ---- 1. Build mut_type -> aachange map from ALL tissues ----
  kras_map = rbindlist(boostdm, idcol = "boostdm_tissue", fill = TRUE) |>
    filter(gene_name == "KRAS",
           position %in% positions,
           consequence != "synonymous") |>
    select(mut_type, aachange, position) |>
    distinct()

  # ---- 2. Count genomic sites per mut_type ----
  site_freqs = kras_map[, .N, by = "mut_type"]

  # ---- 3. Filter donors: tissue + category ----
  donors_keep = metadata |>
    filter(tissue == tissue_select, category == category_select) |>
    pull(donor) |>
    unique()

  # ---- 4. Join signature rates, multiply through by site count & ncells ----
  sig_rates = sig_donor_rates |>
    filter(donor %in% donors_keep) |>
    inner_join(site_freqs, by = "mut_type", relationship = "many-to-many") |>
    filter(!is.na(signature))

  # ---- 5. Use expected_rates total, compute global signature proportions ----
  ratio_val = ratios |>
    filter(tissue == tissue_select, category == category_select,
           gene_name == "KRAS") |>
    pull(ratio)
  if (length(ratio_val) == 0) ratio_val = 1

  # Total expected cells per mut_type (avg donor)
  exp_mut = expected_rates |>
    filter(tissue == tissue_select, category == category_select) |>
    left_join(select(metadata, sampleID, donor, category),
              by = c("sampleID", "category")) |>
    group_by(mut_type, donor) |>
    summarize(mle = mean(mle), .groups = "drop") |>
    group_by(mut_type) |>
    summarize(exp_mle = mean(mle), .groups = "drop")

  # Global signature proportions (averaged across donors)
  sig_props = sig_rates |>
    group_by(mut_type, signature) |>
    summarize(sig_mle = mean(mle), .groups = "drop") |>
    group_by(mut_type) |>
    mutate(sig_prop = sig_mle / sum(sig_mle)) |>
    ungroup() |>
    select(mut_type, signature, sig_prop)

  # Combine: total cells * signature proportion
  df = exp_mut |>
    left_join(site_freqs, by = "mut_type") |>
    inner_join(sig_props, by = "mut_type", relationship = "many-to-many") |>
    mutate(mle = exp_mle * sig_prop * N * ncells * ratio_val) |>
    select(-exp_mle, -sig_prop, -N) |>
    left_join(distinct(kras_map, mut_type, aachange), by = "mut_type") |>
    group_by(aachange, signature) |>
    summarize(mle = sum(mle), .groups = "drop") |>
    mutate(signature = factor(signature))

  # ---- 7. Trinucleotide-context label per aachange ----
  label_map = kras_map |>
    select(aachange, mut_type) |>
    distinct() |>
    group_by(aachange) |>
    summarize(label = paste(mut_type, collapse = "\n"), .groups = "drop")

  # Bar-top height for label positioning
  bar_tops = df |>
    group_by(aachange) |>
    summarize(y_pos = sum(mle), .groups = "drop")
  label_map = left_join(label_map, bar_tops, by = "aachange")

  # ---- 8. Order bars: by position, then alphabetically ----
  aachange_order = kras_map |>
    select(aachange, position) |>
    distinct() |>
    arrange(position, aachange) |>
    pull(aachange)

  df$aachange      = factor(df$aachange,      levels = aachange_order)
  label_map$aachange = factor(label_map$aachange, levels = aachange_order)

  # ---- 9. Plot with controlled legend ----
  # Subset colors to only those present in this specific tissue/category
  current_sigs = levels(df$signature)
  active_colors = sig_colors[names(sig_colors) %in% current_sigs]

  ggplot(df, aes(x = aachange, y = mle, fill = signature)) +
    geom_col() +
    geom_text(data = label_map,
              aes(x = aachange, y = y_pos, label = label),
              vjust = -0.3, size = 2.5, inherit.aes = FALSE) +
    scale_fill_manual(values = active_colors, drop = TRUE) +
    theme_cowplot() +
    theme(legend.position = "right") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25)),
                       labels = scales::label_comma()) +
    coord_cartesian(clip = "off") +
    labs(x = "AA change",
         y = "Number of cells with\nKRAS mutation",
         subtitle = paste(tissue_select, category_select, sep = " - "))
}
