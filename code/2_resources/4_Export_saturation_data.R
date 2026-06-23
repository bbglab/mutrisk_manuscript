### 4_Export_saturation_data.R
# Export saturation analysis results to flat files for mutrisk-web.
#
# Run from the mutrisk_manuscript project root:
#   Rscript code/2_resources/4_Export_saturation_data.R
#
# Produces four files in processed_data/saturation/:
#   exome_sites_sample.tsv.gz    — 100k random exome sites from wintr::RefCDS_WGS (reproducible seed)
#   sorted_rates.tsv.gz          — per-donor sorted per-site rates (Fig 1C / Fig 1D source data)
#   saturation_curves.tsv.gz     — per-donor saturation curves over ncells range (Fig 2A source data)
#   saturation_intersects.tsv.gz — saturation probability at tissue-specific ncells (Fig 2D source data)
#
# Column names match the mutrisk-web database schema (sql/schema.sql).
# See SESSION_NOTES_38.md and mutrisk-web/CLAUDE.md for design decisions.

library(data.table)
library(tidyverse)
library(mutrisk)
source("code/0_functions/analysis_variables.R")
set.seed(1234)

# ---------------------------------------------------------------------------
# 1. Load metadata (same pattern as 1_Figure1_and_2.R)
# ---------------------------------------------------------------------------
md_files <- list.files("processed_data/", recursive = TRUE, pattern = "_metadata",
                       full.names = TRUE)
names(md_files) <- gsub("_.*", "", basename(md_files))
metadata <- lapply(md_files, fread) |>
  rbindlist(idcol = "tissue", use.names = TRUE, fill = TRUE) |>
  dplyr::select(any_of(c("tissue", "sampleID", "category", "age", "donor"))) |>
  dplyr::distinct()

message("Metadata loaded: ", nrow(metadata), " rows across ",
        length(unique(metadata$tissue)), " tissues")

# ---------------------------------------------------------------------------
# 2. Build gene_counts from wintr::RefCDS_WGS  (same 6 lines as Fig 1_and_2.R)
# ---------------------------------------------------------------------------
genes <- lapply(wintr::RefCDS_WGS, \(x) as.data.table(x[["L"]][, 1:4]))
names(genes) <- sapply(wintr::RefCDS_WGS, \(x) x[["gene_name"]])
gene_counts <- rbindlist(genes, idcol = "gene_name")
gene_counts$mut_type <- rep(mutrisk::triplet_match_substmodel$mut_type,
                            nrow(gene_counts) / 192)
colnames(gene_counts)[-c(1, 6)] <- c("syn", "mis", "nonsense", "splice")
gene_counts <- gene_counts |>
  pivot_longer(-c(gene_name, mut_type), names_to = "consequence", values_to = "count") |>
  filter(count != 0) |>
  setDT()
set.seed(1234)
gene_counts_sample <- sample_n(gene_counts, 1e5)

# ---------------------------------------------------------------------------
# 3. Data-only version of analyze_probability (no plots, no cowplot dependency)
#    Uses get_adjusted_rates() and get_prob_mutated_range() from
#    code/0_functions/plot_mut_prob.R (auto-sourced by analysis_variables.R).
# ---------------------------------------------------------------------------
analyze_probability_data <- function(gene_counts, filter_normal = FALSE) {

  prob_rates_tissue        <- list()
  prob_intersects_ci_tissue <- list()
  result_plot_list         <- list()

  for (tissue in c("colon", "lung", "blood")) {
    message("  tissue: ", tissue)

    expected_rates <- fread(paste0("processed_data/", tissue, "/", tissue,
                                   "_expected_rates.tsv.gz"))
    if (isTRUE(filter_normal)) {
      expected_rates <- expected_rates |>
        filter(category %in% c("normal", "non-smoker"))
    }

    expected_rates <- expected_rates |>
      left_join(metadata, by = c("sampleID", "category")) |>
      as.data.frame()

    expected_rates <- expected_rates |>
      mutate(groupid = expected_rates[["donor"]]) |>
      group_by(category, mut_type, groupid) |>
      summarize(across(c(mle, cilow, cihigh), mean), .groups = "drop") |>
      setDT()

    ratios <- fread(paste0("processed_data/", tissue, "/", tissue,
                           "_mut_ratios.tsv.gz"))

    # -----------------------------------------------------------------------
    # 3a. sorted per-site rates matrix → result_plot_df
    # -----------------------------------------------------------------------
    ncells_mid <- tissue_ncells$ncells[tissue_ncells$tissue == tissue]

    category_list <- list()
    for (category_select in unique(expected_rates$category)) {
      message("    category: ", category_select)
      ratios_cat          <- ratios[category == category_select, ]
      expected_rates_cat  <- expected_rates[category == category_select, ]
      expected_rates_list <- split(expected_rates_cat, expected_rates_cat$groupid)

      m <- matrix(NA, ncol = length(expected_rates_list), nrow = nrow(gene_counts))
      colnames(m) <- names(expected_rates_list)

      for (groupid in names(expected_rates_list)) {
        genome_rates     <- get_adjusted_rates(
          expected_rates_sample = expected_rates_list[[groupid]],
          gene_counts           = gene_counts,
          ratios_cat            = ratios_cat)
        m[, groupid] <- genome_rates$adjusted_rates
      }

      m_sorted <- apply(m, 2, sort)
      if (!is.matrix(m_sorted)) m_sorted <- t(m_sorted)

      # Bin into 1000 bins if ≥1e5 rows (same logic as 1_Figure1_and_2.R)
      if (nrow(m_sorted) >= 1e5) {
        binsize <- 1000
        group   <- rep(seq_len(nrow(m_sorted) %/% binsize), each = binsize)
        output_m <- rowsum(m_sorted, group) / binsize
      } else {
        output_m <- m_sorted
      }
      colnames(output_m) <- colnames(m)
      category_list[[category_select]] <- as.data.table(output_m)
    }

    binsize <- nrow(category_list[[1]])

    plot_data   <- lapply(category_list, pivot_longer, everything()) |>
      rbindlist(idcol = "category")
    result_plot <- plot_data |>
      arrange(category, name, value) |>
      group_by(category, name)

    result_plot_df <- result_plot |>
      ungroup() |>
      mutate(n_mutated_cells = value * ncells_mid,
             ncells = ncells_mid,
             x = rep(seq_len(binsize), dplyr::n_groups(result_plot)) / 1e4,
             x = x / max(x))

    result_plot_list[[tissue]] <- result_plot_df

    # -----------------------------------------------------------------------
    # 3b. Saturation curves + intersects at tissue ncells_ci
    # -----------------------------------------------------------------------
    ncells_ci_row  <- tissue_ncells_ci[tissue_ncells$tissue == tissue, ]
    ncells_ci_vals <- as.numeric(ncells_ci_row[-1])   # drops "tissue" column
    ncells_ci_names <- names(ncells_ci_row[-1])       # "high_estimate", "mid_estimate", "low_estimate"

    prob_rates        <- list()
    prob_intersects_ci <- list()

    for (category_select in unique(expected_rates$category)) {
      ratios_cat          <- ratios[category == category_select, ]
      expected_rates_cat  <- expected_rates[category == category_select, ]
      expected_rates_list <- split(expected_rates_cat, expected_rates_cat$groupid)

      list_probs         <- list()
      list_intersects_ci <- list()

      for (name in names(expected_rates_list)) {
        list_probs[[name]] <- get_prob_mutated_range(
          expected_rates_sample = expected_rates_list[[name]],
          gene_counts           = gene_counts,
          ratios_cat            = ratios_cat,
          step                  = 0.2)

        list_intersects_ci[[name]] <- get_prob_mutated_range(
          expected_rates_sample = expected_rates_list[[name]],
          gene_counts           = gene_counts,
          ratios_cat            = ratios_cat,
          range                 = ncells_ci_vals) |>
          mutate(name = ncells_ci_names)
      }

      prob_rates[[category_select]]         <- rbindlist(list_probs, idcol = "groupID")
      prob_intersects_ci[[category_select]] <- rbindlist(list_intersects_ci, idcol = "groupID")
    }

    prob_rates_tissue[[tissue]]         <- rbindlist(prob_rates, idcol = "category")
    prob_intersects_ci_tissue[[tissue]] <- rbindlist(prob_intersects_ci, idcol = "category")
  }

  # Combine across tissues
  bind_dfs <- function(df_list) {
    rbindlist(df_list, idcol = "tissue") |>
      mutate(tissue = factor(tissue, levels = c("colon", "lung", "blood")))
  }

  list(
    result_plot_df  = bind_dfs(result_plot_list),
    rates           = bind_dfs(prob_rates_tissue),
    intersects_ci   = bind_dfs(prob_intersects_ci_tissue)
  )
}

# ---------------------------------------------------------------------------
# 4. Run the analysis (normal-like categories only, same as exome_analysis_normal)
# ---------------------------------------------------------------------------
message("Running saturation analysis (normal-like categories)...")
res <- analyze_probability_data(gene_counts = gene_counts_sample, filter_normal = TRUE)
message("Analysis done.")

# ---------------------------------------------------------------------------
# 5. Post-process and add age to intersects_ci
# ---------------------------------------------------------------------------
metadata_age <- metadata |>
  dplyr::select(donor, age) |>
  dplyr::distinct()

# sorted_rates: rename 'name' → 'donor' for clarity
sorted_rates <- res$result_plot_df |>
  dplyr::rename(donor = name) |>
  dplyr::select(tissue, category, donor, x, value, n_mutated_cells, ncells) |>
  dplyr::mutate(tissue = as.character(tissue))

# saturation_curves: rename 'groupID' → 'donor'
saturation_curves <- res$rates |>
  dplyr::rename(donor = groupID) |>
  dplyr::select(tissue, category, donor, ncells, prob) |>
  dplyr::mutate(tissue = as.character(tissue))

# saturation_intersects: add age, rename columns
saturation_intersects <- res$intersects_ci |>
  dplyr::rename(donor = groupID, estimate = name) |>
  dplyr::left_join(metadata_age, by = "donor") |>
  dplyr::select(tissue, category, donor, age, ncells, prob, estimate) |>
  dplyr::mutate(tissue = as.character(tissue))

# ---------------------------------------------------------------------------
# 6. Export
# ---------------------------------------------------------------------------
out_dir <- "processed_data/saturation"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Exporting exome_sites_sample.tsv.gz ...")
fwrite(gene_counts_sample, file.path(out_dir, "exome_sites_sample.tsv.gz"))

message("Exporting sorted_rates.tsv.gz ...")
fwrite(sorted_rates, file.path(out_dir, "sorted_rates.tsv.gz"))

message("Exporting saturation_curves.tsv.gz ...")
fwrite(saturation_curves, file.path(out_dir, "saturation_curves.tsv.gz"))

message("Exporting saturation_intersects.tsv.gz ...")
fwrite(saturation_intersects, file.path(out_dir, "saturation_intersects.tsv.gz"))

# ---------------------------------------------------------------------------
# 7. Sanity-check row counts
# ---------------------------------------------------------------------------
message("\n=== Row counts ===")
message("  exome_sites_sample:    ", nrow(gene_counts_sample), " rows")
message("  sorted_rates:          ", nrow(sorted_rates), " rows")
message("  saturation_curves:     ", nrow(saturation_curves), " rows")
message("  saturation_intersects: ", nrow(saturation_intersects), " rows")
message("\nColumn names:")
message("  sorted_rates:          ", paste(names(sorted_rates), collapse = ", "))
message("  saturation_curves:     ", paste(names(saturation_curves), collapse = ", "))
message("  saturation_intersects: ", paste(names(saturation_intersects), collapse = ", "))
message("\nEstimate values in saturation_intersects: ",
        paste(unique(saturation_intersects$estimate), collapse = ", "))
message("\n✓ Saturation data exported to ", out_dir)
