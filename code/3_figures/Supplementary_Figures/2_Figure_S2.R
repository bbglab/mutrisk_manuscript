# 2_Figure_S2.R
# Purpose: Generate supplementary figure S2 -- signature refitting analysis
#   across blood, colon, and lung cohorts.
source("code/0_functions/analysis_variables.R")

# load metadata
md_files = list.files("processed_data/", recursive = TRUE, pattern = "_metadata",
                      full.names = TRUE)
names(md_files) = gsub("_.*", "", basename(md_files))
metadata = lapply(md_files, fread) |>
  rbindlist(idcol = "tissue", use.names = TRUE, fill = TRUE) |>
  dplyr::select(any_of(c("tissue", "sampleID", "category", "age", "donor"))) |>
  dplyr::distinct()

sig_contri_files = list.files("processed_data/", recursive = TRUE,
                                     pattern = "signature_contributions",
                                     full.names = TRUE)

names(sig_contri_files) = paste0(
  str_split_i(sig_contri_files,pattern = "/", i = c(3)), "_",
  str_split_i(sig_contri_files,pattern = "/", i = c(4)))



means_sig_contrib = sig_contribution_plot = list()
for (name in names(sig_contri_files)) {

  file = sig_contri_files[name]
  sig_per_sample = fread(file) |>
    left_join(metadata) |>
    group_by(donor) |>
    summarize(across(starts_with("SBS"), mean))

  mean_sig_contrib = sig_per_sample |>
    column_to_rownames("donor") |> as.matrix() |>
    prop.table(1) |>
    as.data.frame() |>
    rownames_to_column("donor") |>
    pivot_longer(starts_with("SBS"), names_to = "signature", values_to = "contribution")

  means_sig_contrib[[name]] = mean_sig_contrib |>
    group_by(signature) |>
    summarize(sd = sd(contribution),
              contribution = mean(contribution))

  sig_contribution_plot[[name]] = means_sig_contrib[[name]] |>
    ggplot(aes(x = signature, y = contribution)) +
    geom_col() +
    geom_point(data = mean_sig_contrib) +
    #geom_errorbar(aes(ymin = contribution, ymax = contribution + sd), width = 0.2) +
    labs(title = name, y = "relative signature contribution") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

sig_contribution_plot = wrap_plots(sig_contribution_plot)
ggsave("manuscript/Supplementary_Figures/Figure_S2/Figure_S2_all.png",
       sig_contribution_plot, width = 15, height = 10)


order(parse_number(means_sig_contrib[[2]]$signature))

mean_signatures = means_sig_contrib |>
  rbindlist(idcol = "name")

sigs = unique(mean_signatures$signature)
levels = sigs[order(parse_number(sigs), decreasing = FALSE)]
mean_signatures_ordered = mean_signatures |>
  mutate(name = gsub("_", " ", name),
         signature = factor(signature, levels = levels),
         name = factor(name, levels = c("blood normal", "colon normal", "colon IBD","colon POLD1",
                                        "colon POLE", "lung non-smoker", "lung ex-smoker", "lung smoker")))

colors = setNames(ggsci::pal_igv()(n_distinct(mean_signatures_ordered$signature)),
                  unique(mean_signatures_ordered$signature))

plot_blood = ggplot(mean_signatures_ordered |> filter(grepl("blood", name)),
                                aes(x = name, y = contribution, fill = signature)) +
  geom_col() +
  labs(x = NULL) +
  theme_cowplot() +
  scale_fill_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  labs(y = "relative signature contribution", title = "Blood cohort")

plot_colon = ggplot(mean_signatures_ordered |> filter(grepl("colon", name)),
                                aes(x = name, y = contribution, fill = signature)) +
  geom_col() +
  labs(x = NULL) +
  theme_cowplot() +
  scale_fill_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  labs(y = "relative signature contribution", title = "Colon cohorts")

plot_lung = ggplot(mean_signatures_ordered |> filter(grepl("lung", name)),
                                aes(x = name, y = contribution, fill = signature)) +
  geom_col() +
  labs(x = NULL) +
  theme_cowplot() +
  scale_fill_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  labs(y = "relative signature contribution", title =   "Lung cohorts")

total_plot = plot_blood + plot_colon + plot_lung + plot_layout(widths = c(1,3,2.3))
ggsave("manuscript/Supplementary_Figures/Figure_S2/Figure_S2_raw.pdf" , total_plot,
       width = 12, height = 5, bg = "white")


##### Make an additional plot with the normalized signature contribution through age:
# additional figure on the signature contribution across different samples.
# For example, normal colon, lung and blood:

colon_metadata = fread("processed_data/colon/colon_metadata.tsv")
colon_sigs = fread("processed_data/colon/normal/signature_contributions.tsv")

colon_mean_sig_age = colon_sigs |>
  left_join(colon_metadata) |>
  group_by(donor, age) |>
  summarize(across(starts_with("SBS"), mean)) |>
  pivot_longer(starts_with("SBS"), names_to = "signature", values_to = "contribution") |>
  mutate(donor_age = paste0(donor, " (", age, ")")) |>
  arrange(age) |> ungroup() |>
  mutate(donor_age = fct_reorder(donor_age, age))

colon_mean_sig_age_rel = colon_mean_sig_age |>
  group_by(donor_age) |>
  mutate(rel_contribution = contribution / sum(contribution),
         signature = factor(signature, levels = levels(mean_signatures_ordered$signature)))

donor_contri_colon = ggplot(colon_mean_sig_age_rel, aes(x = donor_age, y = rel_contribution, fill = signature )) +
  geom_col() +
  scale_fill_manual(values = colors) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none")  +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  labs(x = "Donor (Age in years)", y = "Relative signature contribution")


### Lung
lung_metadata = fread("processed_data/lung/lung_metadata.tsv")
lung_sigs = fread("processed_data/lung/non-smoker/signature_contributions.tsv")

lung_mean_sig_age = lung_sigs |>
  left_join(lung_metadata) |>
  group_by(donor, age) |>
  summarize(across(starts_with("SBS"), mean)) |>
  pivot_longer(starts_with("SBS"), names_to = "signature", values_to = "contribution") |>
  mutate(donor_age = paste0(donor, " (", round(age), ")")) |>
  arrange(age) |> ungroup() |>
  mutate(donor_age = fct_reorder(donor_age, age))

lung_mean_sig_age_rel = lung_mean_sig_age |>
  group_by(donor_age) |>
  mutate(rel_contribution = contribution / sum(contribution),
         signature = factor(signature, levels = levels(mean_signatures_ordered$signature)))

donor_contri_lung = ggplot(lung_mean_sig_age_rel,
                           aes(x = donor_age, y = rel_contribution, fill = signature )) +
  geom_col() +
  scale_fill_manual(values = colors) +
  theme_cowplot() +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none")  +
  labs(x = "Donor (Age in years)", y = "Relative signature contribution")


### Blood
blood_metadata = fread("processed_data/blood/blood_metadata.tsv")
blood_sigs = fread("processed_data/blood/normal/signature_contributions.tsv")

blood_mean_sig_age = blood_sigs |>
  left_join(blood_metadata) |>
  group_by(donor, age) |>
  summarize(across(starts_with("SBS"), mean)) |>
  pivot_longer(starts_with("SBS"), names_to = "signature", values_to = "contribution") |>
  mutate(donor_age = paste0(donor, " (", round(age), ")")) |>
  arrange(age) |> ungroup() |>
  mutate(donor_age = fct_reorder(donor_age, age))

blood_mean_sig_age_rel = blood_mean_sig_age |>
  group_by(donor_age) |>
  mutate(rel_contribution = contribution / sum(contribution),
         signature = factor(signature, levels = levels(mean_signatures_ordered$signature)))

donor_contri_blood = ggplot(blood_mean_sig_age_rel, aes(x = donor_age, y = rel_contribution, fill = signature )) +
  geom_col() +
  scale_fill_manual(values = colors) +
  theme_cowplot() +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none")  +
  labs(x = "Donor (Age in years)", y = "Relative signature contribution")

# Combine the relative mutational signature contributions across donors

figure_S2_bottom = donor_contri_blood + donor_contri_colon + donor_contri_lung + plot_layout(widths = c(1,3,1))
ggsave("manuscript/Supplementary_Figures/Figure_S2/figure_s2_bottom.svg", figure_S2_bottom, width = 16, height = 5)


