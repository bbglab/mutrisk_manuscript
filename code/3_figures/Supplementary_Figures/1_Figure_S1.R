# 1_Figure_S1.R
# Purpose: Generate supplementary figure S1 -- exonic mutation overview, whole-genome mutation rates, and trinucleotide context analysis.
library(readxl)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggpp)
library(ggh4x)
library(plyranges)
library(GenomeInfoDb)
source("code/0_functions/analysis_variables.R")

### define the plotting colors:
obs_pal <- c("#3ca951", "#4269d0", "#ff725c", "#6cc5b0", "#9498a0", "#97bbf5", "#9c6b4e",
             "#a463f2", "#efb118", "#ff8ab7")
# blood_colors = "#ff725c"
# lung_colors = c("#4269d0", "#7c86a1", "#161459")
# colon_colors = c("#3ca951", "#6cc5b0", "#145220", "#222e24")

# Load the Whole Genome mutations:
wgs_files = list.files("processed_data", pattern = "cell_muts.tsv.gz", recursive = TRUE, full.names = TRUE)
names(wgs_files) = gsub("_cell_muts.tsv.gz", "", basename(wgs_files))
list_muts_tissue = lapply(wgs_files, fread)
list_wgs_muts = lapply(list_muts_tissue, dplyr::select, sampleID, chr, pos, ref, alt, category) |>
  rbindlist(idcol = "tissue") |>
  filter(category != "chemotherapy")

Supp_tables = list()

wgs_muts = list_wgs_muts  |>
  group_by(tissue, sampleID) |>
  count() |>
  dplyr::rename(wgs_muts = n) |>
  ungroup()

# Data loading: Exonic mutations
dnds_files = list.files("processed_data/", pattern = "dnds", recursive = TRUE, full.names = TRUE) # double check - make sure to remove the duplex and blood mutations
dnds_files = dnds_files[!grepl("bladder|pancreas|liver|unique|comparison|chemotherapy", dnds_files)]
dnds_files = dnds_files[!grepl("exon", dnds_files)]
names(dnds_files) = paste0(str_split_i(dnds_files, "\\/", 3), "_", gsub("_dnds.rds", "", basename(dnds_files)))
list_muts = lapply(dnds_files, \(x) readRDS(x)[["annotmuts"]]) |>
  lapply(dplyr::select,  c(sampleID, chr, pos, ref, mut, txind, tx_gene, gene, strand, ref_cod, mut_cod)) |>
  rbindlist(idcol = "tissue_category", fill = TRUE) |>
  mutate(tissue = gsub("_.*", "", tissue_category),
         category = gsub(".*_", "", tissue_category)) |>
  filter(gene != "intronic")

sample_data = list_muts |>
  group_by(tissue) |>
  summarize(count = n_distinct(sampleID),
            mean_exonic_muts = dplyr::n()/count)

list_muts = left_join(list_muts, sample_data, by = "tissue")

sample_mut_counts = list_muts |>
  group_by(tissue, sampleID) |>
  count(name = "exome_muts")

# update the mutational count rates for whole-genome mutations:
# TODO: check if there are samples (for instance the cord blood samples - which have more than 1 mutation...)
sample_wgs_exome = left_join(wgs_muts, by = c("sampleID", "tissue"), sample_mut_counts) |>
  mutate(exome_muts = ifelse(is.na(exome_muts), 0, exome_muts)) |>  # fill in NA values with 0
  group_by(tissue) |>
  arrange(tissue, exome_muts) |>
  mutate(id = 1:dplyr::n()) |> ungroup()


# load all metadata:
metadata_files = c("processed_data/blood/blood_metadata.tsv", "processed_data/colon/colon_metadata.tsv",
                   "processed_data/lung/lung_metadata.tsv")
names(metadata_files) = gsub("_metadata.tsv", "", basename(metadata_files))
meta_age = lapply(metadata_files, \(x) fread(x) |>
                    dplyr::select(any_of(c("sampleID", "category", "age", "donor", "sensitivity", "coverage")))) |>
  rbindlist(idcol = "tissue", use.names = TRUE, fill = TRUE) |>
  distinct()
fwrite(meta_age, "processed_data/metadata_all.txt.gz")


Supp_tables$Supplementary_table_1 = sample_wgs_exome |>  left_join(meta_age) |> dplyr::select(-id)
Supp_tables$Supplementary_table_2 = meta_age |> dplyr::select(tissue, category, donor, age) |> distinct()

# Numbers for the manuscript: Take only the samples for which the mutations are also matching all filters
n_clones_n_donors = meta_age |>
  inner_join(sample_wgs_exome) |>
  group_by(tissue) |>
  summarize(number_of_clones = n_distinct(sampleID),
            number_of_donors = n_distinct(donor))

Supp_tables$Supplementary_table_4 = n_clones_n_donors

n_clones_n_donors = meta_age |>
  inner_join(sample_wgs_exome) |>
  group_by(tissue, category) |>
  summarize(number_of_clones = n_distinct(sampleID),
            number_of_donors = n_distinct(donor))
Supp_tables$Supplementary_table_5 = n_clones_n_donors

meta_age$age |> max()
meta_age$age |> min()

sample_wgs_age = inner_join(sample_wgs_exome, meta_age) |>
  filter(!grepl("Skin", tissue)) |>
  mutate(corrected_wgs = wgs_muts / sensitivity) # correct both the  WGS  muts for the relative sensitivity

count_category = sample_wgs_age |>
  filter(!is.na(age)) |>
  mutate(tissue_category = paste0(tissue, "_", category),
    category = case_when(category == "non-smoker" & grepl("Lung", tissue) ~ "normal",
                              category == "colon_normal" ~ "normal",
                              category == "colon_IBD" ~ "IBD",
                              tissue == "skin_keratinocyte" ~ "keratinocyte",
                              tissue == "skin_melanocyte" ~ "melanocyte",
                              tissue == "skin_fibroblast" ~ "fibroblast",
                              grepl("colon_POLE", tissue) ~ "POLE",
                              grepl("colon_POLD1", tissue) ~ "POLD1",
                              tissue == "blood_normal" ~ "normal",
                              category == "NORMAL" ~ "normal",
                              category == "ARLD" ~ "ARLD",
                              category == "NAFLD" ~ "NAFLD",
                              .default = category),
         tissue = case_when(grepl("colon", tissue) ~ "Colon",
                            grepl("skin", tissue) ~ "Skin",
                            grepl("lung", tissue) ~ "Lung",
                            grepl("blood", tissue) ~ "Blood",
                            .default = tissue),
         tissue_category2 = paste(tissue, category)) |>
  mutate(tissue = factor(tissue, levels = c("Blood", "Lung", "Colon", "Skin", "liver")),
         category = factor(category, levels = c("normal", "IBD", "POLD1", "POLE", "fibroblast",
                                                "keratinocyte", "melanocyte", "non-smoker", "ex-smoker",
                                                "smoker", "NAFLD", "ARLD"))) |>
  group_by(tissue_category2, tissue, category) |>
  mutate(corrected_exome = exome_muts / sensitivity) |>
  arrange(corrected_exome) |>
  mutate(id = 1:dplyr::n())

means = count_category |>
  group_by(tissue_category, tissue, category) |>
  summarize(n_samples = dplyr::n(),
            mean = round(mean(corrected_exome),1)) |>
  mutate(label = paste0("n = ", n_samples))

exome_counts_plot = ggplot(count_category, aes(x = id, y = corrected_exome, group = tissue_category)) +
  geom_point(aes(color = tissue_category)) +
  geom_text_npc(data = means, aes(npcx = 0.5, npcy = 0.96, label = label), hjust = 0.5) +
  geom_segment(data = means, aes(x = n_samples/4, xend = (n_samples/4)*3, y = mean), color = "black") +
  geom_text(data = means, aes(x = n_samples/4, y = mean, label = round(mean, 1)), color = "black", vjust = -0.5, hjust = 0.2) +
  facet_nested(~ tissue + category,  scales = "free_x") +
  theme_bw() +
  scale_y_log10(guide = "axis_logticks") +
  scale_color_manual(values = tissue_category_colors[-2]) +
  labs(x = NULL, y = "Exonic mutations/cell") +
  theme(panel.grid = element_blank(),  axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        strip.background = element_blank(), legend.position = "none",
        ggh4x.facet.nestline = element_line())
exome_counts_plot
ggsave("manuscript/Supplementary_Figures/Figure_S1/count_genome_category.png", exome_counts_plot, width = 12, height = 4)
ggsave("manuscript/Supplementary_Figures/Figure_S1/supplementary_figure_1.png", exome_counts_plot, width = 12, height = 4)
# expand the graph: Make level for category, cell type and potential category

#### Whole genome analysis ####

# mean number of mutations across samples
sample_wgs_age |>
  group_by(tissue, category) |>
  summarize(mean_wgs = mean(corrected_wgs),
            sd_wgs = sd(corrected_wgs))


library(lme4)
library(lmerTest)
# consider removing the intercept of 0 model
# run with an intercept of 0 (at the moment not used in the study)
#run_lmer_0 = function(data) { lmer(corrected_wgs ~ 0 + age  +  (0 + age|donor), data = data)}

# run with
run_lmer = function(data) { lme4::lmer(corrected_wgs~ age  + (1|donor), data = data) }
run_lmertest = function(data) { lmerTest::lmer(corrected_wgs~ age  +  sensitivity + (1|donor), data = data) } # lmertest function which returns p-values
# Run linear-mixed-effects models for the different tissues. Set the effect of patients as "random" as we perform repeated sampling
nonexposed_data = sample_wgs_age |>
  group_by(tissue) |>
  filter(!grepl("POLE|POLD1|IBD", tissue)) |>
  filter(category %in% c("-", "non-smoker", "normal", "NORMAL"))

# show summaries of the modeling data
model = run_lmertest(nonexposed_data |> filter(tissue == "blood"))
summary(model)
model = run_lmertest(nonexposed_data |> filter(tissue == "lung"))
summary(model)
model = run_lmertest(nonexposed_data |> filter(tissue == "colon"))
summary(model)

# Use the lmer function from the lme4 package
models = nonexposed_data |>
  group_by(tissue) |>
  group_map(~ as.data.frame(summary(run_lmer(.))[["coefficients"]]))
names(models) = group_keys(nonexposed_data)[[1]]
models = lapply(models, \(x) x[1:2,])
models_nonexposed = rbindlist(models, idcol = "tissue") |>
  mutate(type = rep(c("intercept", "slope"), length(models))) |>
  dplyr::select(tissue, Estimate, type) |>
  pivot_wider(names_from = type, values_from = Estimate)

#plot nonexposed tissue:
nonexposed_colors = c(blood_colors[1], colon_colors[1],  lung_colors[1])
nonexposed_plot = nonexposed_data |>
  mutate(tissue = gsub("_.*", "", tissue)) |>
  ggplot(aes(x = age, y = corrected_wgs, fill = tissue, color = tissue, group = tissue)) +
  annotate("rect", xmin = 0, xmax = 85, ymin = 0, ymax = 5000, alpha = 0, color= "grey30", linetype = "dashed") +
  geom_point(shape = 21,  stroke = .3, size = 2, color = "white")  +
  geom_abline(data = models_nonexposed,
              aes(slope = slope, intercept = intercept, color = tissue)) +
  theme_cowplot() +
  scale_color_manual(values = nonexposed_colors) +
  scale_fill_manual(values = nonexposed_colors) +
  guides(color = "none") +
  scale_y_continuous(labels = scales::comma) +
  theme(legend.position = "inside", legend.position.inside = c(0.1, 0.7), plot.title =element_text(hjust = 0.5), legend.background = element_blank(),
        legend.key = element_rect(fill = NA)) +
  labs(y = "number of whole genome mutations", title = "unexposed tissues", fill = NULL, x = "Age (years)")
nonexposed_plot

# check the difference between corrected for coverage and non corrected:
nonexposed_plot_comparison = nonexposed_data |>
  mutate(tissue = gsub("_.*", "", tissue)) |>
  ggplot(aes(x = age, y= corrected_wgs, fill = tissue, color = tissue, group = donor)) +
  annotate("rect", xmin = 0, xmax = 85, ymin = 0, ymax = 5000, alpha = 0, color= "grey30", linetype = "dashed") +
  geom_point(aes(x = age + 1, y = wgs_muts), alpha = 0.5, color = "black") +
  geom_point(shape = 21,  stroke = .3, size = 2, color = "white")  +
  geom_abline(data = models_nonexposed,
              aes(slope = slope, intercept = intercept, color = tissue)) +
  theme_cowplot() +
  scale_color_manual(values = nonexposed_colors) +
  scale_fill_manual(values = nonexposed_colors) +
  guides(color = "none") +
  facet_grid(. ~ tissue) +
  scale_y_continuous(labels = scales::comma) +
  theme(legend.position = "inside", legend.position.inside = c(0.1, 0.7), plot.title =element_text(hjust = 0.5), legend.background = element_blank(),
        legend.key = element_rect(fill = NA)) +
  labs(y = "number of SNVs in the whole genome", title = "unexposed tissues", fill = NULL, x = "Age (years)")
nonexposed_plot_comparison
# conclusion there is a small effect of the correction for the individual samples. Biggest differences are for some of the colon samples.

# plot to compare the effect of correcting for coverage
data = nonexposed_data |>
  pivot_longer(cols = c(wgs_muts, corrected_wgs), values_to = "muts", names_to = "type")

data_variance = data |>
  group_by(tissue, donor, type) |>
  summarize(IQR = IQR(muts),
            sd = sd(muts))

# Plot mutation accumulation in the lung:
lung_data = sample_wgs_age |>
  filter(grepl("lung", tissue))

exposed_lung_plot = lung_data |>
  mutate(category = factor(category, levels = c("non-smoker", "ex-smoker", "smoker"))) |>
  ggplot(aes(x = age, y = corrected_wgs, fill = category, group = category, color = )) +
  geom_point(shape = 21, color = "white", stroke = 0.3, size = 2)  +
  annotate("rect", xmin = 0, xmax = 85, ymin = 0, ymax = 5000, alpha = 0, color = "grey30", linetype = "dashed") +
  geom_abline(data = models_nonexposed |> filter(grepl("lung", tissue)),
              aes(slope = slope, intercept = intercept, color = lung_colors[[1]])) +
  theme_cowplot() +
  scale_color_manual(values = lung_colors) +
  scale_fill_manual(values = lung_colors) +
  scale_y_continuous(labels = scales::comma) +
  labs(y = NULL, fill = NULL, title = "Lung", x = "Age (years)") +
  theme(legend.position = "inside", legend.position.inside = c(0.1, 0.7), plot.title = element_text(hjust = 0.5)) +
  guides(color = "none")
exposed_lung_plot

# Colon
colon_data = sample_wgs_age |>
  filter(grepl("colon", tissue)) |>
  mutate(category = case_when(grepl("POLD1", category) ~ "POLD1",
                              grepl("POLE", category) ~ "POLE",
                              .default = category),
         category = factor(category, levels = c("normal", "IBD", "POLD1", "POLE")))

exposed_colon_plot = colon_data |>
  ggplot(aes(x = age, y = corrected_wgs, fill = category, group = category, color = category)) +
  geom_abline(data = models_nonexposed |> filter(grepl("colon", tissue)),
              aes(slope = slope, intercept = intercept, color = colon_colors[1])) +
  geom_point(shape = 21, color = "white", stroke = 0.3, size = 2)  +
  annotate("rect", xmin = 0, xmax = 85, ymin = 0, ymax = 5000, alpha = 0, color = "grey30", linetype = "dashed") +
  theme_cowplot() +
  scale_color_manual(values = colon_colors) +
  scale_fill_manual(values = colon_colors) +
  scale_y_continuous(labels = scales::comma) +
  theme(legend.position = "inside", legend.position.inside = c(0.1, 0.7), plot.title = element_text(hjust = 0.5)) +
  labs(y = NULL, fill  = NULL, title = "Colon", x = "Age (years)") +
  guides(color = "none")
exposed_colon_plot

aging_rate_plot = nonexposed_plot + exposed_lung_plot + exposed_colon_plot
ggsave("manuscript/Supplementary_Figures/Figure_S1/wgs_rate_tissues.png", aging_rate_plot,  width = 12, height = 4.5, dpi = 600)
ggsave("manuscript/Supplementary_Figures/Figure_S1/wgs_rate_tissues.pdf", aging_rate_plot,  width = 11, height = 4.2, dpi = 600)

FS1A = prep_plot(aging_rate_plot, "A")
FS1B = prep_plot(exome_counts_plot, "B")

Figure_S1 = FS1A / FS1B
ggsave("manuscript/Supplementary_Figures/Figure_S1/Figure_S1.png", Figure_S1,  width = 15, height = 10, dpi = 600)
ggsave("manuscript/Supplementary_Figures/Figure_S1/Figure_S1.pdf", Figure_S1,  width = 15, height = 10)


# figure S5B - comparison of colon normal and POLD1 mutated cells
figure_S5B = colon_data |>
  filter(category %in% c("normal", "POLD1")) |>
  ggplot(aes(x = age, y = corrected_wgs, fill = category, group = category, color = category)) +
  geom_abline(data = models_nonexposed |> filter(grepl("colon", tissue)),
              aes(slope = slope, intercept = intercept, color = colon_colors[1])) +
  geom_point(shape = 21, color = "white", stroke = 0.3, size = 2)  +
  annotate("rect", xmin = 0, xmax = 85, ymin = 0, ymax = 5000, alpha = 0, color = "grey30", linetype = "dashed") +
  theme_cowplot() +
  scale_color_manual(values = colon_colors[c(1,3)]) +
  scale_fill_manual(values = colon_colors[c(1,3)]) +
  scale_y_continuous(labels = scales::comma) +
  theme(legend.position = "inside", legend.position.inside = c(0.1, 0.7), plot.title = element_text(hjust = 0.5)) +
  labs(y = NULL, fill  = NULL, title = "Colon", x = "Age (years)") +
  guides(color = "none")
figure_S5B
saveRDS(figure_S5B, "manuscript/figure_panels/figure_s5/Figure_S5B.rds")

# TODO - check the plot below - see if it is possible to make it with the lmm-based approaches for correctness.
# make a simple plot for the biorender figure:
aging_rate_plot_simple = nonexposed_data |>
  mutate(tissue = gsub("\n.*", "", tissue)) |>
  ggplot(aes(x = age, y = corrected_wgs, fill = tissue, group = tissue)) +
  ggpmisc::stat_poly_line(aes(color = tissue), fill = "white") +
  theme_cowplot() +
  scale_color_manual(values = nonexposed_colors) +
  guides(fill = "none") +
  labs(y = "mutations", color = NULL) +
  theme(legend.position = "inside", legend.position.inside = c(0.25, 0.8),
        axis.text = element_blank())
ggsave("plots/aging_rate_simple.png", width = 2.2, height = 1.8)

# Introduction slide - 96-trinucleotide context
# split the lung data into smokers and nonsmokers:
lung_muts = meta_age |>
  filter(tissue == "lung") |>
  filter(category %in% c("smoker", "non-smoker"))  |>
  inner_join(list_muts_tissue$lung, by = c("sampleID", "category")) |>
  mutate(sampleID = category) |>
  dplyr::select(sampleID, chr, pos, ref, alt, status = category) |>
  group_by(sampleID)

lung_list = split(lung_muts, lung_muts$sampleID)
non_smoker_context = get_mut_context(lung_list$`non-smoker`)
smoker_context = get_mut_context(lung_list$smoker)
colon_context = get_mut_context(list_muts_tissue$colon) |>
  mutate(sampleID = "Colon")
blood_context = get_mut_context(list_muts_tissue$blood) |>
  mutate(sampleID = "Blood")

muts_list = rbindlist(list(colon_context, smoker_context, non_smoker_context, blood_context), fill = TRUE)
mm = muts_to_mat(muts_list)
colnames(mm)[3:4]  = c("Lung\nnon-smoker", "Lung\nsmoker")

# Plot the 96-trinucleotide profile of smoking:
profiles_tissues = plot_96_profile2(mm, relative = F, free_y = TRUE, horizontal_labels = TRUE) +
  theme(axis.text = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank())
ggsave("plots/trinuc_tissue.png", profiles_tissues, width = 6.5, height = 5)
ggsave("plots/trinuc_tissue.pdf", profiles_tissues, width = 6.5, height = 5)

# check the signature contributions of smoking samples:
signatures = read.delim("raw_data/resources/COSMIC_v3.4_SBS_GRCh37.txt") |> column_to_rownames("Type") |> as.matrix()
signatures = signatures[TRIPLETS_96,]
sig_contribution = MutationalPatterns::fit_to_signatures(mm, as.matrix(signatures))$contribution
rel_contri = prop.table(sig_contribution[,4]) * 100
contri = as.data.frame( t(t(signatures) * sig_contribution[,"Lung\nsmoker"]))
contribution = data.frame(smoking = contri$SBS4 + contri$SBS92 + contri$SBS29,
                          aging = contri$SBS1 + contri$SBS5 + contri$SBS40a ,
                          other = rowSums(contri |> dplyr::select(-c(SBS1, SBS4, SBS5, SBS92, SBS29, SBS40a))))

contri = data.frame(sample = c("cigarette\nsmoke\n", "aging\n", "other\n"),
                    substitution = "C>A",
                    contribution = proportions(colSums(contribution)) * 100) |>
  mutate(contribution = paste0(round(contribution, 1), "%"))

colnames(contribution) = paste(contri$sample, contri$contribution )

# Plot the 96-trinucleotide profile of smoking:
profiles_signatures = plot_96_profile2(contribution, relative = F, free_y = TRUE, horizontal_labels = TRUE) +
  theme(axis.text = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank())
ggsave("plots/trinuc_contribution.png", profiles_signatures, width = 5.4, height = 3.6)

ggsave("plots/trinuc_contirbution.pdf", profiles_signatures, width = 6.5, height = 5)


# plot the contribution in the 6 different mutation types
colnames(contribution) = c("cigarette\nsmoke", "aging", "other")
contri_long = contribution |>
  rownames_to_column("triplet") |>
  left_join(triplet_match_substmodel) |>
  group_by(type) |>
  summarize(across(c(`cigarette\nsmoke`, aging, other), sum)) |>
  pivot_longer(-type, names_to = "signature", values_to = "number of mutations")

ggplot(contri_long, aes(x = type, y = `number of mutations`, fill = type)) +
  geom_col() +
  ggh4x::facet_grid2(. ~ type, scales = "free",  axes = "x",remove_labels = "x",
                     strip = ggh4x::strip_themed(background_x = list(element_rect(fill = "#2EBAED"),
                                                                     element_rect(fill = "#000000"),
                                                                     element_rect(fill = "#DE1C14"),
                                                                     element_rect(fill = "#D4D2D2"),
                                                                     element_rect(fill = "#ADCC54"),
                                                                     element_rect(fill = "#F0D0CE")))) +
  theme_classic() +
  scale_fill_manual(values = COLORS6) +
  scale_y_continuous(expand=expansion(mult=c(0,0.1))) +
  guides(fill = "none") +
  theme_classic() +
  theme(axis.text= element_blank(),
        panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_text(colour = "white", size = 12),
        strip.background = element_rect(colour = "white"),
        axis.ticks = element_blank(), axis.line = element_line(linewidth = 0.3)) +
  labs(x = NULL)

# signature contributions
ggplot(contri_long, aes(x = type, y = `number of mutations`, fill = signature)) +
  geom_col() +
  ggh4x::facet_grid2(. ~ type, scales = "free",  axes = "x",remove_labels = "x",
                     strip = ggh4x::strip_themed(background_x = list(element_rect(fill = "#2EBAED"),
                                                                     element_rect(fill = "#000000"),
                                                                     element_rect(fill = "#DE1C14"),
                                                                     element_rect(fill = "#D4D2D2"),
                                                                     element_rect(fill = "#ADCC54"),
                                                                     element_rect(fill = "#F0D0CE")))) +
  theme_classic() +
  ggsci::scale_fill_igv() +
  scale_y_continuous(expand=expansion(mult=c(0,0.1))) +
  theme_classic() +
  theme(legend.position = "inside", legend.position.inside = c(1.2, 0.65),
        axis.text = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_text(colour = "white", size = 12),
        strip.background = element_rect(colour = "white"),
        axis.ticks = element_blank(), axis.line = element_line(linewidth = 0.3),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        plot.margin = unit(c(1,100,1,1), "pt")) +
  labs(x = NULL, fill = NULL) +
  guides(shape = guide_legend(override.aes = list(size = 2)),
         color = guide_legend(override.aes = list(size = 2)))

# TODO: Compare smoker vs non-smoker input data profiles
# large-scale-genome context plot
get_mb_rates = function(context_file, genome_bed) {

  gr = context_file |>
    dplyr::select(chr, pos) |>
    `colnames<-`(c("seqnames","start")) |>
    mutate(end = start) |>
    as_granges() |>
    `seqlevelsStyle<-`("UCSC")

  bin_overlaps = findOverlaps(gr, genome_bed)

  genome_df = as.data.frame(genome_bed)
  genome_df$counts = tabulate(subjectHits(bin_overlaps), nbins =  nrow(genome_df))
  return(genome_df)
}

seqinfo <- rtracklayer::SeqinfoForUCSCGenome("hg19")
seqinfo <- keepStandardChromosomes(seqinfo)
granges <- tileGenome(seqinfo, tilewidth = 1e5, cut.last.tile.in.chrom = T)

rates_lung_smoker = get_mb_rates(ungroup(smoker_context), granges)
rates_lung_nonsmoker = get_mb_rates(ungroup(non_smoker_context), granges)
rates_colon = get_mb_rates(colon_context, granges)
rates_blood = get_mb_rates(blood_context, granges)
df_rates = tibble(chr = as.character(seqnames(granges)),
                  bin = 1:length(granges),
                  `Lung non-smoker` = rates_lung_nonsmoker$counts,
                  `Lung smoker` = rates_lung_smoker$counts,
                  Colon =rates_colon$counts,
                  Blood =rates_blood$counts)

df_rate = df_rates |>
  filter(bin < 1.5e3) |>
  mutate(across(-c(chr, bin), \(x) x/sum(x))) |>
  filter(chr == "chr1") |>
  pivot_longer(-c(chr, bin))

label_data = df_rate |>
  dplyr::select(name) |> distinct()

genomic_rate_plot = df_rate |>
  ggplot( aes(x = bin, y = value, color = name)) +
  geom_line() +
  facet_grid(name ~ . , scales = "free_y",  axes = "all", axis.labels = "margins") +
  cowplot::theme_cowplot() +
  geom_text_npc(data = label_data, aes(label = name, npcx = 0.05, npcy = 1.3), size = 4.5) +
  scale_color_manual(values = c(blood_colors, colon_colors[1], lung_colors[c(1,3)])) +
  theme(legend.position = "none", axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
        strip.text = element_blank()) +
  labs(x = "Chromosome 1 - 100kb bins", y = NULL, color = NULL)
genomic_rate_plot
ggsave("manuscript/Schematic_poster_presentations/100Kb_genomic_rate.png", genomic_rate_plot, width = 7, height = 4, dpi = 600)


# Save supplementary tables 1-3
# also update this for the exposed conditions:
openxlsx::write.xlsx(Supp_tables, file = "manuscript/Supplementary_Tables/Supplementary_Tables_1-5.xlsx")