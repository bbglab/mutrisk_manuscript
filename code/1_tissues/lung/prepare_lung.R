### Processe Lung data
# 1. Prepare the script to match the prerequisites:
# metadata: Data frame with 4 columns: "sampleID", "category", "age" and "donor" and "sensitivity"
# cell_muts: Data frame with 6 columns: "sampleID", "chr", "pos", "ref", "alt", "category", "donor"
library(dndscv)
library(MutationalPatterns)
library(readxl)
source("code/0_functions/analysis_variables.R")

# load the dndscv RefCDS
tissue = "lung"

outdir = paste0("processed_data/", tissue, "/")

metadata = read.delim("raw_data/lung/Lung_organoids_telomeres_with_contamination_20190408.txt") |>
  dplyr::select(Sample, Smoking, Age, Patient, seq.X) |>
  dplyr::rename(sampleID = Sample, category = Smoking, age = Age, donor = Patient, coverage = seq.X)

# function sourced from: https://stackoverflow.com/questions/12945687/read-all-worksheets-in-an-excel-workbook-into-an-r-list-with-data-frames
read_excel_allsheets <- function(filename) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X, col_types = "text"))
  names(x) <- sheets
  return(x)
}

input_sbs_list = read_excel_allsheets("raw_data/lung/Bronchial epithelium subs.xlsx") |>
  rbindlist(idcol = "donor") |>
  mutate(across(c(Pos, DEP, MTR, VAF), as.numeric))
input_indel_list = read_excel("raw_data/lung/Bronchial epithelium indels.xlsx") |>
  mutate(donor = str_sub(Sample, 1,7)) |>
  dplyr::select(all_of(colnames(input_sbs_list)))
input_muts = bind_rows(input_sbs_list, input_indel_list) |>
  arrange(donor) |>
  dplyr::rename(sampleID = Sample, chr = Chrom, pos = Pos, ref = Ref, alt = Alt) |>
  mutate(across(c(pos, DEP, MTR, VAF), as.numeric))

# plot the variant allele frequencies of a single sample to determine the clonality (and quality) of the input data
VAF_lung_samples = input_muts |>
  filter(grepl("PD26988", sampleID)) |>
  ggplot(aes(x = VAF, fill = sampleID)) + ggridges::geom_density_ridges_gradient(aes(y = sampleID)) +
  theme(legend.position = "none")
ggsave("plots/lung/VAF_plot_PD26988.png", VAF_lung_samples, height = 8, width = 4, dpi = 400)

# get the underlying VAFS for each sample:
vaf_estimates = input_muts |>
  group_by(sampleID, donor) |>
  summarize(vaf_estimate = estimate_vaf(VAF))

# update the metadata file with the sensitivity data for the different individuals
metadata = metadata |>
  full_join(vaf_estimates) |>
  filter(!is.na(category)) |>
  mutate(sensitivity = get_sensitivity(coverage, vaf_estimate, minalt = 2),
         category = factor(category, levels = c("non-smoker", "ex-smoker", "smoker"))) |>
  select(sampleID, category, donor, age, sensitivity, coverage)

select_sigs = c("SBS1", "SBS4", "SBS5", "SBS2", "SBS13", "SBS92", "SBS16") # check the signatures present in the lung paper

# join the mutations with the expected mutation data
cell_muts = inner_join(input_muts, metadata) |>
  dplyr::select(sampleID, chr, pos, ref, alt, category, donor, VAF) |>
  filter(nchar(ref) == 1 & nchar(alt) == 1)


vaf_overview = create_vaf_overview(cell_muts = cell_muts |> dplyr::rename(vaf = VAF), sample_names = c("PD26988", "PD37453"))
ggsave("manuscript/Supplementary_notes/Supplementary_Note_1/figure_lung.png", vaf_overview, width = 10, height = 10)

# calculate the estimated coverage for each sample
list_results = list()
for (i in unique(metadata$category)) {
  list_results[[i]] = effect_coverage_vaf(cell_muts |> filter(category %in% i),
                                          metadata = metadata |> filter(category %in% i))
  list_results[[i]]$plot = list_results[[i]]$plot + ggtitle(i)
}

list_plots = lapply(list_results, \(x) x[["plot"]])
ggsave(paste0("manuscript/Supplementary_notes/Supplementary_Note_1/", tissue, "_coverage_var_correction.png"),
       wrap_plots(list_plots, ncol = 2), width = 12, height = 7.5)


# save the data:
fwrite(metadata, paste0(outdir, tissue, "_metadata.tsv"))
fwrite(cell_muts, file = paste0("processed_data/", tissue, "/", tissue, "_cell_muts.tsv.gz"))

output_path = paste0("processed_data/", tissue, "/")
if (!dir.exists(output_path)) {dir.create(output_path)}

input_signatures = c("SBS1", "SBS4", "SBS5", "SBS2", "SBS13", "SBS92", "SBS16") # check the signatures present in the lung paper

list_results = list()
for (i in unique(metadata$category)) {

  mutrisk_pipeline(output_path = output_path,
                   cell_muts = cell_muts |> filter(category == i),
                   metadata = metadata |> filter(category == i),
                   name = i, input_signatures = input_signatures,
                   multiple_refit_methods = TRUE,
                   sensitivity_correction =  TRUE)
}

# summarize the results from the tissue-specific analysis (see following code from and adapt:)
sig_rate_files = list.files(paste0("processed_data/", tissue), pattern = 'sig_rate_per_sample',
                            full.names = TRUE, recursive = TRUE)
rates = lapply(sig_rate_files, fread)
names(rates) = gsub("_sig_rate_per_sample.tsv.gz", "", basename(sig_rate_files))
sig_donor_rates = rbindlist(rates, use.names = TRUE, fill = TRUE) |>
  mutate(category = factor(category, levels = levels(metadata$category)))  |>
  inner_join(metadata) |>
  group_by(donor, signature) |>
  summarize(across(contains(">"), mean)) |>
  pivot_longer(contains(">"))
fwrite(sig_donor_rates, file = paste0("processed_data/", tissue, "/", tissue, "_sig_donor_rates.tsv.gz"))

# load the mutation rates
mrates_files = list.files(paste0("processed_data/", tissue), pattern = "rate_per_sample.tsv.gz", full.names = TRUE, recursive = TRUE)
mrates_files = mrates_files[!grepl("sig", mrates_files)] # exclude the signature-specific variants
rates = lapply(mrates_files, fread)
names(rates) = gsub("_rate_per_sample.tsv.gz", "", basename(mrates_files))
expected_rates = rbindlist(rates) |>
  select(-sensitivity, -coverage) |>
  mutate(category = factor(category, levels = levels(metadata$category)))  |>
  inner_join(metadata) |>
  dplyr::select(-c(donor, age, sensitivity))
fwrite(expected_rates, file = paste0("processed_data/", tissue, "/", tissue, "_expected_rates.tsv.gz"))

# load the relative mutation ratio for each sample for the APC gene:
dnds_files = list.files(paste0("processed_data/", tissue), pattern = "^.*dnds.rds", full.names = TRUE, recursive = TRUE)
dnds_files = dnds_files[!grepl("unique|exon", dnds_files)]
names(dnds_files) = names(rates)
ratios = lapply(dnds_files, \(x) {
  readRDS(x)$genemuts |>
    mutate(ratio = exp_syn_cv / exp_syn) |>
    dplyr::select(gene_name, ratio)}) |>
  rbindlist(idcol = "category") |>
  mutate(category = factor(category, levels = levels(metadata$category)))
fwrite(ratios, file = paste0("processed_data/", tissue, "/", tissue, "_mut_ratios.tsv.gz"))