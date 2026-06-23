# Preprocess IntOGen mutation data for KRAS across colon/lung/blood cohorts.
#
# Reads the pan-cancer IntOGen mutation dump and the cohort metadata table,
# produces a pivoted table with one row per amino-acid position and one
# column per substitution class, values = percentage of cohort samples
# carrying that mutation.
#
# Intended for the KRAS all-mutation mirror plot (left column of the
# composite figure), where IntOGen observations sit above MutRisk predictions.
#
# Run from the mutrisk_manuscript project root:
#   Rscript code/2_resources/6_intogen_kras_distribution.R
#
# Inputs:
#   raw_data/intogen/mutations.tsv   — pan-cancer mutation dump
#   raw_data/intogen/cohorts.tsv     — COHORT → CANCER_TYPE mapping + SAMPLES
#
# Output:
#   processed_data/intogen/kras_distribution.tsv.gz

library(data.table)
library(tidyverse)

source("code/0_functions/analysis_variables.R")

# ---- 1. Load data ----
cat("Loading IntOGen data …\n")
intogen  = fread("raw_data/intogen/mutations.tsv")
cohorts  = fread("raw_data/intogen/cohorts.tsv")

# ---- 2. Join cohort metadata (CANCER_TYPE, platform, total SAMPLES) ----
intogen = merge(intogen, cohorts[, .(COHORT, CANCER_TYPE, PLATFORM, SAMPLES)],
                by = "COHORT", all.x = TRUE)

# ---- 3. Filter to KRAS in the three target cancer types ----
target_cancers = c("COADREAD", "LUSC", "AML")

kras = intogen[SYMBOL == "KRAS" & CANCER_TYPE %in% target_cancers]

if (nrow(kras) == 0) {
  stop(
    "No KRAS rows found for CANCER_TYPE ", paste(target_cancers, collapse = "/"),
    ".\n  Available CANCER_TYPE values in file: ",
    paste(sort(unique(intogen$CANCER_TYPE)), collapse = ", "),
    "\n  Available symbols matching 'KRAS': ",
    paste(unique(intogen[grepl("KRAS", SYMBOL, ignore.case = TRUE)]$SYMBOL), collapse = ", "),
    call. = FALSE
  )
}

cat(sprintf("  Found %d KRAS rows across cancer types %s\n",
            nrow(kras), paste(unique(kras$CANCER_TYPE), collapse = ", ")))

# ---- 4. Remove rows without a valid amino-acid position ----
# Keep only positions that parse cleanly as a single integer
kras = kras[grepl("^[0-9]+$", Protein_position)]
kras[, Protein_position := as.integer(Protein_position)]
setnames(kras, "Protein_position", "position")

# Drop multi-nucleotide mutations (REF or ALT length > 1)
kras = kras[nchar(REF) == 1 & nchar(ALT) == 1]

cat(sprintf("  After filtering: %d rows\n", nrow(kras)))

# ---- 5. Map tissue name from CANCER_TYPE ----
kras[, tissue := fcase(
  CANCER_TYPE == "COADREAD", "colon",
  CANCER_TYPE == "LUSC",     "lung",
  CANCER_TYPE == "AML",      "blood"
)]

# ---- 6. Map REF>ALT to the six substitution classes (COLORS6) ----
kras[, type := paste0(REF, ">", ALT)]
kras[, type := fcase(
  type %in% c("G>T", "C>A"), "C>A",
  type %in% c("G>C", "C>G"), "C>G",
  type %in% c("G>A", "C>T"), "C>T",
  type %in% c("A>T", "T>A"), "T>A",
  type %in% c("A>G", "T>C"), "T>C",
  type %in% c("A>C", "T>G"), "T>G"
)]

# Report any mutations that couldn't be mapped
unmapped = unique(kras[is.na(type), .(REF, ALT)])
if (nrow(unmapped) > 0) {
  cat("  WARNING: Unmapped REF>ALT pairs (these rows will be dropped):\n")
  print(unmapped)
}
kras = kras[!is.na(type)]

# ---- 7. Get cohort-level sample denominator ----
# SAMPLES from cohorts.tsv = total samples in the cohort.
# Summing across all cohorts within a cancer type gives the denominator.
cohort_totals = cohorts[CANCER_TYPE %in% target_cancers,
                         .(total_samples = sum(SAMPLES)),
                         by = CANCER_TYPE]

cat("  Cancer type sample denominators:\n")
print(cohort_totals)

# ---- 8. Aggregate: sum SAMPLES per (tissue, cancer_type, position, class) ----
# SAMPLES.x from intogen = number of samples with this specific mutation
kras_agg = kras[, .(samples = sum(SAMPLES.x)),
                 by = .(tissue, CANCER_TYPE, position, type)]

# ---- 9. Join cohort totals and compute percentage ----
kras_agg = merge(kras_agg, cohort_totals, by = "CANCER_TYPE")
kras_agg[, samples_percent := 100 * samples / total_samples]

# ---- 10. Pivot to wide format: position × substitution class ----
intogen_pivot = dcast(
  kras_agg,
  tissue + CANCER_TYPE + position ~ type,
  value.var = "samples_percent",
  fill = 0
)

# Ensure all six substitution classes are present (even if 0 mutations)
for (cls in names(COLORS6)) {
  if (!cls %in% names(intogen_pivot)) {
    intogen_pivot[, (cls) := 0]
  }
}
setcolorder(intogen_pivot, c("tissue", "CANCER_TYPE", "position", names(COLORS6)))

# ---- 11. Save ----
dir.create("processed_data/intogen", showWarnings = FALSE, recursive = TRUE)
fwrite(intogen_pivot, "processed_data/intogen/kras_distribution.tsv.gz")

cat(sprintf("Wrote processed_data/intogen/kras_distribution.tsv.gz (%d rows × %d cols)\n",
            nrow(intogen_pivot), ncol(intogen_pivot)))
