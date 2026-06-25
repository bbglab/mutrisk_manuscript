# Script to classify individual mutations in cancer genes as drivers
# Definition of cancer driver mutation: Being classified in BoostDM predictions and in either COSMIC/GENIE mutation databases.
library(data.table)
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(patchwork)
library(cowplot)
library(eulerr)

output_folder = "processed_data/boostdm/boostdm_genie_cosmic/"

nucs = c("C", "A", "G", "T")
# Load data
# GENIE
genie_data = fread("raw_data/GENIE_17/data_mutations_extended.txt")
genie_data = genie_data |>
  filter(Reference_Allele %in% nucs & Tumor_Seq_Allele2 %in% nucs) |>
  filter(Consequence %in% c("missense_variant", "synonymous_variant", "stop_gained", "start_lost", "stop_lost")) |>
  mutate(aachange = sub("p.", "", HGVSp_Short),
    type = paste0(Reference_Allele, ">", Tumor_Seq_Allele2),
    type = case_when(type == "A>T" ~ "T>A", type == "A>C" ~ "T>G",
                     type == "A>G" ~ "T>C", type == "G>C" ~ "C>G",
                     type == "G>A" ~ "C>T", type == "G>T" ~ "C>A",
                     .default = type)) |>
  select(Hugo_Symbol, Protein_position, aachange, Tumor_Sample_Barcode, type) |>
  dplyr::rename(gene_name = Hugo_Symbol, position = Protein_position, SAMPLE_ID = Tumor_Sample_Barcode)

# merge the data to obtain sample specificity
genie_sample_info = fread("raw_data/GENIE_17/GENIE_release_17.0/data_clinical_sample.txt", skip = 4) |>
  select(SAMPLE_ID, CANCER_TYPE, ONCOTREE_CODE)
genie_cancer_type = inner_join(genie_data, genie_sample_info) |>
  arrange(ONCOTREE_CODE, gene_name) |>
  mutate(genie_present = TRUE)
fwrite(genie_cancer_type, "processed_data/GENIE_17/GENIE_17_processed.txt.gz")
genie_cancer_type = genie_cancer_type |> select(-type, -SAMPLE_ID) |> distinct()

# COSMIC mutations
COSMIC_raw = fread("raw_data/COSMIC/Cosmic_MutantCensus_Tsv_v101_GRCh37/Cosmic_MutantCensus_v101_GRCh37.tsv.gz")
COSMIC = COSMIC_raw |>
  filter(grepl("missense_variant|stop_gained|frameshift_variant", MUTATION_DESCRIPTION) &
           nchar(GENOMIC_MUT_ALLELE) ==1) |>
  mutate(aachange = sub("p.", "", MUTATION_AA),
         position = as.numeric(parse_number(gsub("\\.", "", MUTATION_AA))),
         pos = GENOME_START, alt = GENOMIC_MUT_ALLELE) |>
  select(GENE_SYMBOL, position, aachange, COSMIC_PHENOTYPE_ID) |>
  filter(!aachange %in% "?")  |>
  dplyr::rename(gene_name = GENE_SYMBOL)

COSMIC_phenotype_link = fread("raw_data/COSMIC/Cosmic_Classification_Tsv_v101_GRCh37/Cosmic_Classification_v101_GRCh37.tsv.gz") |>
  select(COSMIC_PHENOTYPE_ID, PRIMARY_SITE)
COSMIC_cancer_type = left_join(COSMIC, COSMIC_phenotype_link) |>
  select(-COSMIC_PHENOTYPE_ID) |>
  mutate(cosmic_present = TRUE) |>
  distinct()

boostdm_hg19 = fread("processed_data/boostdm/boostdm_hg19.txt.gz") # Generated in the BoostDM_RefCDS.R script

# Save and filter mutation-site state
combine_boostdm_cancer = function(cosmic_select, genie_select, boostdm_muts) {

  cols = c("gene_name", "position", "pos", "alt", "aachange", "genie_present", "cosmic_present")

  cosmic_select = cosmic_select |> select(any_of(cols)) |> distinct()
  genie_select = genie_select |> select(any_of(cols)) |> distinct()
  obs_muts = full_join(cosmic_select, genie_select) |>
    distinct() |>
    mutate(cancer_present = TRUE)

  boostdm_cancer = left_join(boostdm_muts, obs_muts) |>
    mutate(across(c(cosmic_present, genie_present, cancer_present), ~ !is.na(.)),
           driver = (cosmic_present | genie_present) & boostDM_class)

  return(boostdm_cancer)
}


intersect_boostDM_cancer = function(boostdm_cohort, cosmic_cohort, genie_cohort) {

  genes = unique(boostdm_cohort[, gene_name])
  intersect_gene_list = list()
  for (gene in genes) {
    print(gene)
    cosmic_select = cosmic_cohort |> filter(gene_name == gene)
    genie_select =  genie_cohort |> filter(gene_name == gene)
    boostdm_select = boostdm_cohort |> filter(gene_name == gene)

    intersect_gene_list[[gene]] = combine_boostdm_cancer(cosmic_select =  cosmic_select, genie_select = genie_select,
                                                         boostdm_muts = boostdm_select)
  }
  rbindlist(intersect_gene_list)
}

euler_boostdm_cancer = function(boostdm_cancer, title) {
 mat = boostdm_cancer |> select("cosmic_present", "genie_present", "boostDM_class")*1 # multiply by 1 to convert to numeric
 colnames(mat) = c("Present in\nCosmic", "Mutation in GENIE", "Mutation is\nBoostDM Driver")
 set = eulerr::venn(mat)
 plot(set, main = title)
}

euler_boostdm_cancer_present = function(boostdm_cancer, title) {
  mat = boostdm_cancer |> select("cancer_present", "boostDM_class")*1 # multiply by 1 to convert to numeric
  colnames(mat) = c("Present in\nCosmic/GENIE", "Mutation is\nBoostDM Driver")
  set = eulerr::venn(mat)
  ggplotify::as.ggplot(plot(set, main = title))
}


# first get the observed mutations for the single samples. Merge on the change on the amino acid level
tissue = "pancancer"
boostdm_pancancer = boostdm_hg19[cohort == "CANCER", ]
bDM_cancer_intersect = intersect_boostDM_cancer(boostdm_cohort = boostdm_pancancer,
                                                cosmic_cohort = COSMIC_cancer_type,  genie_cohort = genie_cancer_type)
bdmFile = paste0("processed_data/boostdm/boostdm_genie_cosmic/", tissue, "_boostDM_intersect.txt")
fwrite(bDM_cancer_intersect, bdmFile)
R.utils::gzip(bdmFile)

  ##### Colon ######
tissue = "colon"
boostdm_colon = boostdm_hg19[cohort == "COADREAD", ]
cosmic_cohort = COSMIC_cancer_type |> filter(PRIMARY_SITE == "large_intestine")
genie_cohort =  genie_cancer_type |>   filter(ONCOTREE_CODE %in% c("COAD", "READ", "COADREAD"))

bDM_colon_intersect = intersect_boostDM_cancer(boostdm_colon, cosmic_cohort, genie_cohort)

bdmFile = paste0("processed_data/boostdm/boostdm_genie_cosmic/", tissue, "_boostDM_cancer.txt")
fwrite(bDM_colon_intersect, bdmFile)
R.utils::gzip(bdmFile, overwrite = TRUE)

#### Cosmic/Genie Intersection Colon ####
tissue = "colon_pancancer"
bDM_pancancer_colon_intersect = intersect_boostDM_cancer(boostdm_pancancer, cosmic_cohort, genie_cohort)
fwrite(bDM_pancancer_colon_intersect, paste0("processed_data/boostdm/boostdm_genie_cosmic/", tissue, "_boostDM_intersect.txt.gz"))


# make supplementary figure S3
euler_boostdm_cancer_present(bDM_colon_intersect |> filter(gene_name == "TP53"),
                             title = "TP53 Colon Driver")

euler_cohort_level = euler_boostdm_cancer(bDM_colon_intersect |> filter(gene_name == "TP53"),
                                          title = "TP53 Colon Driver")
ggsave("manuscript/Supplementary_Figures/Figure_S5/Figure_S5.png", euler_cohort_level, width = 8, height = 6)
ggsave("manuscript/Supplementary_Figures/Figure_S5/Figure_S5.svg", euler_cohort_level, width = 8, height = 6)


##### LUNG ######
tissue = "lung"
boostdm_lung = boostdm_hg19[cohort == "LUNG", ]
cosmic_cohort = COSMIC_cancer_type |> filter(PRIMARY_SITE == "lung")
genie_cohort =  genie_cancer_type |>   filter(CANCER_TYPE %in% c("Non-Small Cell Lung Cancer", "Small Cell Lung Cancer", "Lung Cancer"))
bDM_lung_intersect = intersect_boostDM_cancer(boostdm_lung, cosmic_cohort, genie_cohort)

bdmFile = paste0("processed_data/boostdm/boostdm_genie_cosmic/", tissue, "_boostDM_cancer.txt")
fwrite(bDM_lung_intersect, bdmFile)
R.utils::gzip(bdmFile)




# #### Blood ##### ## DEPRECATED NOW USING CH DATA
# tissue = "NON_SOLID"
# boostdm_ns = boostdm_hg19[cohort == "NON_SOLID", ]
# cosmic_cohort = COSMIC_cancer_type |> filter(PRIMARY_SITE == "haematopoietic_and_lymphoid_tissue")
# genie_cohort =  genie_cancer_type |>   filter(CANCER_TYPE %in%  c("Leukemia", "B-Lymphoblastic Leukemia/Lymphoma", "Myeloproliferative Neoplasms",
#                                                                   "Myelodysplastic Syndromes", "Mature T and NK Neoplasms",
#                                                                   "Myelodysplastic/Myeloproliferative Neoplasms"))
# bDM_ns_intersect = intersect_boostDM_cancer(boostdm_ns, cosmic_cohort, genie_cohort)
#
# bdmFile = paste0("processed_data/boostdm/boostdm_genie_cosmic/", tissue, "_boostDM_cancer.txt")
# fwrite(bDM_ns_intersect, bdmFile)
# R.utils::gzip(bdmFile)

fwrite(bDM_lung_intersect, paste0("processed_data/boostdm/boostdm_genie_cosmic/", tissue, "_boostDM_cancer.txt.gz"))

#### Cosmic/Genie Intersection Colon ####
tissue = "lung_pancancer"
bDM_pancancer_lung_intersect = intersect_boostDM_cancer(boostdm_pancancer, cosmic_cohort, genie_cohort)
fwrite(bDM_pancancer_lung_intersect, paste0("processed_data/boostdm/boostdm_genie_cosmic/", tissue, "_boostDM_intersect.txt.gz"))


#### Blood #####
tissue = "NON_SOLID"
boostdm_ns = boostdm_hg19[cohort == "NON_SOLID", ]
cosmic_cohort = COSMIC_cancer_type |> filter(PRIMARY_SITE == "haematopoietic_and_lymphoid_tissue")
genie_cohort =  genie_cancer_type |>   filter(CANCER_TYPE %in%  c("Leukemia", "B-Lymphoblastic Leukemia/Lymphoma", "Myeloproliferative Neoplasms",
                                                                  "Myelodysplastic Syndromes", "Mature T and NK Neoplasms",
                                                                  "Myelodysplastic/Myeloproliferative Neoplasms"))
bDM_ns_intersect = intersect_boostDM_cancer(boostdm_ns, cosmic_cohort, genie_cohort)
fwrite(bDM_ns_intersect, paste0("processed_data/boostdm/boostdm_genie_cosmic/", tissue, "_boostDM_cancer.txt.gz"))

#### Cosmic/Genie Intersection Colon ####
tissue = "blood_pancancer"
bDM_pancancer_blood_intersect = intersect_boostDM_cancer(boostdm_pancancer, cosmic_cohort, genie_cohort)
fwrite(bDM_pancancer_blood_intersect, paste0("processed_data/boostdm/boostdm_genie_cosmic/", tissue, "_boostDM_intersect.txt.gz"))


# Using Clonal hematopoiesiss - Blood CH data
boostdm_ch_hg19 = fread("processed_data/boostdm_ch/boostdm_ch_hg19.txt.gz") # BoostDM_ch_RefCDS.R script

#### TP53 in Clonal hematopoiesis
tissue = "CH"
cosmic_cohort = COSMIC_cancer_type |> filter(PRIMARY_SITE == "haematopoietic_and_lymphoid_tissue")
genie_cohort =  genie_cancer_type |>   filter(CANCER_TYPE %in%  c("Leukemia", "B-Lymphoblastic Leukemia/Lymphoma", "Myeloproliferative Neoplasms",
                                                                                      "Myelodysplastic Syndromes", "Mature T and NK Neoplasms", "Myelodysplastic/Myeloproliferative Neoplasms"))
bDM_ch_intersect = intersect_boostDM_cancer(boostdm_ch_hg19, cosmic_cohort, genie_cohort)

bdmFile = paste0("processed_data/boostdm/boostdm_genie_cosmic/", tissue, "_boostDM_cancer.txt")
fwrite(bDM_ch_intersect, bdmFile)
R.utils::gzip(bdmFile)

#### Generate new table
bDM_lung_intersect = bDM_lung_intersect |>
  select(-cosmic_present, -genie_present) |>
  mutate(driver = boostDM_class)
bDM_colon_intersect = bDM_colon_intersect |>
  select(-cosmic_present, -genie_present) |>
  mutate(driver = boostDM_class)
bDM_ch_intersect = bDM_ch_intersect |>
  select(-cosmic_present, -genie_present) |>
  mutate(driver = boostDM_class)

list = list(
  "Lung cancer" = bDM_lung_intersect,
  "Colon cancer" = bDM_colon_intersect,
  "Clonal hematopoiesis" = bDM_ch_intersect)
writexl::write_xlsx(list, "manuscript/Supplementary_Tables/Supplementary_Table_6.xlsx")