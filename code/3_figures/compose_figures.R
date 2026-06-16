library(figpatch)
library(gridExtra)
source("code/0_functions/analysis_variables.R")

# Figure 1:
# load figure 1 graphic:
mg = 7 # set margins
figure_1A <- fig("manuscript/figure_panels/figure_1/Figure 1A.png") |> prep_plot(label = "A", all_margin = 0)
figure_1B = readRDS("manuscript/figure_panels/figure_1/figure_1B.rds") |> prep_plot(label = "B", all_margin = mg)
figure_1C = readRDS("manuscript/figure_panels/figure_1/figure_1C.rds") |> prep_plot(label = "C", all_margin = mg)
figure_1D = readRDS("manuscript/figure_panels/figure_1/figure_1D.rds") |> prep_plot(label = "D", all_margin = mg)

figure_1_middle = figure_1B # plot_spacer() + plot_layout(widths = c(2.8, 0.7))
figure_1_bottom = figure_1C + figure_1D + plot_layout(widths = c(2.5, 1))
figure_1  = figure_1A / figure_1_middle / figure_1_bottom + plot_layout(heights =  c(1,1, 1))

ggsave("manuscript/Figure_1/figure_1.png", figure_1, width = 16, height = 13)
ggsave("manuscript/Figure_1/figure_1.pdf", figure_1, width = 16, height = 13)
ggsave("manuscript/Figure_1/figure_1.svg", figure_1, width = 16, height = 13)

#### Figure 2
mg = 10 # set margins
figure_2A = readRDS("manuscript/figure_panels/figure_2/figure_2A.rds") |> prep_plot(label = "A", all_margin = mg)
figure_2B = fig("manuscript/figure_panels/figure_2/Figure_2B_small.png") |> prep_plot(label = "B", all_margin = 0)
figure_2C = readRDS("manuscript/figure_panels/figure_2/figure_2C.rds") |> prep_plot(label = "C", all_margin = mg)
figure_2D = readRDS("manuscript/figure_panels/figure_2/figure_2D.rds") |> prep_plot(label = "D", all_margin = mg)

figure_2_middle =  figure_2B + figure_2C
figure_2 = figure_2A / figure_2_middle / figure_2D + plot_layout(heights = c(0.8, 1,1))
ggsave("manuscript/Figure_2/figure_2.png", figure_2, width = 12, height = 12)
ggsave("manuscript/Figure_2/figure_2.pdf", figure_2, width = 12, height = 12)

#### Figure 3
mg = 8

# Load TP53 site counts for table under figure 3B
tp53_site_counts = readRDS("manuscript/figure_panels/figure_3/figure_3B_table.rds")
# Ensure consistent tissue order: colon, blood, lung
tp53_site_counts = tp53_site_counts[match(c("colon", "blood", "lung"), tissue)]

# Colors: colon = #3ca951, blood = #ff725c, lung = #4269d0
tp53_table_data = data.frame(
  Tissue = c("Colon", "Blood", "Lung"),
  Driver = tp53_site_counts$driver_sites,
  NonDriver = tp53_site_counts$nondriver_sites,
  stringsAsFactors = FALSE
)
colnames(tp53_table_data) = c("Tissue", "TP53 BoostDM\ndriver SNV sites", "TP53 BoostDM\nnon-driver SNV sites")

# Per-cell styling: tissue column gets colored background + white bold text
cell_fill = c("#3ca951", "#ff725c",   "#4269d0",
              "white", "white",   "white",
              "white", "white",   "white")
cell_font  = c("bold",    "bold",   "bold",
               "plain",    "plain",   "plain",
               "plain",    "plain",   "plain")
cell_color = c("white",   "white",   "white",
               "black",   "black",   "black",
               "black",   "black",   "black")

table_theme = ttheme_default(
  base_size = 9,
  padding = unit(c(8, 4), "mm"),
  core = list(
    bg_params = list(fill = cell_fill, col = "grey80"),
    fg_params = list(fontface = cell_font, col = cell_color)
  ),
  colhead = list(
    bg_params = list(fill = "grey90"),
    fg_params = list(fontface = "bold"),
    padding = unit(c(10, 4), "mm")
  )
)

tp53_table = tableGrob(tp53_table_data, rows = NULL, theme = table_theme)

figure_3A = readRDS("manuscript/figure_panels/figure_3/figure_3A.rds") |> prep_plot(label = "A", all_margin = mg)
figure_3B = fig("manuscript/figure_panels/figure_3/figure_3B_only_boostDM.png") |> prep_plot(label = "B", all_margin = mg)
figure_3B_column = (figure_3B / wrap_elements(tp53_table)) + plot_layout(heights = c(2, 1))
figure_3_top = figure_3A + figure_3B_column + plot_layout(widths = c(2.5, 1))

tissue_plots_raw = readRDS("manuscript/figure_panels/figure_3/figure_3CDE.rds")
tissue_plots = mapply(prep_plot, tissue_plots_raw, label = c("C", "D", "E"), all_margin = mg)

figure_3_bottom = wrap_plots(tissue_plots, nrow = 1, widths = c(4, 3, 1.2))
figure_3 = figure_3_top / figure_3_bottom + plot_layout(heights = c(1.6, 1))

ggsave("manuscript/Figure_3/figure_3.bdm.v2.6.png", figure_3, width = 16, height = 11)
ggsave("manuscript/Figure_3/figure_3.bdm.v2.6.pdf", figure_3, width = 16, height = 11)


##### Figure 4
mg = 5

# perform these operations later in the compose-figures part:
list_figure_4AB = readRDS("manuscript/figure_panels/figure_4/figures_AB.rds")
figure_4A = list_figure_4AB[[1]] |> prep_plot(label = 'A', all_margin = mg)
figure_4B = list_figure_4AB[[2]]  |> prep_plot(label = 'B', all_margin = mg)
figure_4_top = figure_4A + figure_4B

figs = readRDS("manuscript/figure_panels/figure_4/figures_C_F_H.rds")
names(figs) = c("C", "F", "H")
annotated_figs = lapply(names(figs), \(x) prep_plot(figs[[x]], label = x, all_margin = mg))
figure_4_middle = wrap_plots(c(annotated_figs[1],  annotated_figs[2:3]), nrow = 1)

figure_4D = readRDS("manuscript/figure_panels/figure_4/figures_adenoma.rds") |>
  prep_plot(label = 'D', all_margin = mg)

figures_4_CRC = readRDS("manuscript/figure_panels/figure_4/figure_4_CRC.rds")
figures_4_CRC = figures_4_CRC[-c(2,4)]
annotated_figs = mapply(prep_plot, figures_4_CRC, c("E", "G", "I"), all_margin = 2)
figure_4_bottom = wrap_plots(c(figure_4D, annotated_figs), nrow = 1)

figure_4 = figure_4_top / figure_4_middle / figure_4_bottom# + plot_layout(heights = c(1.5,1, 1))

ggsave("manuscript/Figure_4/figure_4.png", figure_4, width = 18, height = 11)
ggsave("manuscript/Figure_4/figure_4.pdf", figure_4, width = 18, height = 11)

##### Figure 4 A-E: APC/KRAS progression (BoostDM-only analysis)
# This figure is generated by running:
#   Rscript code/3_figures/41_figure_apc_progression_AE_boostdm.R
# The output is saved in manuscript/figure_panels/figure_4/figure_apc_progression_AE_boostdm.{png,pdf,svg}
# Panels:
#   A: Adenoma incidence with APC mutation
#   B: Single APC driver cells vs CRC cases
#   C: Double APC driver cells vs CRC cases
#   D: Double APC + KRAS driver cells vs CRC cases
#   E: Cumulative CRC vs carrier cells (existing figure)

# If you need to compose it here manually, load the saved panels:
# panels_ae = readRDS("manuscript/figure_panels/figure_4/figure_apc_progression_panels.rds")
# figure_4ae = (panels_ae$A + panels_ae$B) / (panels_ae$C + panels_ae$D) / panels_ae$E
# ggsave("manuscript/Figure_4/figure_4_AE.png", figure_4ae, width = 14, height = 16, dpi = 300)

##### Figure 5
figure_5A = readRDS("manuscript/figure_panels/figure_5/figure_5A.rds") |> prep_plot("A", all_margin = 3)
figure_5B = readRDS("manuscript/figure_panels/figure_5/figure_5B.rds")
figure_5C = readRDS("manuscript/figure_panels/figure_5/figure_5C.rds")
figure_5D = readRDS("manuscript/figure_panels/figure_5/figure_5D.rds")

# save final completed plot
figure_5 = figure_5A / figure_5B / figure_5C / figure_5D
ggsave("manuscript/Figure_5/Figure_5.pdf", width = 13, height = 15)
ggsave("manuscript/Figure_5/Figure_5.png", width = 13, height = 15)

# Figure S6
figure_S6A = readRDS("manuscript/Supplementary_Figures/Figure_S6/figure_S6A.rds") |> prep_plot("A")
figure_S6B = readRDS("manuscript/Supplementary_Figures/Figure_S6/Figure_S6B.rds") |> prep_plot("B")
figure_S6C = readRDS("manuscript/Supplementary_Figures/Figure_S6/Figure_S6C.rds") |> prep_plot("C")
figure_S6 = (figure_S6A | figure_S6B) / figure_S6C
<<<<<<< Updated upstream
ggsave("manuscript/Supplementary_Figures/Figure_S6/Figure_S6.svg", figure_S6, width = 14, height = 12)
=======
ggsave("manuscript/Supplementary_Figures/Figure_S6/Figure_S6.png", figure_S6, width = 14, height = 12)
>>>>>>> Stashed changes
