library(figpatch)
library(gridExtra)
library(cowplot)
source("code/0_functions/analysis_variables.R")

# Figure 1:
# load figure 1 graphic:
print("Composing Figure1")
mg = 7 # set margins

# Load Figure 1A if present, otherwise use an empty placeholder
if (file.exists("manuscript/figure_panels/figure_1/Figure 1A.png")) {
  figure_1A <- fig("manuscript/figure_panels/figure_1/Figure 1A.png") |> prep_plot(label = "A", all_margin = 0)
} else {
  print("WARNING: Figure 1A.png not found — using empty placeholder")
  figure_1A <- ggplot() + theme_void()
}

figure_1B = readRDS("manuscript/figure_panels/figure_1/figure_1B.rds") |> prep_plot(label = "B", all_margin = mg)
figure_1C = readRDS("manuscript/figure_panels/figure_1/figure_1C.rds") |> prep_plot(label = "C", all_margin = mg)
figure_1D = readRDS("manuscript/figure_panels/figure_1/figure_1D.rds") |> prep_plot(label = "D", all_margin = mg)

figure_1_middle = figure_1B # plot_spacer() + plot_layout(widths = c(2.8, 0.7))
figure_1_bottom = figure_1C + figure_1D + plot_layout(widths = c(2.5, 1))
figure_1  = figure_1A / figure_1_middle / figure_1_bottom + plot_layout(heights =  c(1,1, 1))

ggsave("manuscript/Figure_1/figure_1.png", figure_1, width = 16, height = 13)
ggsave("manuscript/Figure_1/figure_1.pdf", figure_1, width = 16, height = 13)
ggsave("manuscript/Figure_1/figure_1.svg", figure_1, width = 16, height = 13)

print("Figure1 saved in manuscript/Figure_1")

#### Figure 2

print("Composing Figure 2")
mg = 10 # set margins
figure_2A = readRDS("manuscript/figure_panels/figure_2/figure_2A.rds") |> prep_plot(label = "A", all_margin = mg)

# Load Figure_2B_small.png if present, otherwise use an empty placeholder
if (file.exists("manuscript/figure_panels/figure_2/Figure_2B_small.png")) {
  figure_2B = fig("manuscript/figure_panels/figure_2/Figure_2B_small.png") |> prep_plot(label = "B", all_margin = 0)
} else {
  print("WARNING: Figure_2B_small.png not found — using empty placeholder")
  figure_2B <- ggplot() + theme_void()
}

figure_2C = readRDS("manuscript/figure_panels/figure_2/figure_2C.rds") |> prep_plot(label = "C", all_margin = mg)
figure_2D = readRDS("manuscript/figure_panels/figure_2/figure_2D.rds") |> prep_plot(label = "D", all_margin = mg)

figure_2_middle =  figure_2B + figure_2C
figure_2 = figure_2A / figure_2_middle / figure_2D + plot_layout(heights = c(0.8, 1,1))
ggsave("manuscript/Figure_2/figure_2.png", figure_2, width = 12, height = 12)
ggsave("manuscript/Figure_2/figure_2.pdf", figure_2, width = 12, height = 12)
ggsave("manuscript/Figure_2/figure_2.svg", figure_2, width = 12, height = 12)

print("Figure2 saved in manuscript/Figure_2/")

#### Figure 3
print("Composing Figure 3")
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
cell_color = c("white",    "white",   "white",
                "black",    "black",   "black",
                "black",    "black",   "black")

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

# Load Figure 1A if present, otherwise use an empty placeholder
if (file.exists("manuscript/figure_panels/figure_3/figure_3B_only_boostDM.png")) {
figure_3B <- fig("manuscript/figure_panels/figure_3/figure_3B_only_boostDM.png") |> prep_plot(label = "B", all_margin = 4)
} else {
  print("WARNING: Figure 3B_only_boostDM.png not found — using empty placeholder")
  figure_3B <- ggplot() + theme_void()
}

figure_3B_column = (figure_3B / wrap_elements(tp53_table)) + plot_layout(heights = c(2, 1))
figure_3_top = figure_3A + figure_3B_column + plot_layout(widths = c(2.5, 1))

tissue_plots_raw = readRDS("manuscript/figure_panels/figure_3/figure_3CDE.rds")
tissue_plots = mapply(prep_plot, tissue_plots_raw, label = c("C", "D", "E"), all_margin = mg)

figure_3_bottom = wrap_plots(tissue_plots, nrow = 1, widths = c(4, 3, 1.2))
figure_3 = figure_3_top / figure_3_bottom + plot_layout(heights = c(1.6, 1))

ggsave("manuscript/Figure_3/figure_3.png", figure_3, width = 16, height = 11)
ggsave("manuscript/Figure_3/figure_3.pdf", figure_3, width = 16, height = 11)
ggsave("manuscript/Figure_3/figure_3.svg", figure_3, width = 16, height = 11)

print("Figure3 saved in manuscript/Figure_3/")

##### Figure 4

print("Composing Figure 4")
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
ggsave("manuscript/Figure_4/figure_4.pdf", figure_4, width = 18, height = 11)

print("Figure4 saved in manuscript/Figure_4/")

##### Figure 5: APC/KRAS progression

print("Composing Figure 5")
mg <- 5

figure_5A <- readRDS("manuscript/figure_panels/figure_5/figure_5A.rds") |> prep_plot("A", all_margin = mg)
figure_5B <- readRDS("manuscript/figure_panels/figure_5/figure_5B.rds") |> prep_plot("B", all_margin = mg)
figure_5C <- readRDS("manuscript/figure_panels/figure_5/figure_5C.rds") |> prep_plot("C", all_margin = mg)
figure_5D <- readRDS("manuscript/figure_panels/figure_5/figure_5D.rds") |> prep_plot("D", all_margin = mg)

figure_5 <- cowplot::plot_grid(
  figure_5A, figure_5B,
  figure_5C, figure_5D,
  labels = NULL, ncol = 2, align = "hv"
)

dir.create("manuscript/Figure_5", showWarnings = FALSE, recursive = TRUE)
ggsave("manuscript/Figure_5/figure_5.png", figure_5, width = 13, height = 10.5, dpi = 300, bg = "white")
ggsave("manuscript/Figure_5/figure_5.pdf", figure_5, width = 13, height = 10.5)
ggsave("manuscript/Figure_5/figure_5.svg", figure_5, width = 13, height = 10.5)

print("Figure4 saved in manuscript/Figure_5/")

##### Figure 6

print("Composing Figure 6")
figure_6A = readRDS("manuscript/figure_panels/figure_6/figure_6A.rds") |> prep_plot("A", all_margin = 3)
# new row under panel A: Figure S8 panels A/B + binned Figure S9B
figure_6B = readRDS("manuscript/figure_panels/figure_6/figure_6B.rds")
figure_6C = readRDS("manuscript/figure_panels/figure_6/figure_6C.rds")
figure_6D = readRDS("manuscript/figure_panels/figure_6/figure_6D.rds")

# save final completed plot
figure_6 = figure_6A / figure_6B / figure_6C / figure_6D
ggsave("manuscript/Figure_6/Figure_6.pdf", figure_6, width = 16, height = 19)
ggsave("manuscript/Figure_6/Figure_6.png", figure_6, width = 16, height = 19)
ggsave("manuscript/Figure_6/Figure_6.svg", figure_6, width = 16, height = 19)
print("Figure6 saved in manuscript/Figure_6/")

# Figure S5

print("Composing Figure S5")
figure_S5A = readRDS("manuscript/figure_panels/figure_s5/Figure_S5A.rds") |> prep_plot("A")
figure_S5B = readRDS("manuscript/figure_panels/figure_s5/Figure_S5B.rds") |> prep_plot("B")
figure_S5C = readRDS("manuscript/figure_panels/figure_s5/Figure_S5C.rds") |> prep_plot("C")
figure_S5 = (figure_S5A | figure_S5B) / figure_S5C
ggsave("manuscript/Supplementary_Figures/Figure_S5/Figure_S5.png", figure_S5, width = 14, height = 12)
ggsave("manuscript/Supplementary_Figures/Figure_S5/Figure_S5.svg", figure_S5, width = 14, height = 12)
ggsave("manuscript/Supplementary_Figures/Figure_S5/Figure_S5.pdf", figure_S5, width = 14, height = 12)

print("Figure S5 saved in manuscript/Supplementary_Figures/Figure_S5/")

# Figure S6

mg = 8

print("Composing Figure S6")
figure_S6A = readRDS("manuscript/figure_panels/figure_s6/Figure_S6A.rds") |> prep_plot("A", all_margin = mg)
figure_S6_top = plot_spacer() + figure_S6A + plot_spacer() + plot_layout(widths = c(1, 3, 1)) 

figure_S6B = readRDS("manuscript/figure_panels/figure_s6/Figure_S6B.rds") |> prep_plot("B", all_margin = mg)
figure_S6 = figure_S6_top / figure_S6B + plot_layout(heights = c(1.5, 2))

ggsave("manuscript/Supplementary_Figures/Figure_S6/Figure_S6.png", figure_S6, width = 11, height = 10)
ggsave("manuscript/Supplementary_Figures/Figure_S6/Figure_S6.svg", figure_S6, width = 11, height = 10)
ggsave("manuscript/Supplementary_Figures/Figure_S6/Figure_S6.pdf", figure_S6, width = 11, height = 10)

print("Figure S6 saved in manuscript/Supplementary_Figures/Figure_S6/")