# Dual-axis plotting functions for CRC vs cell count figures

library(ggplot2)
library(scales)

# Colors matching Python script
COLOR_OUTCOME = "#4a4a4a"  # Dark gray for CRC cases
COLOR_CELLS = "#2e8b57"    # Sea green for cell counts

#' Style axes for dual-axis plots
#'
#' @param ax ggplot object
#' @param y_color Color for y-axis (right side)
#' @return Modified ggplot object
style_axes_dual = function(ax, y_color = COLOR_OUTCOME) {
  ax = ax +
    theme(
      axis.line.y = element_line(color = y_color),
      axis.text.y = element_text(color = y_color),
      axis.title.y = element_text(color = y_color),
      axis.line.y.right = element_line(color = COLOR_CELLS),
      axis.text.y.right = element_text(color = COLOR_CELLS),
      axis.title.y.right = element_text(color = COLOR_CELLS),
      axis.line.x = element_line(color = "#404040"),
      axis.text.x = element_text(color = "#404040"),
      axis.title.x = element_text(color = "#404040"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.margin = margin(10, 35, 10, 10, "pt")  # Extra right margin for second y-axis
    )

  return(ax)
}

#' Create dual-axis plot: CRC cases vs carrier cells
#'
#' @param crc_data Data frame with age and crc_cases_per_million
#' @param cells_data Data frame with age, cell_count, cell_count_lo, cell_count_hi
#' @param left_ylabel Label for left y-axis (CRC cases)
#' @param right_ylabel Label for right y-axis (cell count)
#' @param x_max Maximum x-axis value (default 85)
#' @param y_max_override Override for maximum y-axis value (optional)
#' @param scientific_right Use scientific notation for right axis (default FALSE)
#' @return ggplot object with dual axes
plot_crc_cells_dual_axis = function(
    crc_data,
    cells_data,
    left_ylabel = "CRC cases per 1,000,000 people",
    right_ylabel = "Carrier cells",
    x_max = 85,
    y_max_override = NULL,
    scientific_right = FALSE) {

  # Ensure data is sorted by age
  crc_data = crc_data %>% arrange(age)
  cells_data = cells_data %>% arrange(age)

  # Calculate y-axis limits
  crc_max = max(crc_data$crc_cases_per_million, na.rm = TRUE)
  cells_max = max(cells_data$cell_count_hi, na.rm = TRUE)

  if (!is.null(y_max_override)) {
    crc_max = y_max_override
  } else {
    crc_max = crc_max * 1.15
  }
  cells_max = cells_max * 1.08

  # Build plot
  p = ggplot() +
    # CRC cases (left axis)
    geom_line(
      data = crc_data,
      aes(x = age, y = crc_cases_per_million),
      color = COLOR_OUTCOME,
      linewidth = 1.1
    ) +
    geom_point(
      data = crc_data,
      aes(x = age, y = crc_cases_per_million),
      color = COLOR_OUTCOME,
      size = 2
    ) +
    # Cell counts with error bars (right axis)
    geom_errorbar(
      data = cells_data,
      aes(x = age, ymin = cell_count_lo, ymax = cell_count_hi),
      color = COLOR_CELLS,
      width = 0,
      linewidth = 0.5
    ) +
    geom_point(
      data = cells_data,
      aes(x = age, y = cell_count),
      color = COLOR_CELLS,
      size = 2.5,
      alpha = 0.9
    ) +
    # Axis settings
    scale_x_continuous(
      limits = c(0, x_max),
      breaks = c(0, 20, 40, 60, 80),
      name = "Age (years)"
    ) +
    scale_y_continuous(
      limits = c(0, crc_max),
      breaks = pretty_breaks(n = 5),
      labels = if (scientific_right) scientific_format() else comma,
      name = left_ylabel,
      sec.axis = sec_axis(
        transform = ~ .,
        breaks = pretty_breaks(n = 5),
        labels = if (scientific_right) scientific_format() else comma,
        name = right_ylabel
      )
    ) +
    # Color styling
    scale_color_manual(
      values = c(COLOR_OUTCOME, COLOR_CELLS),
      name = "",
      labels = c(left_ylabel, right_ylabel)
    ) +
    # Theme
    theme_cowplot() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      axis.line.y = element_line(color = COLOR_OUTCOME),
      axis.text.y = element_text(color = COLOR_OUTCOME),
      axis.title.y = element_text(color = COLOR_OUTCOME, face = "bold"),
      axis.line.y.right = element_line(color = COLOR_CELLS),
      axis.text.y.right = element_text(color = COLOR_CELLS),
      axis.title.y.right = element_text(color = COLOR_CELLS, face = "bold"),
      axis.title.x = element_text(color = "#404040", face = "bold"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    )

  return(p)
}

#' Simplified dual-axis plot for adenoma incidence (single axis)
#'
#' @param adenoma_data Data frame with age and fraction_adenoma_apc
#' @param y_max Maximum y-axis value (default 0.75)
#' @return ggplot object
plot_adenoma_incidence = function(adenoma_data, y_max = 0.75) {
  p = ggplot(adenoma_data, aes(x = age, y = fraction_adenoma_apc)) +
    geom_point(color = COLOR_OUTCOME, size = 2) +
    geom_line(color = COLOR_OUTCOME, linewidth = 1.1) +
    scale_x_continuous(
      limits = c(20, 85),
      breaks = c(20, 40, 60, 80),
      name = "Age (years)"
    ) +
    scale_y_continuous(
      limits = c(0, y_max),
      breaks = pretty_breaks(n = 5),
      labels = percent_format(accuracy = 1),
      name = "Lifetime risk for Adenoma\nwith APC mutation"
    ) +
    theme_cowplot() +
    theme(
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    )

  return(p)
}
