# Supplement figures

run_supplement_figures <- function(project_root) {
  suppressPackageStartupMessages({
    library(grid)
    library(ggplot2)
    library(patchwork)
    library(dplyr)
    library(cowplot)
  })
  root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  sys.source(file.path(root, "R/01_figure_style.R"), envir = environment())
  output_folder <- file.path(root, "Manuscript/Supplement/Supplement_figures")
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

  png(file.path(output_folder, "Figure_S1_desiccation_design.png"), width = 2700, height = 810,
    res = 300, bg = "white")
  grid.newpage()
  column_widths <- c(0.19, 0.19, 0.23, 0.19, 0.2)
  column_edges <- c(0, cumsum(column_widths))
  column_headers <- c("Experiment", "Preconditioning\ncomparison", "Desiccation\ncomparison",
    "Desiccation\nconditions", "Rehydration")
  experiment_rows <- list(c("Desiccation\nexperiment 1 (D1)", "None", "Non-desiccated control\nvs desiccated treatment",
    "8 h\n<10% RH", "M9 buffer\n17 h"), c("Desiccation\nexperiment 2 (D2)", "0 h vs 96 h\n97–98% RH",
    "Both groups\ndesiccated", "48 h\n<10% RH", "dH2O, 24 h\nshaking at 50 rpm"))
  column_colours <- c("#eeeeee", "#fff2cc", "#f4cccc", "#cfe8f3", "#d9d2e9")
  for (column_index in 1:5) {
    for (row_index in 0:2) {
      row_centre <- if (row_index == 0)
        0.83 else if (row_index == 1)
        0.515 else 0.185
      row_height <- if (row_index == 0)
        0.28 else 0.33
      grid.rect(x = column_edges[column_index] + column_widths[column_index]/2, y = row_centre,
        width = column_widths[column_index], height = row_height, gp = gpar(fill = if (row_index ==
          0)
          column_colours[column_index] else adjustcolor(column_colours[column_index], alpha.f = 0.25), col = "#888", lwd = 0.7))
      grid.text(if (row_index == 0)
        column_headers[column_index] else experiment_rows[[row_index]][column_index], x = column_edges[column_index] + column_widths[column_index]/2,
        y = row_centre, gp = gpar(fontsize = if (row_index == 0)
          11.5 else 10.5, fontface = if (row_index == 0)
          "bold" else "plain", fontfamily = "sans"))
    }
  }
  dev.off()

  statistics_folder <- file.path(root, "Results/statistics/salinity")
  plate_recovery <- read.csv(file.path(statistics_folder, "m9_recovery_analysis_ready.csv"))
  pooled_recovery <- read.csv(file.path(statistics_folder, "m9_recovery_by_salinity.csv"))
  tested_salinities <- read.csv(file.path(statistics_folder, "m9_recovery_tested_salinity_grid.csv"))
  species_labels <- c(ALT.22.08 = "bolditalic(Panagrolaimus)~bold(\"sp.\")~plain(\"(WWL002)\")",
    PAP.22.17 = "bolditalic(Acrobeloides~tricornis)~plain(\"(WWL709)\")", `A. max` = "bolditalic(Acrobeloides~maximus)~plain(\"(A. max)\")")

  make_recovery_panel <- function(id, strain) {
    plate_points <- filter(plate_recovery, Experiment_id == id)
    pooled_points <- filter(pooled_recovery, Experiment_id == id)
    empty_salinity_levels <- filter(tested_salinities, Experiment_id == id, !has_cryptobiotic_candidates)
    pooled_points <- mutate(pooled_points, label = ifelse(recovered_animals > 0, paste0(recovered_animals,
      "/", cryptobiotic_candidates), ""), label_y = ifelse(recovered_fraction_candidate >
      0.95, 100 * recovered_fraction_candidate - 9, 100 * recovered_fraction_candidate + 8))
    ggplot() +
      geom_point(data = plate_points, aes(salinity_PSU, 100 * recovered_fraction_candidate),
      colour = "#999999", alpha = 0.5, size = 1.7, position = position_jitter(width = 0.15,
        height = 0, seed = 1211)) +
      geom_point(data = pooled_points, aes(salinity_PSU, 100 *
      recovered_fraction_candidate), shape = 21, fill = figure_style$blue, colour = "#333333",
      size = 2.8, stroke = 0.65) +
      geom_point(data = empty_salinity_levels, aes(salinity_PSU,
      -1.5), shape = 4, colour = "#666666", size = 2.8, stroke = 0.9) +
      geom_text(data = pooled_points,
      aes(salinity_PSU, label_y, label = label), size = 3.4, colour = "#333333") +
      facet_wrap(~factor(species_labels[[strain]]),
      labeller = label_parsed) +
      scale_x_continuous(breaks = seq(0, 40, 10)) +
      scale_y_continuous(breaks = seq(0,
      100, 25)) +
      coord_cartesian(xlim = c(0, 40), ylim = c(-4, 104)) +
      labs(x = NULL, y = NULL,
      title = if (grepl("ALT", id))
        paste0(unique(plate_points$temperature_C), " °C") else NULL) +
      figure_theme() +
      theme(plot.margin = margin(4, 5, 4, 5), plot.title = element_text(hjust = 0.5,
      size = 16, face = "bold"), plot.title.position = "panel")
  }

  recovery_panels <- ((make_recovery_panel("ALT.22.08_22C", "ALT.22.08") + make_recovery_panel("ALT.22.08_16C",
    "ALT.22.08"))/(make_recovery_panel("A.max_22C", "A. max") + make_recovery_panel("PAP.22.17_16C",
    "PAP.22.17"))) + plot_layout(heights = c(1, 1))
  recovery_plot <- ggdraw() +
      draw_plot(recovery_panels, x = 0.065, y = 0.065, width = 0.935,
    height = 0.935) +
      draw_label("Salinity (PSU)", x = 0.55, y = 0.06, size = 14) +
      draw_label("Recovery (%)",
    x = 0.06, y = 0.54, angle = 90, size = 14)
  recovery_temp <- file.path(root, ".tmp", "Figure_S2_recovery_render.png")
  ggsave(recovery_temp, recovery_plot, width = 9, height = 7.5, dpi = 300, bg = "white")
  if (!file.copy(recovery_temp, file.path(output_folder, "Figure_S2_recovery.png"), overwrite = TRUE))
    warning("Rendered figure must be copied to final destination")


  dir.create(file.path(root, "Results/cache/supplement"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(points = plate_recovery, pooled = pooled_recovery, tested = tested_salinities),
    file.path(root, "Results/cache/supplement/recovery_plot_inputs.rds"))
  invisible(output_folder)
}