# Load packages. ----
library(export)
library(gridExtra)
library(mosaic)
library(readxl)
library(tidyverse)

# Read and clean results. ----
result <- readxl::read_xlsx(
  path = "results_constant_rates_n500_reformatted_for_r.xlsx",
  sheet = "main"
)
result <- result |>
  dplyr::mutate(
    scenario_informal = as.factor(scenario_informal),
    scenario = factor(
      scenario,
      levels = c("Random", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8",
                 "Q9", "Q10", "Q11", "Q12")
    ),
    xy = as.factor(xy),
    is_max = as.factor(is_max)
  )

# Process results for Figure 2. ----
result <- result |>
  dplyr::select(
    scenario, log_loss_bb, log_loss_nb, log_loss_ind_rj, log_loss_dep_rj
  ) |>
  tidyr::pivot_longer(
    cols = 2:5,
    names_to = "method",
    values_to = "log_loss"
  ) |>
  dplyr::mutate(method = as.factor(method)) |>
  dplyr::mutate(
    method = forcats::fct_relevel(
      method,
      c("log_loss_bb", "log_loss_nb", "log_loss_ind_rj", "log_loss_dep_rj")
    )
  ) |>
  dplyr::mutate(
    method = forcats::fct_recode(
      method,
      "Beta-\nBinomial" = "log_loss_bb",
      "Naive\nBayes" = "log_loss_nb",
      "RJ-\nIndependent" = "log_loss_ind_rj",
      "RJ-\nDependent" = "log_loss_dep_rj"
    )
  ) |>
  dplyr::arrange(scenario, method, log_loss)

mosaic::favstats(log_loss ~ method, data = result)
#>             method     min      Q1  median      Q3      max      mean        sd
#> 1  Beta-\nBinomial 0.01511 0.31334 0.57270 0.79851  3.38139 0.6265445 0.4323811
#> 2     Naive\nBayes 0.00001 0.23081 0.48583 0.74533 11.51293 0.5778947 0.6473719
#> 3 RJ-\nIndependent 0.00001 0.01309 0.03874 0.13239  6.90776 0.2381866 0.5699230
#> 4   RJ-\nDependent 0.00001 0.00904 0.03459 0.12897 11.51293 0.2420802 0.6304005
#>       n missing
#> 1 13000       0
#> 2 13000       0
#> 3 13000       0
#> 4 13000       0

# Create the individual plots. ----
result %>%
  dplyr::filter(scenario %in% "Random") %>%
  ggplot2::ggplot(mapping = aes(x = method, y = log_loss, color = method)) +
    ylim(0, max(result$log_loss)) +
    geom_jitter(alpha = 0.25, size = 0.75) +
    geom_violin(fill = "black", color = "black") +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 21,
      size = 1.25,
      color = "black",
      fill = "white"
    ) +
    scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
    theme_minimal(base_size = 6) +
    theme(legend.position = "none") +
    labs(x = NULL, y = "Log loss")
export::graph2svg(
  file = "figure2_plot_random",
  width = 1.8,
  height = 1.75,
  font = "Arial",
  bg = "transparent",
)
result %>%
  dplyr::filter(scenario %in% "Q1") %>%
  ggplot2::ggplot(mapping = aes(x = method, y = log_loss, color = method)) +
    ylim(0, max(result$log_loss)) +
    geom_jitter(alpha = 0.25, size = 0.75) +
    geom_violin(fill = "black", color = "black") +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 21,
      size = 1.25,
      color = "black",
      fill = "white"
    ) +
    scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
    theme_minimal(base_size = 6) +
    theme(legend.position = "none") +
    labs(x = NULL, y = "Log loss")
export::graph2svg(
  file = "figure2_plot_Q1",
  width = 1.8,
  height = 1.75,
  font = "Arial",
  bg = "transparent",
)
result %>%
  dplyr::filter(scenario %in% "Q2") %>%
  ggplot2::ggplot(mapping = aes(x = method, y = log_loss, color = method)) +
    ylim(0, max(result$log_loss)) +
    geom_jitter(alpha = 0.25, size = 0.75) +
    geom_violin(fill = "black", color = "black") +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 21,
      size = 1.25,
      color = "black",
      fill = "white"
    ) +
    scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
    theme_minimal(base_size = 6) +
    theme(legend.position = "none") +
    labs(x = NULL, y = "Log loss")
export::graph2svg(
  file = "figure2_plot_Q2",
  width = 1.8,
  height = 1.75,
  font = "Arial",
  bg = "transparent",
)
result %>%
  dplyr::filter(scenario %in% "Q3") %>%
  ggplot2::ggplot(mapping = aes(x = method, y = log_loss, color = method)) +
    ylim(0, max(result$log_loss)) +
    geom_jitter(alpha = 0.25, size = 0.75) +
    geom_violin(fill = "black", color = "black") +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 21,
      size = 1.25,
      color = "black",
      fill = "white"
    ) +
    scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
    theme_minimal(base_size = 6) +
    theme(legend.position = "none") +
    labs(x = NULL, y = "Log loss")
export::graph2svg(
  file = "figure2_plot_Q3",
  width = 1.8,
  height = 1.75,
  font = "Arial",
  bg = "transparent",
)
result %>%
  dplyr::filter(scenario %in% "Q4") %>%
  ggplot2::ggplot(mapping = aes(x = method, y = log_loss, color = method)) +
    ylim(0, max(result$log_loss)) +
    geom_jitter(alpha = 0.25, size = 0.75) +
    geom_violin(fill = "black", color = "black") +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 21,
      size = 1.25,
      color = "black",
      fill = "white"
    ) +
    scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
    theme_minimal(base_size = 6) +
    theme(legend.position = "none") +
    labs(x = NULL, y = "Log loss")
export::graph2svg(
  file = "figure2_plot_Q4",
  width = 1.8,
  height = 1.75,
  font = "Arial",
  bg = "transparent",
)
result %>%
  dplyr::filter(scenario %in% "Q5") %>%
  ggplot2::ggplot(mapping = aes(x = method, y = log_loss, color = method)) +
    ylim(0, max(result$log_loss)) +
    geom_jitter(alpha = 0.25, size = 0.75) +
    geom_violin(fill = "black", color = "black") +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 21,
      size = 1.25,
      color = "black",
      fill = "white"
    ) +
    scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
    theme_minimal(base_size = 6) +
    theme(legend.position = "none") +
    labs(x = NULL, y = "Log loss")
export::graph2svg(
  file = "figure2_plot_Q5",
  width = 1.8,
  height = 1.75,
  font = "Arial",
  bg = "transparent",
)
result %>%
  dplyr::filter(scenario %in% "Q6") %>%
  ggplot2::ggplot(mapping = aes(x = method, y = log_loss, color = method)) +
    ylim(0, max(result$log_loss)) +
    geom_jitter(alpha = 0.25, size = 0.75) +
    geom_violin(fill = "black", color = "black") +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 21,
      size = 1.25,
      color = "black",
      fill = "white"
    ) +
    scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
    theme_minimal(base_size = 6) +
    theme(legend.position = "none") +
    labs(x = NULL, y = "Log loss")
export::graph2svg(
  file = "figure2_plot_Q6",
  width = 1.8,
  height = 1.75,
  font = "Arial",
  bg = "transparent",
)
result %>%
  dplyr::filter(scenario %in% "Q7") %>%
  ggplot2::ggplot(mapping = aes(x = method, y = log_loss, color = method)) +
    ylim(0, max(result$log_loss)) +
    geom_jitter(alpha = 0.25, size = 0.75) +
    geom_violin(fill = "black", color = "black") +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 21,
      size = 1.25,
      color = "black",
      fill = "white"
    ) +
    scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
    theme_minimal(base_size = 6) +
    theme(legend.position = "none") +
    labs(x = NULL, y = "Log loss")
export::graph2svg(
  file = "figure2_plot_Q7",
  width = 1.8,
  height = 1.75,
  font = "Arial",
  bg = "transparent",
)
result %>%
  dplyr::filter(scenario %in% "Q8") %>%
  ggplot2::ggplot(mapping = aes(x = method, y = log_loss, color = method)) +
    ylim(0, max(result$log_loss)) +
    geom_jitter(alpha = 0.25, size = 0.75) +
    geom_violin(fill = "black", color = "black") +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 21,
      size = 1.25,
      color = "black",
      fill = "white"
    ) +
    scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
    theme_minimal(base_size = 6) +
    theme(legend.position = "none") +
    labs(x = NULL, y = "Log loss")
export::graph2svg(
  file = "figure2_plot_Q8",
  width = 1.8,
  height = 1.75,
  font = "Arial",
  bg = "transparent",
)
result %>%
  dplyr::filter(scenario %in% "Q9") %>%
  ggplot2::ggplot(mapping = aes(x = method, y = log_loss, color = method)) +
    ylim(0, max(result$log_loss)) +
    geom_jitter(alpha = 0.25, size = 0.75) +
    geom_violin(fill = "black", color = "black") +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 21,
      size = 1.25,
      color = "black",
      fill = "white"
    ) +
    scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
    theme_minimal(base_size = 6) +
    theme(legend.position = "none") +
    labs(x = NULL, y = "Log loss")
export::graph2svg(
  file = "figure2_plot_Q9",
  width = 1.8,
  height = 1.75,
  font = "Arial",
  bg = "transparent",
)
result %>%
  dplyr::filter(scenario %in% "Q10") %>%
  ggplot2::ggplot(mapping = aes(x = method, y = log_loss, color = method)) +
    ylim(0, max(result$log_loss)) +
    geom_jitter(alpha = 0.25, size = 0.75) +
    geom_violin(fill = "black", color = "black") +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 21,
      size = 1.25,
      color = "black",
      fill = "white"
    ) +
    scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
    theme_minimal(base_size = 6) +
    theme(legend.position = "none") +
    labs(x = NULL, y = "Log loss")
export::graph2svg(
  file = "figure2_plot_Q10",
  width = 1.8,
  height = 1.75,
  font = "Arial",
  bg = "transparent",
)
result %>%
  dplyr::filter(scenario %in% "Q11") %>%
  ggplot2::ggplot(mapping = aes(x = method, y = log_loss, color = method)) +
    ylim(0, max(result$log_loss)) +
    geom_jitter(alpha = 0.25, size = 0.75) +
    geom_violin(fill = "black", color = "black") +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 21,
      size = 1.25,
      color = "black",
      fill = "white"
    ) +
    scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
    theme_minimal(base_size = 6) +
    theme(legend.position = "none") +
    labs(x = NULL, y = "Log loss")
export::graph2svg(
  file = "figure2_plot_Q11",
  width = 1.8,
  height = 1.75,
  font = "Arial",
  bg = "transparent",
)
result %>%
  dplyr::filter(scenario %in% "Q12") %>%
  ggplot2::ggplot(mapping = aes(x = method, y = log_loss, color = method)) +
    ylim(0, max(result$log_loss)) +
    geom_jitter(alpha = 0.25, size = 0.75) +
    geom_violin(fill = "black", color = "black") +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 21,
      size = 1.25,
      color = "black",
      fill = "white"
    ) +
    scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
    theme_minimal(base_size = 6) +
    theme(legend.position = "none") +
    labs(x = NULL, y = "Log loss")
export::graph2svg(
  file = "figure2_plot_Q12",
  width = 1.8,
  height = 1.75,
  font = "Arial",
  bg = "transparent",
)
