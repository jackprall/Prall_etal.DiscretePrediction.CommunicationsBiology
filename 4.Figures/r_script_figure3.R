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

# Process results for Figure 3. ----
result <- result |>
  dplyr::select(scenario, xy, is_max, terminal_std, log_loss_dep_rj) |>
  dplyr::mutate(is_max = forcats::fct_recode(is_max, "No" = "0", "Yes" = "1"))

# Create the individual plots. ----
result %>%
  dplyr::filter(scenario %in% "Random") %>%
  ggplot2::ggplot(
    mapping = aes(x = terminal_std, y = log_loss_dep_rj, color = is_max)
  ) +
    xlim(0, 1) +
    ylim(0, max(result$log_loss_dep_rj)) +
    geom_hline(yintercept = 0.693, color = "black", linetype = "dotted") +
    geom_point(alpha = 0.25, size = 0.75) +
    geom_smooth(method = "lm", formula = y ~ log(x), se = FALSE) +
    theme_minimal(base_size = 6) +
    theme(panel.grid.minor = element_blank()) +
    labs(
      x = "Standardized terminal branch length\n(terminal branch length / tree height)",
      y = "Log loss",
      color = "\"Match\" with sister"
    )
export::graph2svg(
  file = "figure3_plot_random_temp",
  width = 1.8,
  height = 1.75,
  font = "Arial",
  bg = "transparent",
)
result %>%
  dplyr::filter(scenario %in% "Random") %>%
  ggplot2::ggplot(
    mapping = aes(x = terminal_std, y = log_loss_dep_rj, color = is_max)
  ) +
    xlim(0, 1) +
    ylim(0, max(result$log_loss_dep_rj)) +
    geom_hline(yintercept = 0.693, color = "black", linetype = "dotted") +
    geom_point(alpha = 0.25, size = 0.75) +
    geom_smooth(method = "lm", formula = y ~ log(x), se = FALSE) +
    theme_minimal(base_size = 6) +
    theme(legend.position = "none", panel.grid.minor = element_blank()) +
    labs(
      x = "Standardized terminal branch length\n(terminal branch length / tree height)",
      y = "Log loss",
      color = "\"Match\" with sister"
    )
export::graph2svg(
  file = "figure3_plot_random",
  width = 1.8,
  height = 1.75,
  font = "Arial",
  bg = "transparent",
)
result %>%
  dplyr::filter(scenario %in% "Q1") %>%
  ggplot2::ggplot(
    mapping = aes(x = terminal_std, y = log_loss_dep_rj, color = is_max)
  ) +
    xlim(0, 1) +
    ylim(0, max(result$log_loss_dep_rj)) +
    geom_hline(yintercept = 0.693, color = "black", linetype = "dotted") +
    geom_point(alpha = 0.25, size = 0.75) +
    geom_smooth(method = "lm", formula = y ~ log(x), se = FALSE) +
    theme_minimal(base_size = 6) +
    theme(legend.position = "none", panel.grid.minor = element_blank()) +
    labs(
      x = "Standardized terminal branch length\n(terminal branch length / tree height)",
      y = "Log loss",
      color = "\"Match\" with sister"
    )
export::graph2svg(
  file = "figure3_plot_Q1",
  width = 1.8,
  height = 1.75,
  font = "Arial",
  bg = "transparent",
)
result %>%
  dplyr::filter(scenario %in% "Q2") %>%
  ggplot2::ggplot(
    mapping = aes(x = terminal_std, y = log_loss_dep_rj, color = is_max)
  ) +
    xlim(0, 1) +
    ylim(0, max(result$log_loss_dep_rj)) +
    geom_hline(yintercept = 0.693, color = "black", linetype = "dotted") +
    geom_point(alpha = 0.25, size = 0.75) +
    geom_smooth(method = "lm", formula = y ~ log(x), se = FALSE) +
    theme_minimal(base_size = 6) +
    theme(legend.position = "none", panel.grid.minor = element_blank()) +
    labs(
      x = "Standardized terminal branch length\n(terminal branch length / tree height)",
      y = "Log loss",
      color = "\"Match\" with sister"
    )
export::graph2svg(
  file = "figure3_plot_Q2",
  width = 1.8,
  height = 1.75,
  font = "Arial",
  bg = "transparent",
)
result %>%
  dplyr::filter(scenario %in% "Q3") %>%
  ggplot2::ggplot(
    mapping = aes(x = terminal_std, y = log_loss_dep_rj, color = is_max)
  ) +
    xlim(0, 1) +
    ylim(0, max(result$log_loss_dep_rj)) +
    geom_hline(yintercept = 0.693, color = "black", linetype = "dotted") +
    geom_point(alpha = 0.25, size = 0.75) +
    geom_smooth(method = "lm", formula = y ~ log(x), se = FALSE) +
    theme_minimal(base_size = 6) +
    theme(legend.position = "none", panel.grid.minor = element_blank()) +
    labs(
      x = "Standardized terminal branch length\n(terminal branch length / tree height)",
      y = "Log loss",
      color = "\"Match\" with sister"
    )
export::graph2svg(
  file = "figure3_plot_Q3",
  width = 1.8,
  height = 1.75,
  font = "Arial",
  bg = "transparent",
)
result %>%
  dplyr::filter(scenario %in% "Q4") %>%
  ggplot2::ggplot(
    mapping = aes(x = terminal_std, y = log_loss_dep_rj, color = is_max)
  ) +
    xlim(0, 1) +
    ylim(0, max(result$log_loss_dep_rj)) +
    geom_hline(yintercept = 0.693, color = "black", linetype = "dotted") +
    geom_point(alpha = 0.25, size = 0.75) +
    geom_smooth(method = "lm", formula = y ~ log(x), se = FALSE) +
    theme_minimal(base_size = 6) +
    theme(legend.position = "none", panel.grid.minor = element_blank()) +
    labs(
      x = "Standardized terminal branch length\n(terminal branch length / tree height)",
      y = "Log loss",
      color = "\"Match\" with sister"
    )
export::graph2svg(
  file = "figure3_plot_Q4",
  width = 1.8,
  height = 1.75,
  font = "Arial",
  bg = "transparent",
)
result %>%
  dplyr::filter(scenario %in% "Q5") %>%
  ggplot2::ggplot(
    mapping = aes(x = terminal_std, y = log_loss_dep_rj, color = is_max)
  ) +
    xlim(0, 1) +
    ylim(0, max(result$log_loss_dep_rj)) +
    geom_hline(yintercept = 0.693, color = "black", linetype = "dotted") +
    geom_point(alpha = 0.25, size = 0.75) +
    geom_smooth(method = "lm", formula = y ~ log(x), se = FALSE) +
    theme_minimal(base_size = 6) +
    theme(legend.position = "none", panel.grid.minor = element_blank()) +
    labs(
      x = "Standardized terminal branch length\n(terminal branch length / tree height)",
      y = "Log loss",
      color = "\"Match\" with sister"
    )
export::graph2svg(
  file = "figure3_plot_Q5",
  width = 1.8,
  height = 1.75,
  font = "Arial",
  bg = "transparent",
)
result %>%
  dplyr::filter(scenario %in% "Q6") %>%
  ggplot2::ggplot(
    mapping = aes(x = terminal_std, y = log_loss_dep_rj, color = is_max)
  ) +
    xlim(0, 1) +
    ylim(0, max(result$log_loss_dep_rj)) +
    geom_hline(yintercept = 0.693, color = "black", linetype = "dotted") +
    geom_point(alpha = 0.25, size = 0.75) +
    geom_smooth(method = "lm", formula = y ~ log(x+0.000001), se = FALSE) +
    theme_minimal(base_size = 6) +
    theme(legend.position = "none", panel.grid.minor = element_blank()) +
    labs(
      x = "Standardized terminal branch length\n(terminal branch length / tree height)",
      y = "Log loss",
      color = "\"Match\" with sister"
    )
export::graph2svg(
  file = "figure3_plot_Q6",
  width = 1.8,
  height = 1.75,
  font = "Arial",
  bg = "transparent",
)
result %>%
  dplyr::filter(scenario %in% "Q7") %>%
  ggplot2::ggplot(
    mapping = aes(x = terminal_std, y = log_loss_dep_rj, color = is_max)
  ) +
    xlim(0, 1) +
    ylim(0, max(result$log_loss_dep_rj)) +
    geom_hline(yintercept = 0.693, color = "black", linetype = "dotted") +
    geom_point(alpha = 0.25, size = 0.75) +
    geom_smooth(method = "lm", formula = y ~ log(x), se = FALSE) +
    theme_minimal(base_size = 6) +
    theme(legend.position = "none", panel.grid.minor = element_blank()) +
    labs(
      x = "Standardized terminal branch length\n(terminal branch length / tree height)",
      y = "Log loss",
      color = "\"Match\" with sister"
    )
export::graph2svg(
  file = "figure3_plot_Q7",
  width = 1.8,
  height = 1.75,
  font = "Arial",
  bg = "transparent",
)
result %>%
  dplyr::filter(scenario %in% "Q8") %>%
  ggplot2::ggplot(
    mapping = aes(x = terminal_std, y = log_loss_dep_rj, color = is_max)
  ) +
    xlim(0, 1) +
    ylim(0, max(result$log_loss_dep_rj)) +
    geom_hline(yintercept = 0.693, color = "black", linetype = "dotted") +
    geom_point(alpha = 0.25, size = 0.75) +
    geom_smooth(method = "lm", formula = y ~ log(x), se = FALSE) +
    theme_minimal(base_size = 6) +
    theme(legend.position = "none", panel.grid.minor = element_blank()) +
    labs(
      x = "Standardized terminal branch length\n(terminal branch length / tree height)",
      y = "Log loss",
      color = "\"Match\" with sister"
    )
export::graph2svg(
  file = "figure3_plot_Q8",
  width = 1.8,
  height = 1.75,
  font = "Arial",
  bg = "transparent",
)
result %>%
  dplyr::filter(scenario %in% "Q9") %>%
  ggplot2::ggplot(
    mapping = aes(x = terminal_std, y = log_loss_dep_rj, color = is_max)
  ) +
    xlim(0, 1) +
    ylim(0, max(result$log_loss_dep_rj)) +
    geom_hline(yintercept = 0.693, color = "black", linetype = "dotted") +
    geom_point(alpha = 0.25, size = 0.75) +
    geom_smooth(method = "lm", formula = y ~ log(x), se = FALSE) +
    theme_minimal(base_size = 6) +
    theme(legend.position = "none", panel.grid.minor = element_blank()) +
    labs(
      x = "Standardized terminal branch length\n(terminal branch length / tree height)",
      y = "Log loss",
      color = "\"Match\" with sister"
    )
export::graph2svg(
  file = "figure3_plot_Q9",
  width = 1.8,
  height = 1.75,
  font = "Arial",
  bg = "transparent",
)
result %>%
  dplyr::filter(scenario %in% "Q10") %>%
  ggplot2::ggplot(
    mapping = aes(x = terminal_std, y = log_loss_dep_rj, color = is_max)
  ) +
    xlim(0, 1) +
    ylim(0, max(result$log_loss_dep_rj)) +
    geom_hline(yintercept = 0.693, color = "black", linetype = "dotted") +
    geom_point(alpha = 0.25, size = 0.75) +
    geom_smooth(method = "lm", formula = y ~ log(x), se = FALSE) +
    theme_minimal(base_size = 6) +
    theme(legend.position = "none", panel.grid.minor = element_blank()) +
    labs(
      x = "Standardized terminal branch length\n(terminal branch length / tree height)",
      y = "Log loss",
      color = "\"Match\" with sister"
    )
export::graph2svg(
  file = "figure3_plot_Q10",
  width = 1.8,
  height = 1.75,
  font = "Arial",
  bg = "transparent",
)
result %>%
  dplyr::filter(scenario %in% "Q11") %>%
  ggplot2::ggplot(
    mapping = aes(x = terminal_std, y = log_loss_dep_rj, color = is_max)
  ) +
    xlim(0, 1) +
    ylim(0, max(result$log_loss_dep_rj)) +
    geom_hline(yintercept = 0.693, color = "black", linetype = "dotted") +
    geom_point(alpha = 0.25, size = 0.75) +
    geom_smooth(method = "lm", formula = y ~ log(x), se = FALSE) +
    theme_minimal(base_size = 6) +
    theme(legend.position = "none", panel.grid.minor = element_blank()) +
    labs(
      x = "Standardized terminal branch length\n(terminal branch length / tree height)",
      y = "Log loss",
      color = "\"Match\" with sister"
    )
export::graph2svg(
  file = "figure3_plot_Q11",
  width = 1.8,
  height = 1.75,
  font = "Arial",
  bg = "transparent",
)
result %>%
  dplyr::filter(scenario %in% "Q12") %>%
  ggplot2::ggplot(
    mapping = aes(x = terminal_std, y = log_loss_dep_rj, color = is_max)
  ) +
    xlim(0, 1) +
    ylim(0, max(result$log_loss_dep_rj)) +
    geom_hline(yintercept = 0.693, color = "black", linetype = "dotted") +
    geom_point(alpha = 0.25, size = 0.75) +
    geom_smooth(method = "lm", formula = y ~ log(x), se = FALSE) +
    theme_minimal(base_size = 6) +
    theme(legend.position = "none", panel.grid.minor = element_blank()) +
    labs(
      x = "Standardized terminal branch length\n(terminal branch length / tree height)",
      y = "Log loss",
      color = "\"Match\" with sister"
    )
export::graph2svg(
  file = "figure3_plot_Q12",
  width = 1.8,
  height = 1.75,
  font = "Arial",
  bg = "transparent",
)
