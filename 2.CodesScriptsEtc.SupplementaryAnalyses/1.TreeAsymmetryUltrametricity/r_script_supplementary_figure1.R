# Load packages. ----
library(export)
library(mosaic)
library(readxl)
library(tidyverse)

# Read and clean results. ----
result <- readxl::read_xlsx(
  path = "pilot_tree_asymmetry_ultrametricity.xlsx",
  sheet = "Results"
)
result <- result |>
  dplyr::select(Comb, P.Correct) |>
  dplyr::mutate(Comb = as.factor(Comb)) |>
  dplyr::mutate(
    Comb = forcats::fct_recode(
      Comb,
      "Non-ultrametric\nAsymmetric" = "NA",
      "Non-ultrametric\nSymmetric" = "NS",
      "Ultrametric\nAsymmetric" = "UA",
      "Ultrametric\nSymmetric" = "US"
    )
  )

# Generate summary statistics. ----
summ <- mosaic::favstats(P.Correct ~ Comb, data = result)
summ$IQR <- summ$Q3 - summ$Q1
summ
#>                          Comb   min    Q1 median    Q3 max     mean        sd    n missing   IQR
#> 1 Non-ultrametric\nAsymmetric 0.004 0.969  0.986 0.994   1 0.952353 0.1440105 1000       0 0.025
#> 2  Non-ultrametric\nSymmetric 0.007 0.970  0.985 0.994   1 0.957262 0.1275071 1000       0 0.024
#> 3     Ultrametric\nAsymmetric 0.004 0.973  0.989 0.996   1 0.959767 0.1276653 1000       0 0.023
#> 4      Ultrametric\nSymmetric 0.013 0.974  0.987 0.996   1 0.955067 0.1502178 1000       0 0.022

# Create the plot. ----
result %>%
  ggplot2::ggplot(mapping = aes(x = Comb, y = P.Correct, color = Comb)) +
    geom_violin() +
    geom_jitter(alpha = 0.05, width = 0.25) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(
      x = "\nPhylogenetic tree",
      y = "Probability of a correct prediction\n"
    )
export::graph2pdf(
  file = "Prall_etal.ComBio.DiscretePrediction.SupplementaryFigure1.pdf",
  width = 6,
  height = 6 / 1.618,
  font = "Arial",
  bg = "transparent",
)
export::graph2png(
  file = "Prall_etal.ComBio.DiscretePrediction.SupplementaryFigure1.png",
  width = 6,
  height = 6 / 1.618,
  font = "Arial",
  bg = "transparent",
)
export::graph2svg(
  file = "Prall_etal.ComBio.DiscretePrediction.SupplementaryFigure1.svg",
  width = 6,
  height = 6 / 1.618,
  font = "Arial",
  bg = "transparent",
)
