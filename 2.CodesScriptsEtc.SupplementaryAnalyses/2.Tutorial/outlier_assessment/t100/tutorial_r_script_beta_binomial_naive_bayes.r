# Tutorial.

# Load packages. ----
library(mosaic)
library(tidyverse)

# Load and prepare the trait dataset. ----
traits <- readr::read_tsv(
  file = "tutorial_edited_traits.tsv",
  col_names = FALSE
)
colnames(traits) <- c("tip", "x", "y")
traits_train <- traits |>
  dplyr::filter(tip != "t100")

# Predict using the Bayesian Beta-Binomial method. ----
alpha <- 1 + sum(traits_train$y)
beta <- 1 + length(traits_train$y) - sum(traits_train$y)
posterior_samp <- rbeta(n = 1e+05, shape1 = alpha, shape2 = beta)
posterior_pred <- rbinom(n = 1e+05, size = 1, prob = posterior_samp)
point_pred_prob <- mean(posterior_pred)
point_pred_prob
#> [1] 0.64703

# Naive Bayes. =================================================================
conf_mat <- tally(traits_train$y ~ traits_train$x)
prior <- (conf_mat[2] + conf_mat[4]) / sum(conf_mat)
likelihood <- conf_mat[4] / (conf_mat[2] + conf_mat[4])
marginal <- (conf_mat[3] + conf_mat[4]) / sum(conf_mat)
posterior <- prior * likelihood / marginal
posterior
#> [1] 0.9661017
