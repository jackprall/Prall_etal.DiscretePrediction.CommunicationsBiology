# Tutorial.

# Load packages. ----
library(ape)
library(phytools)
library(tidyverse)

# Load the phylogenetic tree. ----
tree <- ape::read.nexus(file = "tutorial_tree.nex")

# Load the trait dataset. ----
traits <- readr::read_tsv(file = "tutorial_traits.tsv", col_names = FALSE)
colnames(traits) <- c("tip", "y", "x")

# Specify the test tip (i.e., the tip to be predicted). ----
test_tip <- "t100"

# Create the two trait files for `BayesTraits`. ----
traits |>
  dplyr::mutate(
    y = ifelse(tip == test_tip, "-", as.character(y)),
    x = ifelse(tip == test_tip, "-", as.character(x))
  ) |>
  readr::write_tsv(file = "tutorial_traits_save.tsv", col_names = FALSE)
traits |>
  dplyr::mutate(y = ifelse(tip == test_tip, "?", as.character(y))) |>
  readr::write_tsv(file = "tutorial_traits_load.tsv", col_names = FALSE)
