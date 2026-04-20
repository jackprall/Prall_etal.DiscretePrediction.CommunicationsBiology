# Tutorial.

# Load packages. ----
library(ape)
library(phytools)
library(tidyverse)

# Load the phylogenetic tree. ----
tree <- ape::read.nexus(file = "tutorial_tree.nex")

# Specify the transition rate matrix. ----
a <- 0.65
b <- 0.05
Q <- matrix(
  data = c(
    0, a, b, 0,
    a, 0, 0, a,
    a, 0, 0, a,
    0, a, b, 0
  ),
  nrow = 4,
  ncol = 4,
  byrow = TRUE
)
rownames(Q) <- colnames(Q) <- c("00", "01", "10", "11")
diag(Q) <- -rowSums(Q)
Q
#>       00    01    10    11
#> 00 -0.70  0.65  0.05  0.00
#> 01  0.65 -1.30  0.00  0.65
#> 10  0.65  0.00 -1.30  0.65
#> 11  0.00  0.65  0.05 -0.70

# Simulate traits X and Y on the phylogenetic tree. ----
set.seed(2)
traits_temp <- phytools::sim.history(
  tree = tree,
  Q = Q,
  direction = "row_to_column"
)
x <- as.numeric(substr(x = traits_temp$states, start = 1, stop = 1))
y <- as.numeric(substr(x = traits_temp$states, start = 2, stop = 2))
traits <- data.frame(
  tip = tree$tip.label,
  x = x,
  y = y
)
head(traits)
#>    tip x y
#> 1   t1 0 1
#> 2  t10 0 0
#> 3 t100 1 1
#> 4 t101 1 1
#> 5 t102 1 1
#> 6 t103 0 0

# Save the simulated traits. ----
readr::write_tsv(x = traits, file = "tutorial_traits.tsv", col_names = FALSE)
