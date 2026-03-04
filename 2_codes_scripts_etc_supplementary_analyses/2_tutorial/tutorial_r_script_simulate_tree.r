# Tutorial.

# Load packages. ----
library(ape)
library(phytools)

# Simulate a phylogenetic tree. ----
set.seed(1)
tree <- ape::rphylo(n = 500, birth = 1.5, death = 0.5)
tree
#>
#> Phylogenetic tree with 500 tips and 499 internal nodes.
#>
#> Tip labels:
#>   t1, t2, t3, t4, t5, t6, ...
#>
#> Rooted; includes branch lengths.
tree$edge.length <- tree$edge.length / sum(tree$edge.length) * 700
sum(tree$edge.length)
#> [1] 700

# Save the simulated phylogenetic tree. ----
phytools::writeNexus(tree = tree, file = "tutorial_tree.nex")
