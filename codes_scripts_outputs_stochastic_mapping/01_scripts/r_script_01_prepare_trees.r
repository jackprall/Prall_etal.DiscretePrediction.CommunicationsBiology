# Load packages. ----
library(ape)

# Prepare trees. ----
n_trees <- 1000
trees <- vector(mode = "list", length = n_trees)
trees_height <- vector(length = n_trees)
setwd("50Tips_tree")
for (i in 1:n_trees) {
  trees[[i]] <- ape::read.nexus(
    file = paste0("Full_tree.", i, ".tre")
  )
  vcv <- ape::vcv.phylo(phy = trees[[i]])
  trees_height[i] <- max(diag(vcv))
}
setwd("..")
save(trees, file = "trees.RData")
write.table(
  x = trees_height,
  file = "temp_trees_height.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
