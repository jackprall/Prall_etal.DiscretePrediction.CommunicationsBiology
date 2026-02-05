# Load packages. ----
library(ape)

# Prepare training trees. ----
# ER.L - DEP2.H.
scenarios <- c("DEP1.H", "DEP1.L", "DEP1.M", "DEP2.H", "DEP2.L", "DEP2.M",
               "DR.LH", "DR.LM", "DR.MH", "ER.H", "ER.L", "ER.M")
load(file = "trees.RData")
n_trees <- 1000
for (i in scenarios) {
  load(file = paste0("test_taxa_", i, ".RData"))
  trees_train <- vector(mode = "list", length = n_trees)
  for (j in 1:n_trees) {
    trees_train[[j]] <- ape::drop.tip(phy = trees[[j]], tip = test_taxa[j])
  }
  save(trees_train, file = paste0("trees_train_", i,".RData"))
}
# Random.
load(file = "test_taxa_random.RData")
trees_train <- vector(mode = "list", length = n_trees)
for (j in 1:n_trees) {
  trees_train[[j]] <- ape::drop.tip(phy = trees[[j]], tip = test_taxa[j])
}
save(trees_train, file = "trees_train_random.RData")
