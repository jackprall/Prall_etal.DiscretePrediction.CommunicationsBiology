# Load packages. ----
library(castor)
library(miceadds)
library(phytools)
library(tidyverse)

# Load custom functions. ----
source("functions.r")

# Load tree, trait data, and estimated transition rates. ----
setwd("inputs")
load("trees.RData")
trees_full_list <- trees
rm(trees)
scenarios <- c("ER.L", "ER.M", "ER.H", "DR.LM", "DR.LH", "DR.MH", "DEP1.L",
               "DEP1.M", "DEP1.H", "DEP2.L", "DEP2.M", "DEP2.H", "random")
trees_train_list <- setNames(
  object = vector(mode = "list", length = length(scenarios)),
  nm = scenarios
)
for (i in 1:length(scenarios)) {
  miceadds::load.Rdata(
    filename = paste0("trees_train_", scenarios[i], ".RData"),
    objname = "temp"
  )
  trees_train_list[[i]] <- temp
}
load("traits_test_list.RData")
for (i in 1:length(scenarios)) {
  traits_test_list[[i]] <- purrr::list_rbind(traits_test_list[[i]])
  traits_test_list[[i]]$tip_states <- dplyr::case_when(
    traits_test_list[[i]]$xy == "00" ~ 1,
    traits_test_list[[i]]$xy == "01" ~ 2,
    traits_test_list[[i]]$xy == "10" ~ 3,
    traits_test_list[[i]]$xy == "11" ~ 4
  )
}
n_trials <- 1000
traits_train_list <- setNames(
  object = vector(mode = "list", length = length(scenarios)),
  nm = scenarios
)
for (i in 1:length(scenarios)) {
  miceadds::load.Rdata(
    filename = paste0("traits_train_", scenarios[i], ".RData"),
    objname = "temp"
  )
  for (j in 1:n_trials) {
    temp2 <- temp[[j]]
    temp2 <- temp2[match(trees_train_list[[i]][[j]]$tip.label, temp2$taxon), ]
    temp2$tip_states <- dplyr::case_when(
      temp2$xy == "00" ~ 1,
      temp2$xy == "01" ~ 2,
      temp2$xy == "10" ~ 3,
      temp2$xy == "11" ~ 4
    )
    temp2$tip_states <- as.factor(temp2$tip_states)
    temp2 <- temp2$tip_states
    temp[[j]] <- temp2
  }
  traits_train_list[[i]] <- temp
}
setwd("../outputs")
load("qmat_median_2_dep.RData")
q_final <- q_median_dep
for (i in 1:length(scenarios)) {
  for (j in 1:n_trials) {
    temp <- q_median_dep[[i]][[j]]
    diag(temp) <- 0
    diag(temp) <- -rowSums(temp)
    q_final[[i]][[j]] <- temp
  }
}

# Calculate scaled conditional likelihoods of the states 00, 01, 10, and 11 at
#   the grandparent node of the test taxon (SIMMAP: Dependent). ----
temp_results <- vector(mode = "list", length = length(scenarios))
for (i in 1:length(scenarios)) {
  temp_results[[i]] <- vector(mode = "list", length = n_trials)
}
for (i in 1:length(scenarios)) {
  for (j in 1:n_trials) {
    print(paste0("starting scenario ", scenarios[i], ", trial ", j))
    grandparent_node <- find_grandparent_node(
      tree_full = trees_full_list[[j]],
      tree_train = trees_train_list[[i]][[j]],
      test_tip = traits_test_list[[i]][j, "taxon"]
    )
    scaled_condlh <- get_scaled_condlh(
      tree_train = trees_train_list[[i]][[j]],
      tip_states = as.integer(traits_train_list[[i]][[j]]),
      n_states = 4,
      q_median = q_final[[i]][[j]],
      node = grandparent_node
    )
    temp_results[[i]][[j]] <- scaled_condlh
  }
}

# Calculate the difference in scaled conditional likelihoods at the state that
#   the test taxon has. ----
results <- vector(mode = "list", length = length(scenarios))
for (i in 1:length(scenarios)) {
  results[[i]] <- matrix(
    data = NA,
    nrow = n_trials,
    ncol = 2,
    dimnames = list(1:n_trials, c("anc.lh", "anc.diff"))
  )
}
for (i in 1:length(scenarios)) {
  for (j in 1:n_trials) {
    results[[i]][j, "anc.lh"] <-
      temp_results[[i]][[j]][as.numeric(traits_test_list[[i]][j, "tip_states"])]
    results[[i]][j, "anc.diff"] <- 1 - results[[i]][j, "anc.lh"]
  }
}
for (i in 1:length(scenarios)) {
  results[[i]] |>
    as.data.frame() |>
    readr::write_csv(
      file = paste0("temp_anc_diff_", scenarios[i], ".csv")
    )
}
