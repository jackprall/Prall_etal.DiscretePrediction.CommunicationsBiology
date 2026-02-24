# Load packages. ----
library(doSNOW)
library(phytools)

# Load custom functions. ----
source("functions.r")

# Read trees and trait datasets. ----
load(file = "trees.RData")
load(file = "test_taxa_random.RData")
load(file = "traits_test_random.RData")
load(file = "traits_train_random.RData")

# Specify settings. ----
n_trials <- 1000
n_cores <- 16

# Prepare the trait matrix, adding the test taxon, and specify flat prior
#   probabilities for Y = 0 and Y = 1. ----
simmap_traits <- vector(mode = "list", length = n_trials)
for (i in 1:n_trials) {
  trait <- rbind(traits_test[[i]], traits_train[[i]])
  trait$xy <- as.factor(trait$xy)
  levels(trait$xy) <- c("00", "01", "10", "11")
  temp <- model.matrix(~trait$xy - 1)
  rownames(temp) <- trait$taxon
  colnames(temp) <- c("00", "01", "10", "11")
  if (traits_test[[i]]$x == 0) {
    temp[rownames(temp) %in% test_taxa[i], ] <- c(0.5, 0.5, 0, 0)
  } else {
    temp[rownames(temp) %in% test_taxa[i], ] <- c(0, 0, 0.5, 0.5)
  }
  simmap_traits[[i]] <- temp
}

# Read the posterior distributions of the Q matrix for each of the 1,000
#   trials. ----
load(file = "simmap_posteriors_random_2_dep.RData")
post_Q <- vector(mode = "list", length = n_trials)
for (i in 1:n_trials) {
  post_Q[[i]] <- vector(
    mode = "list",
    length = length(simmap_posteriors_random_2_dep[[i]])
  )
  for (j in 1:length(post_Q[[i]])) {
    post_Q[[i]][[j]] <- simmap_posteriors_random_2_dep[[i]][[j]]$Q
  }
}

# Generate a posterior predictive distribution for the test taxon in each of
#   the 1,000 trials. ----
cluster <- snow::makeCluster(spec = n_cores, type = "SOCK", outfile = "")
doSNOW::registerDoSNOW(cl = cluster)
set.seed(2504)
simmap_ppds_random_2_dep <-
  foreach(
    i = iterators::icount(n_trials)
  ) %dopar% {
  simmap_ppd_random_2_dep <- post_pred_simmap(
    tree = trees[[i]],
    trait = simmap_traits[[i]],
    test_taxon = test_taxa[i],
    post_Q = post_Q[[i]]
  )
  print(paste0("finished trial ", i, "..."))
  return(simmap_ppd_random_2_dep)
}
snow::stopCluster(cl = cluster)
save(simmap_ppds_random_2_dep, file = "simmap_ppds_random_2_dep.RData")
