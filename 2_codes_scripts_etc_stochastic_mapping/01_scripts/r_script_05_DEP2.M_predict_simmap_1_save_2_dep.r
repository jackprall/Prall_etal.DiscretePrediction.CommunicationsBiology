# Load packages. ----
library(doSNOW)
library(phytools)

# Load custom functions. ----
source("functions.r")

# Read trees and trait datasets. ----
load(file = "trees_train_DEP2.M.RData")
load(file = "traits_train_DEP2.M.RData")

# Specify settings. ----
n_trials <- 1000
n_cores <- 16

# Specify constant inputs for the `make.simmap` function. ----
Q_dep <- matrix(
  data = c(
    0, 1, 2, 0,
    3, 0, 0, 4,
    5, 0, 0, 6,
    0, 7, 8, 0
  ),
  nrow = 4,
  ncol = 4,
  byrow = TRUE
)
rownames(Q_dep) <- colnames(Q_dep) <- c("00", "01", "10", "11")

# Fit a phylogenetic continuous-time Markov chain (CTMC) model on each of the
#   1,000 training datasets. ----
cluster <- snow::makeCluster(spec = n_cores, type = "SOCK", outfile = "")
doSNOW::registerDoSNOW(cl = cluster)
set.seed(2502)
simmap_posteriors_DEP2.M_2_dep <-
  foreach(
    i = iterators::icount(n_trials)
  ) %dopar% {
  trait_train <- traits_train[[i]]
  trait_train$xy <- as.factor(trait_train$xy)
  levels(trait_train$xy) <- c("00", "01", "10", "11")
  trait_train <- model.matrix(~trait_train$xy - 1)
  rownames(trait_train) <- traits_train[[i]]$taxon
  colnames(trait_train) <- c("00", "01", "10", "11")
  simmap_posterior_DEP2.M_2_dep <- phytools::make.simmap(
    tree = trees_train[[i]],
    x = trait_train,
    model = Q_dep,
    Q = "mcmc",
    pi = "equal",
    prior = list(alpha = 1.0, beta = 1.0),  # Gamma(1, 1) = Exponential(1)
    nsim = 400,
    burnin = 8000,
    samplefreq = 100,
    message = FALSE
  )
  print(paste0("finished trial ", i, "..."))
  return(simmap_posterior_DEP2.M_2_dep)
}
snow::stopCluster(cl = cluster)
save(
  simmap_posteriors_DEP2.M_2_dep,
  file = "simmap_posteriors_DEP2.M_2_dep.RData"
)
