# Load packages. ----
library(tidyverse)

# Save BayesTraits posteriors. ----
scenarios <- c("ER.L", "ER.M", "ER.H", "DR.LM", "DR.LH", "DR.MH", "DEP1.L",
               "DEP1.M", "DEP1.H", "DEP2.L", "DEP2.M", "DEP2.H", "random")
n_trials <- 1000
n_samps <- 1001
# Dependent.
for (i in 1:length(scenarios)) {
  Q_post <- vector(mode = "list", length = n_trials)
  for (j in 1:n_trials) {
    Q_post[[j]] <- vector(mode = "list", length = n_samps)
  }
  setwd("posterior_bayestraits")
  setwd(scenarios[i])
  for (j in 1:n_trials) {
    post <- readr::read_tsv(
      file = paste0(scenarios[i], ".Dep.MCMC.", j, ".Rates.Log.txt"),
      skip = 53,
      name_repair = "minimal",
      show_col_types = FALSE
    )
    for (k in 1:n_samps) {
      mat <- matrix(
        data = c(
          0, post$q12[k], post$q13[k], 0,
          post$q21[k], 0, 0, post$q24[k],
          post$q31[k], 0, 0, post$q34[k],
          0, post$q42[k], post$q43[k], 0
        ),
        nrow = 4,
        ncol = 4,
        byrow = TRUE
      )
      diag(mat) <- -rowSums(mat)
      print(mat)
      Q_post[[j]][[k]] <- mat
    }
  }
  setwd("../..")
  save(
    Q_post,
    file = paste0("bayestraits_posteriors_", scenarios[i], "_2_dep.RData")
  )
}
