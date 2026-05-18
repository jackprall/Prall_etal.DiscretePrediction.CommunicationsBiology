# Load packages. ----
library(tidyverse)

# Process SIMMAP predictions. ----
source("functions.r")
load(file = "simmap_ppds_DEP2.H_2_dep.RData")
load(file = "simmap_ppds_DEP2.L_2_dep.RData")
load(file = "simmap_ppds_DEP2.M_2_dep.RData")
load(file = "simmap_ppds_DR.LH_2_dep.RData")
ppds_simmap_list_dep <- list(
  simmap_ppds_DR.LH_2_dep,
  simmap_ppds_DEP2.L_2_dep,
  simmap_ppds_DEP2.M_2_dep,
  simmap_ppds_DEP2.H_2_dep
)
scenarios <- c("DR.LH", "DEP2.L", "DEP2.M", "DEP2.H")
traits_test_list <- setNames(
  object = vector(mode = "list", length = length(scenarios)),
  nm = scenarios
)
for (i in scenarios) {
  load(file = paste0("traits_test_", i, ".RData"))
  traits_test_list[[i]] <- traits_test
  rm(traits_test)
}
names(traits_test_list) <- NULL
save(ppds_simmap_list_dep, file = "ppds_simmap_list_2_dep.RData")
save(traits_test_list, file = "traits_test_list.RData")
n_trials <- 1000
results <- vector(mode = "list", length = length(ppds_simmap_list_dep))
for (i in 1:length(scenarios)) {
  results[[i]] <- matrix(
    data = NA,
    nrow = n_trials,
    ncol = 3,
    dimnames = list(1:n_trials, c("p.dep.s2", "c.dep.s2", "ll.dep.s2"))
  )
}
for (i in 1:length(scenarios)) {
  traits_test <- purrr::list_rbind(traits_test_list[[i]])
  ppds_dep <- ppds_simmap_list_dep[[i]]
  pred_probs_dep <- as.numeric(unlist(lapply(ppds_dep, "mean")))
  results[[i]][, 1] <- pred_probs_dep
  preds_dep <- ifelse(pred_probs_dep > 0.5, 1, 0)
  c_dep <- ifelse(traits_test$y == preds_dep, 100, 0)
  results[[i]][, 2] <- c_dep
  for (j in 1:n_trials) {
    results[[i]][j, 3] <- calc_log_loss(
      true = traits_test$y[j],
      pred = pred_probs_dep[j]
    )
  }
}
for (i in 1:length(scenarios)) {
  results[[i]] |>
    as.data.frame() |>
    readr::write_csv(
      file = paste0("temp_results_", scenarios[i], ".csv")
    )
}
