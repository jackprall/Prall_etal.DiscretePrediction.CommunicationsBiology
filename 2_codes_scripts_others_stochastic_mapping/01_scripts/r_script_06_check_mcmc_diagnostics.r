# Load packages. ----
library(export)
library(miceadds)
library(LaplacesDemon)
library(tidyverse)

# Check MCMC diagnostics. ----
scenarios <- c("ER.L", "ER.M", "ER.H", "DR.LM", "DR.LH", "DR.MH", "DEP1.L",
               "DEP1.M", "DEP1.H", "DEP2.L", "DEP2.M", "DEP2.H", "random")
n_trials <- 1000
n_samps <- 400
# 1. Independent.
loglh_ind <- ess_ind <- setNames(
  object = vector(mode = "list", length = length(scenarios)),
  nm = scenarios
)
for (i in 1:length(scenarios)) {
  loglh_ind[[i]] <- vector(mode = "list", length = n_trials)
  ess_ind[[i]] <- vector(mode = "list", length = n_trials)
}
for (i in 1:length(scenarios)) {
  setwd("outputs")
  miceadds::load.Rdata(
    filename = paste0("simmap_posteriors_", scenarios[i], "_1_ind.RData"),
    objname = "temp2"
  )
  setwd("..")
  print("Finished loading posterior samples...")
  for (j in 1:n_trials) {
    temp1 <- vector(length = n_samps)
    for (k in 1:n_samps) {
      temp1[k] <- as.numeric(temp2[[j]][[k]]$logL)
    }
    loglh_ind[[i]][[j]] <- temp1
    ess_ind[[i]][j] <- LaplacesDemon::ESS(x = temp1)
    rm(temp1)
  }
  rm(temp2)
  print(paste0("Finished ", scenarios[i]))
}
save(loglh_ind, file = "mcmc_diag_loglh_1_ind.RData")
save(ess_ind, file = "mcmc_diag_ess_1_ind.RData")
dat_ind <- data.frame(
  scenario = rep(scenarios, each = 1000),
  ess = unlist(ess_ind)
)
summary_ind <- dat_ind |>
  dplyr::group_by(scenario) |>
  dplyr::summarize(
    median_ess = median(ess)
  )
summary_ind
#> # A tibble: 13 × 2
#>    scenario median_ess
#>    <chr>         <dbl>
#>  1 DEP1.H        270.
#>  2 DEP1.L         31.3
#>  3 DEP1.M        147.
#>  4 DEP2.H        273.
#>  5 DEP2.L         52.6
#>  6 DEP2.M        154.
#>  7 DR.LH         200.
#>  8 DR.LM          95.6
#>  9 DR.MH         228.
#> 10 ER.H          272.
#> 11 ER.L           30.6
#> 12 ER.M          154.
#> 13 random         32.1
dat_ind |>
  dplyr::mutate(trial = rep(1:n_trials, length(scenarios))) |>
  tidyr::pivot_wider(
    names_from = scenario,
    values_from = ess
  ) |>
  readr::write_csv(file = "temp_ess_1_ind.csv")
dat_ind |>
  dplyr::mutate(
    scenario = forcats::fct_relevel(
      scenario,
      c("ER.L", "ER.M", "ER.H", "DR.LM", "DR.LH", "DR.MH", "DEP1.L",
        "DEP1.M", "DEP1.H", "DEP2.L", "DEP2.M", "DEP2.H", "random")
    )
  ) |>
  ggplot2::ggplot(mapping = aes(x = ess, color = scenario)) +
    facet_wrap(facets = vars(scenario), ncol = 3) +
    geom_histogram(color = "white", fill = "dark gray") +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(
      title = "MCMC Diagnostics",
      subtitle = "SIMMAP: Independent",
      x = "Effective Sample Size",
      y = "Frequency"
    )
export::graph2pdf(
  file = "plot_mcmc_diag_1_ind.pdf",
  width = 6,
  height = 8
)
# 2. Dependent.
loglh_dep <- ess_dep <- setNames(
  object = vector(mode = "list", length = length(scenarios)),
  nm = scenarios
)
for (i in 1:length(scenarios)) {
  loglh_dep[[i]] <- vector(mode = "list", length = n_trials)
  ess_dep[[i]] <- vector(mode = "list", length = n_trials)
}
for (i in 1:length(scenarios)) {
  setwd("outputs")
  miceadds::load.Rdata(
    filename = paste0("simmap_posteriors_", scenarios[i], "_2_dep.RData"),
    objname = "temp2"
  )
  setwd("..")
  print("Finished loading posterior samples...")
  for (j in 1:n_trials) {
    temp1 <- vector(length = n_samps)
    for (k in 1:n_samps) {
      temp1[k] <- as.numeric(temp2[[j]][[k]]$logL)
    }
    loglh_dep[[i]][[j]] <- temp1
    ess_dep[[i]][j] <- LaplacesDemon::ESS(x = temp1)
    rm(temp1)
  }
  rm(temp2)
  print(paste0("Finished ", scenarios[i]))
}
save(loglh_dep, file = "mcmc_diag_loglh_2_dep.RData")
save(ess_dep, file = "mcmc_diag_ess_2_dep.RData")
dat_dep <- data.frame(
  scenario = rep(scenarios, each = 1000),
  ess = unlist(ess_dep)
)
summary_dep <- dat_dep |>
  dplyr::group_by(scenario) |>
  dplyr::summarize(
    median_ess = median(ess)
  )
summary_dep
#> # A tibble: 13 × 2
#>    scenario median_ess
#>    <chr>         <dbl>
#>  1 DEP1.H        312.
#>  2 DEP1.L         42.3
#>  3 DEP1.M        196.
#>  4 DEP2.H        295.
#>  5 DEP2.L         47.4
#>  6 DEP2.M        189.
#>  7 DR.LH         230.
#>  8 DR.LM         123.
#>  9 DR.MH         271.
#> 10 ER.H          306.
#> 11 ER.L           45.6
#> 12 ER.M          195.
#> 13 random         42.7
dat_dep |>
  dplyr::mutate(trial = rep(1:n_trials, length(scenarios))) |>
  tidyr::pivot_wider(
    names_from = scenario,
    values_from = ess
  ) |>
  readr::write_csv(file = "temp_ess_2_dep.csv")
dat_dep |>
  dplyr::mutate(
    scenario = forcats::fct_relevel(
      scenario,
      c("ER.L", "ER.M", "ER.H", "DR.LM", "DR.LH", "DR.MH", "DEP1.L",
        "DEP1.M", "DEP1.H", "DEP2.L", "DEP2.M", "DEP2.H", "random")
    )
  ) |>
  ggplot2::ggplot(mapping = aes(x = ess, color = scenario)) +
    facet_wrap(facets = vars(scenario), ncol = 3) +
    geom_histogram(color = "white", fill = "dark gray") +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(
      title = "MCMC Diagnostics",
      subtitle = "SIMMAP: Dependent",
      x = "Effective Sample Size",
      y = "Frequency"
    )
export::graph2pdf(
  file = "plot_mcmc_diag_2_dep.pdf",
  width = 6,
  height = 8
)
