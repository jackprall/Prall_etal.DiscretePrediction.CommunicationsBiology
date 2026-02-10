# Load packages. ----
library(export)
library(LaplacesDemon)
library(tidyverse)

# Check MCMC diagnostics. ----
scenarios <- c("ER.L", "ER.M", "ER.H", "DR.LM", "DR.LH", "DR.MH", "DEP1.L",
               "DEP1.M", "DEP1.H", "DEP2.L", "DEP2.M", "DEP2.H", "random")
n_trials <- 1000
n_samps <- 1001
# Dependent.
ess_dep <- setNames(
  object = vector(mode = "list", length = length(scenarios)),
  nm = scenarios
)
for (i in 1:length(scenarios)) {
  ess_dep[[i]] <- vector(mode = "list", length = n_trials)
}
setwd("posterior_bayestraits")
for (i in 1:length(scenarios)) {
  setwd(scenarios[i])
  for (j in 1:n_trials) {
    temp1 <- readr::read_tsv(
      file = paste0(scenarios[i], ".Dep.MCMC.", j, ".Rates.Log.txt"),
      skip = 53,
      name_repair = "minimal",
      show_col_types = FALSE
    )
    ess_dep[[i]][j] <- LaplacesDemon::ESS(x = temp1$Lh)
    rm(temp1)
  }
  setwd("..")
  print(paste0("Finished ", scenarios[i]))
}
setwd("..")
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
#>  1 DEP1.H        1001
#>  2 DEP1.L         881.
#>  3 DEP1.M        1001
#>  4 DEP2.H        1001
#>  5 DEP2.L         958.
#>  6 DEP2.M        1001
#>  7 DR.LH         1001
#>  8 DR.LM         1001
#>  9 DR.MH         1001
#> 10 ER.H          1001
#> 11 ER.L           882.
#> 12 ER.M          1001
#> 13 random        1001
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
      subtitle = "BayesTraits: Dependent",
      x = "Effective Sample Size",
      y = "Frequency"
    )
export::graph2pdf(
  file = "plot_mcmc_diag_2_dep.pdf",
  width = 6,
  height = 8
)
