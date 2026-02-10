# Load packages. ----
library(export)
library(LaplacesDemon)
library(tidyverse)

# Calculate median transition rates, averaged across 1,000 trials. ----
scenarios <- c("ER.L", "ER.M", "ER.H", "DR.LM", "DR.LH", "DR.MH", "DEP1.L",
               "DEP1.M", "DEP1.H", "DEP2.L", "DEP2.M", "DEP2.H", "random")
n_trials <- 1000
n_samps <- 1001
# Dependent.
q_median_dep <- setNames(
  object = vector(mode = "list", length = length(scenarios)),
  nm = scenarios
)
for (i in 1:length(scenarios)) {
  q_median_dep[[i]] <- vector(mode = "list", length = n_trials)
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
    temp2 <- matrix(
      data = c(
        0, median(temp1$q12), median(temp1$q13), 0,
        median(temp1$q21), 0, 0, median(temp1$q24),
        median(temp1$q31), 0, 0, median(temp1$q34),
        0, median(temp1$q42), median(temp1$q43), 0
      ),
      nrow = 4,
      ncol = 4,
      byrow = TRUE
    )
    diag(temp2) <- -rowSums(temp2)
    q_median_dep[[i]][[j]] <- temp2
    rm(temp1)
    rm(temp2)
  }
  setwd("..")
  print(paste0("Finished ", scenarios[i]))
}
setwd("..")
save(q_median_dep, file = "qmat_median_2_dep.RData")
check_dep <- setNames(
  object = vector(mode = "list", length = length(scenarios)),
  nm = scenarios
)
for (i in 1:length(scenarios)) {
  temp <- simplify2array(q_median_dep[[i]])
  temp <- apply(X = temp, MARGIN = 1:2, FUN = median)
  check_dep[[i]] <- temp
  rm(temp)
}
print(check_dep, digits = 2)
#> $ER.L
#>      [,1] [,2] [,3] [,4]
#> [1,] -2.5  1.3  1.2  0.0
#> [2,]  1.6 -2.9  0.0  1.3
#> [3,]  1.5  0.0 -2.9  1.3
#> [4,]  0.0  1.7  1.7 -3.4
#>
#> $ER.M
#>      [,1] [,2] [,3] [,4]
#> [1,] -3.3  1.7  1.7  0.0
#> [2,]  1.7 -3.3  0.0  1.7
#> [3,]  1.6  0.0 -3.3  1.7
#> [4,]  0.0  1.7  1.7 -3.3
#>
#> $ER.H
#>      [,1] [,2] [,3] [,4]
#> [1,] -3.5  1.7  1.7  0.0
#> [2,]  1.7 -3.5  0.0  1.7
#> [3,]  1.7  0.0 -3.5  1.8
#> [4,]  0.0  1.7  1.7 -3.5
#>
#> $DR.LM
#>      [,1]  [,2]  [,3] [,4]
#> [1,] -4.9  2.48  2.42  0.0
#> [2,]  1.0 -3.28  0.00  2.2
#> [3,]  1.0  0.00 -3.26  2.2
#> [4,]  0.0  0.82  0.82 -1.7
#>
#> $DR.LH
#>       [,1]  [,2]  [,3]  [,4]
#> [1,] -5.76  2.93  2.90  0.00
#> [2,]  0.84 -3.75  0.00  2.95
#> [3,]  0.85  0.00 -3.78  2.92
#> [4,]  0.00  0.48  0.47 -0.98
#>
#> $DR.MH
#>      [,1] [,2] [,3] [,4]
#> [1,] -4.8  2.4  2.4  0.0
#> [2,]  1.2 -3.5  0.0  2.2
#> [3,]  1.2  0.0 -3.5  2.3
#> [4,]  0.0  1.0  1.1 -2.1
#>
#> $DEP1.L
#>      [,1] [,2]  [,3] [,4]
#> [1,] -1.8  1.1  0.67  0.0
#> [2,]  1.3 -2.5  0.00  1.1
#> [3,]  2.3  0.0 -4.24  2.0
#> [4,]  0.0  1.7  1.04 -2.7
#>
#> $DEP1.M
#>      [,1] [,2]  [,3] [,4]
#> [1,] -2.6  1.6  0.97  0.0
#> [2,]  1.6 -3.2  0.00  1.6
#> [3,]  2.5  0.0 -4.90  2.4
#> [4,]  0.0  1.6  1.01 -2.6
#>
#> $DEP1.H
#>      [,1] [,2] [,3] [,4]
#> [1,] -2.9  1.7  1.1  0.0
#> [2,]  1.7 -3.4  0.0  1.7
#> [3,]  2.5  0.0 -5.1  2.6
#> [4,]  0.0  1.7  1.1 -2.7
#>
#> $DEP2.L
#>      [,1]  [,2]  [,3]  [,4]
#> [1,] -1.1  0.88  0.18  0.00
#> [2,]  1.2 -2.05  0.00  0.94
#> [3,]  3.0  0.00 -5.76  2.72
#> [4,]  0.0  1.72  0.38 -2.19
#>
#> $DEP2.M
#>      [,1] [,2]  [,3] [,4]
#> [1,] -1.9  1.5  0.36  0.0
#> [2,]  1.5 -2.8  0.00  1.4
#> [3,]  3.3  0.0 -6.50  3.2
#> [4,]  0.0  1.5  0.38 -1.9
#>
#> $DEP2.H
#>      [,1] [,2]  [,3] [,4]
#> [1,] -2.1  1.6  0.42  0.0
#> [2,]  1.6 -3.1  0.00  1.6
#> [3,]  3.4  0.0 -6.79  3.4
#> [4,]  0.0  1.6  0.43 -2.1
#>
#> $random
#>      [,1] [,2] [,3] [,4]
#> [1,] -4.2  2.1  2.1  0.0
#> [2,]  2.1 -4.3  0.0  2.1
#> [3,]  2.1  0.0 -4.2  2.1
#> [4,]  0.0  2.1  2.1 -4.2
