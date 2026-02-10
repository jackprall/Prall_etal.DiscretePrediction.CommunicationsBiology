# Load packages. ----
library(export)
library(miceadds)
library(LaplacesDemon)
library(tidyverse)

# Calculate median transition rates, averaged across 1,000 trials. ----
scenarios <- c("ER.L", "ER.M", "ER.H", "DR.LM", "DR.LH", "DR.MH", "DEP1.L",
               "DEP1.M", "DEP1.H", "DEP2.L", "DEP2.M", "DEP2.H", "random")
n_trials <- 1000
n_samps <- 400
# 1. Independent.
q_median_ind <- setNames(
  object = vector(mode = "list", length = length(scenarios)),
  nm = scenarios
)
for (i in 1:length(scenarios)) {
  q_median_ind[[i]] <- vector(mode = "list", length = n_trials)
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
    temp1 <- vector(mode = "list", length = n_samps)
    for (k in 1:n_samps) {
      temp1[[k]] <- temp2[[j]][[k]]$Q
    }
    temp1 <- simplify2array(temp1)
    temp1 <- apply(X = temp1, MARGIN = 1:2, FUN = median)
    q_median_ind[[i]][[j]] <- temp1
    rm(temp1)
  }
  rm(temp2)
  print(paste0("Finished ", scenarios[i]))
}
save(q_median_ind, file = "qmat_median_1_ind.RData")
check_ind <- setNames(
  object = vector(mode = "list", length = length(scenarios)),
  nm = scenarios
)
for (i in 1:length(scenarios)) {
  temp <- simplify2array(q_median_ind[[i]])
  temp <- apply(X = temp, MARGIN = 1:2, FUN = median)
  check_ind[[i]] <- temp
  rm(temp)
}
print(check_ind, digits = 2)  # on-diagonals have not been adjusted
#> $ER.L
#>        00     01     10     11
#> 00 -0.273  0.070  0.065  0.000
#> 01  0.089 -0.300  0.000  0.065
#> 10  0.077  0.000 -0.300  0.070
#> 11  0.000  0.077  0.089 -0.330
#>
#> $ER.M
#>       00    01    10    11
#> 00 -0.88  0.34  0.34  0.00
#> 01  0.34 -0.87  0.00  0.34
#> 10  0.34  0.00 -0.89  0.34
#> 11  0.00  0.34  0.34 -0.86
#>
#> $ER.H
#>       00    01    10    11
#> 00 -1.27  0.54  0.55  0.00
#> 01  0.54 -1.28  0.00  0.55
#> 10  0.55  0.00 -1.27  0.54
#> 11  0.00  0.55  0.54 -1.28
#>
#> $DR.LM
#>       00    01    10    11
#> 00 -1.09  0.35  0.37  0.00
#> 01  0.15 -0.72  0.00  0.37
#> 10  0.15  0.00 -0.70  0.35
#> 11  0.00  0.15  0.15 -0.42
#>
#> $DR.LH
#>       00    01    10    11
#> 00 -1.71  0.62  0.86  0.00
#> 01  0.21 -1.38  0.00  0.86
#> 10  0.19  0.00 -0.92  0.62
#> 11  0.00  0.19  0.21 -0.49
#>
#> $DR.MH
#>       00    01    10    11
#> 00 -1.60  0.68  0.68  0.00
#> 01  0.32 -1.18  0.00  0.68
#> 10  0.31  0.00 -1.17  0.68
#> 11  0.00  0.31  0.32 -0.77
#>
#> $DEP1.L
#>        00     01     10     11
#> 00 -0.244  0.078  0.051  0.000
#> 01  0.075 -0.243  0.000  0.051
#> 10  0.102  0.000 -0.368  0.078
#> 11  0.000  0.102  0.075 -0.339
#>
#> $DEP1.M
#>       00    01    10    11
#> 00 -0.84  0.39  0.27  0.00
#> 01  0.29 -0.73  0.00  0.27
#> 10  0.39  0.00 -1.01  0.39
#> 11  0.00  0.39  0.29 -0.87
#>
#> $DEP1.H
#>       00    01    10    11
#> 00 -1.25  0.64  0.44  0.00
#> 01  0.43 -1.04  0.00  0.44
#> 10  0.64  0.00 -1.49  0.64
#> 11  0.00  0.64  0.43 -1.26
#>
#> $DEP2.L
#>       00    01    10    11
#> 00 -0.31  0.13  0.06  0.00
#> 01  0.18 -0.37  0.00  0.06
#> 10  0.20  0.00 -0.62  0.13
#> 11  0.00  0.20  0.18 -0.76
#>
#> $DEP2.M
#>       00    01    10    11
#> 00 -0.79  0.40  0.21  0.00
#> 01  0.38 -0.76  0.00  0.21
#> 10  0.45  0.00 -1.14  0.40
#> 11  0.00  0.45  0.38 -1.07
#>
#> $DEP2.H
#>       00    01    10    11
#> 00 -1.11  0.60  0.36  0.00
#> 01  0.48 -0.99  0.00  0.36
#> 10  0.71  0.00 -1.56  0.60
#> 11  0.00  0.71  0.48 -1.41
#>
#> $random
#>        00     01     10     11
#> 00 -0.248  0.062  0.070  0.000
#> 01  0.074 -0.275  0.000  0.070
#> 10  0.086  0.000 -0.280  0.062
#> 11  0.000  0.086  0.074 -0.307
# 2. Dependent.
q_median_dep <- setNames(
  object = vector(mode = "list", length = length(scenarios)),
  nm = scenarios
)
for (i in 1:length(scenarios)) {
  q_median_dep[[i]] <- vector(mode = "list", length = n_trials)
}
for (i in 1:length(scenarios)) {
  #setwd("outputs")
  miceadds::load.Rdata(
    filename = paste0("simmap_posteriors_", scenarios[i], "_2_dep.RData"),
    objname = "temp2"
  )
  #setwd("..")
  print("Finished loading posterior samples...")
  for (j in 1:n_trials) {
    temp1 <- vector(mode = "list", length = n_samps)
    for (k in 1:n_samps) {
      temp1[[k]] <- temp2[[j]][[k]]$Q
    }
    temp1 <- simplify2array(temp1)
    temp1 <- apply(X = temp1, MARGIN = 1:2, FUN = median)
    q_median_dep[[i]][[j]] <- temp1
    rm(temp1)
  }
  rm(temp2)
  print(paste0("Finished ", scenarios[i]))
}
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
#>       00    01    10    11
#> 00 -0.83  0.32  0.31  0.00
#> 01  0.40 -0.92  0.00  0.34
#> 10  0.38  0.00 -0.96  0.36
#> 11  0.00  0.45  0.48 -1.19
#>
#> $ER.M
#>       00    01    10    11
#> 00 -1.31  0.56  0.56  0.00
#> 01  0.56 -1.30  0.00  0.55
#> 10  0.56  0.00 -1.31  0.55
#> 11  0.00  0.56  0.56 -1.32
#>
#> $ER.H
#>       00    01    10    11
#> 00 -1.48  0.65  0.65  0.00
#> 01  0.65 -1.48  0.00  0.66
#> 10  0.67  0.00 -1.50  0.65
#> 11  0.00  0.65  0.65 -1.48
#>
#> $DR.LM
#>       00    01    10    11
#> 00 -1.93  0.84  0.82  0.00
#> 01  0.36 -1.21  0.00  0.64
#> 10  0.34  0.00 -1.16  0.62
#> 11  0.00  0.24  0.24 -0.60
#>
#> $DR.LH
#>       00    01    10    11
#> 00 -2.29  0.94  1.05  0.00
#> 01  0.43 -1.43  0.00  0.81
#> 10  0.31  0.00 -1.25  0.76
#> 11  0.00  0.22  0.21 -0.48
#>
#> $DR.MH
#>       00    01    10    11
#> 00 -2.07  0.90  0.91  0.00
#> 01  0.47 -1.42  0.00  0.79
#> 10  0.44  0.00 -1.41  0.79
#> 11  0.00  0.36  0.36 -0.85
#>
#> $DEP1.L
#>       00    01    10    11
#> 00 -0.55  0.25  0.17  0.00
#> 01  0.31 -0.81  0.00  0.28
#> 10  0.67  0.00 -1.53  0.59
#> 11  0.00  0.47  0.32 -1.00
#>
#> $DEP1.M
#>       00    01    10    11
#> 00 -1.01  0.52  0.35  0.00
#> 01  0.53 -1.21  0.00  0.51
#> 10  0.86  0.00 -1.97  0.83
#> 11  0.00  0.53  0.35 -1.04
#>
#> $DEP1.H
#>       00    01    10    11
#> 00 -1.21  0.65  0.43  0.00
#> 01  0.62 -1.40  0.00  0.60
#> 10  0.96  0.00 -2.17  0.95
#> 11  0.00  0.63  0.42 -1.18
#>
#> $DEP2.L
#>       00    01    10    11
#> 00 -0.47  0.24  0.11  0.00
#> 01  0.31 -0.61  0.00  0.12
#> 10  0.74  0.00 -1.39  0.32
#> 11  0.00  1.02  0.63 -2.10
#>
#> $DEP2.M
#>       00    01    10    11
#> 00 -0.87  0.52  0.21  0.00
#> 01  0.51 -0.97  0.00  0.28
#> 10  0.95  0.00 -2.27  0.95
#> 11  0.00  0.73  0.29 -1.13
#>
#> $DEP2.H
#>       00    01    10    11
#> 00 -1.02  0.63  0.24  0.00
#> 01  0.61 -1.17  0.00  0.44
#> 10  1.13  0.00 -2.64  1.12
#> 11  0.00  0.76  0.27 -1.13
#>
#> $random
#>       00    01    10    11
#> 00 -0.80  0.31  0.32  0.00
#> 01  0.37 -0.95  0.00  0.35
#> 10  0.40  0.00 -0.99  0.35
#> 11  0.00  0.43  0.41 -1.07
