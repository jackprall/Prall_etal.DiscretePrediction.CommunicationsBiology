# Prepare training datasets. ----
# ER.L - DEP2.H.
scenarios <- c("DEP1.H", "DEP1.L", "DEP1.M", "DEP2.H", "DEP2.L", "DEP2.M",
               "DR.LH", "DR.LM", "DR.MH", "ER.H", "ER.L", "ER.M")
n_trials <- 1000
colnames_temp <- c("taxon", "x", "y", "xy")
coltypes_temp <- c("character", "numeric", "numeric", "factor")
for (i in scenarios) {
  load(file = paste0("test_taxa_", i, ".RData"))
  setwd(paste0("ConstantRates/", i))
  traits_train <- vector(mode = "list", length = n_trials)
  for (j in 1:n_trials) {
    temp <- read.table(
      file = paste0(i, ".", j, ".Full_data.txt"),
      header = FALSE,
      sep = "\t",
      skip = 1,
      col.names = colnames_temp,
      colClasses = coltypes_temp
    )
    temp <- temp[!temp$taxon %in% test_taxa[j], ]
    traits_train[[j]] <- temp
  }
  setwd("../..")
  save(traits_train, file = paste0("traits_train_", i,".RData"))
}
# Random.
load(file = "test_taxa_random.RData")
setwd("Random")
traits_train <- vector(mode = "list", length = n_trials)
for (j in 1:n_trials) {
  temp <- read.table(
    file = paste0("Random.", j, ".Full_data.txt"),
    header = FALSE,
    sep = "\t",
    skip = 1,
    col.names = colnames_temp,
    colClasses = coltypes_temp
  )
  temp <- temp[!temp$taxon %in% test_taxa[j], ]
  traits_train[[j]] <- temp
}
setwd("..")
save(traits_train, file = "traits_train_random.RData")
