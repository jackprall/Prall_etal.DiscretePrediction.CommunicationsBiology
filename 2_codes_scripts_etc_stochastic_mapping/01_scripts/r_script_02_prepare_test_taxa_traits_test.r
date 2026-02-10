# Prepare lists of test taxa and test datasets. ----
# ER.L - DEP2.H.
scenarios <- c("DEP1.H", "DEP1.L", "DEP1.M", "DEP2.H", "DEP2.L", "DEP2.M",
               "DR.LH", "DR.LM", "DR.MH", "ER.H", "ER.L", "ER.M")
for (i in scenarios) {
  setwd("Single")
  temp <- read.table(
    file = paste0(i, ".Single.ResultsFull.txt"),
    header = TRUE,
    sep = "\t"
  )
  setwd("..")
  test_taxa <- temp[, 2]
  save(test_taxa, file = paste0("test_taxa_", i, ".RData"))
  temp <- temp[, 2:4]
  colnames(temp) <- c("taxon", "x", "y")
  temp$taxon <- as.character(temp$taxon)
  temp$xy <- as.factor(paste0(temp$x, temp$y))
  traits_test <- split(x = temp, f = 1:nrow(temp))
  save(traits_test, file = paste0("traits_test_", i, ".RData"))
}
# Random.
temp <- read.table(
  file = "Random.Single.ResultsFull.txt",
  header = TRUE,
  sep = "\t"
)
test_taxa <- temp[, 2]
save(test_taxa, file = paste0("test_taxa_random.RData"))
temp <- temp[, 2:4]
colnames(temp) <- c("taxon", "x", "y")
temp$taxon <- as.character(temp$taxon)
temp$xy <- as.factor(paste0(temp$x, temp$y))
traits_test <- split(x = temp, f = 1:nrow(temp))
save(traits_test, file = paste0("traits_test_random.RData"))
