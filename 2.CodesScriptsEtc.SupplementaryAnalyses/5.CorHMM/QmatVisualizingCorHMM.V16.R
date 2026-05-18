###


### This script was written to visualize the results of 

### First, get simulation bash file and any useful functions for this script
bash_file <- Sys.glob("Discrete_Simulation.V*")
version_str <- sub(".*(V[0-9]+).*", "\\1", bash_file)
source(paste0("Scripts/DiscreteFunctions.", version_str, ".R"))

# Get the types and trials arrays
types <- grep("^types=", readLines(bash_file), value = TRUE)
types <- gsub("types=\\(|\\)", "", types)
types <- gsub("\"", "", types)
types <- strsplit(types, " ")[[1]]

# Get the parameters we need to make the trees
pop_size <- grep("^pop_size=", readLines(bash_file), value = TRUE)
pop_size <- as.numeric(gsub("[^0-9.-]", "", pop_size))

# Get the special cases you may run
variable_rates <- read_logical("variable_rates", bash_file)

# Number of iterations when running as batch job/ on HPC
trial_amount <- grep("^num_iterations=", readLines(bash_file), value = TRUE)
trial_amount <- as.numeric(gsub("[^0-9.-]", "", trial_amount))



### Install and call the necessary 
# First, a list of the packages we need
required_packages <- c("ape", "phytools", "corHMM")

# Install any packages that aren't already installed
installed <- required_packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(required_packages[!installed])
}

# Load the libraries
lapply(required_packages, library, character.only = TRUE)




### First display the accuracies and LL scores



### Now, compare the estimated matrices

# Establish the difference between corHMM, MCMC, and RJ tests
labels <- c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43")
labels1 <- paste0("corHMM_", labels)
cor_labels <- paste0("Est_", labels)
labels2 <- paste0("MCMC_", labels)
labels3 <- paste0("RJ_", labels)
summary_row_labels <- c(labels1, labels2, labels3, "MCMC_ESS_Dep", "RJ_ESS_Dep")
BT_labels <- c(labels2, labels3, "MCMC_ESS_Dep", "RJ_ESS_Dep")

# Gather the column names
summary_col_labels <- c(types, "Random")

# Establish the summary matrix
summary <- matrix(data = NA, nrow = length(summary_row_labels), ncol = length(summary_col_labels), dimnames = list(summary_row_labels, summary_col_labels))



### Loop through each matrix type, and pull the relevant information
for (type in summary_col_labels) {
  # Align the types and their data
  if (type == "Random") {
    cor_data <- read.delim("Results/Random/corHMM/Random.corHmmResults.txt")
    BT_data <- read.delim("Results/Random/QmatInfo/Random.QmatInfo.txt")
  } else {
    cor_name <- paste0("Results/ConstantRates/Single/corHMM/", type, ".corHmmResults.txt")
    data_name <- paste0("Results/ConstantRates/Single/QmatInfo/", type, ".QmatInfo.txt")
    cor_data <- read.delim(cor_name)
    BT_data <- read.delim(data_name)
  }
  
  # Establish the number of trials to summarize
  trial_amount <- min(nrow(cor_data), nrow(BT_data))
  
  # Pull the transition rate info, and store it in the summary table
  for (i in 1:length(labels)) {
    summary[labels1[i], type] <- median(cor_data[, cor_labels[i]])
    summary[labels2[i], type] <- median(BT_data[1:trial_amount, labels2[i]])
    summary[labels3[i], type] <- median(BT_data[1:trial_amount, labels3[i]])
  }
  
  # Then pull the ESS numbers for BayesTraits
  summary["MCMC_ESS_Dep", type] <- median(BT_data[1:trial_amount, "MCMC_ESS_Dep"])
  summary["RJ_ESS_Dep", type] <- median(BT_data[1:trial_amount, "RJ_ESS_Dep"])
}


# Write the summary table
SummaryName <- paste0("Results/ConstantRates/Single/corHMM/Summary.corHMM.QmatEstimates.txt")
write.table(summary, file = SummaryName, quote = F, row.names = F, col.names = T, sep = "\t")
