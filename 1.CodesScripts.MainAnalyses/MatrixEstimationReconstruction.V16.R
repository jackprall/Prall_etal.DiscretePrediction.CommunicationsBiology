### This Script checks how effective BayesTraits was at recapturing the evolutionary model
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

# Get the MCMC type
RJmodel <- grep("^RJmodel=", readLines(bash_file), value = TRUE)
RJmodel <- gsub(".*=(.*)", "\\1", RJmodel)

# Number of iterations when running as batch job/ on HPC
trial_amount <- grep("^num_iterations=", readLines(bash_file), value = TRUE)
trial_amount <- as.numeric(gsub("[^0-9.-]", "", trial_amount))

# Keep track of the original directory
original_dir <- getwd()





### Install and call the necessary 
# First, a list of the packages we need
required_packages <- c("ape", "phytools", "castor", "LaplacesDemon")

# Install any packages that aren't already installed
installed <- required_packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(required_packages[!installed])
}

# Load the libraries
lapply(required_packages, library, character.only = TRUE)





# Set up summary matrix, first establish the row names
if (RJmodel == "MCMC" || RJmodel == "RJMCMC") {
  summary_full_labels <- c("Alpha_1", "Alpha_2", "Beta_1", "Beta_2", "q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43", "ESS_Ind", "ESS_Dep")
}
if (RJmodel == "BOTH") {
  # Establish the difference between MCMC and RJ tests
  labels <- c("Alpha_1", "Alpha_2", "Beta_1", "Beta_2", "q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43")
  labels1 <- paste0("MCMC_", labels)
  labels2 <- paste0("RJ_", labels)
  summary_row_labels <- c(labels1, labels2, "MCMC_ESS_Ind", "MCMC_ESS_Dep", "RJ_ESS_Ind", "RJ_ESS_Dep")
}
# Gather the column names
summary_col_labels <- c(types, "Random")

# Establish the summary matrix
summary <- matrix(data = NA, nrow = length(summary_row_labels), ncol = length(summary_col_labels), dimnames = list(summary_row_labels, summary_col_labels))



for (type in types) {
  # Print a progress check
  print(paste0("Beginning to summarize estimated rates matrices from ", type, " tests with Constant Rates."))
  
  # Call the necessary paths
  results_name <- paste0("Results/ConstantRates/Single/", type, ".Single.ResultsFull.txt")
  anccon_name <- paste0("Results/ConstantRates/Single/AncestralStateReconstruction/", type, ".AncReconAndTreeHeight.txt")
  
  # Call the results table  
  colnames_results <- strsplit(readLines(results_name, n = 1), "\t")[[1]]
  results <- read.table(results_name, skip = 1, sep = "\t", col.names = colnames_results)
  
  # Call the Ancestral State Reconstruction Table
  colnames_anccon <- strsplit(readLines(anccon_name, n = 1), "\t")[[1]]
  AncCon <- read.table(anccon_name, skip = 1, sep = "\t", col.names = colnames_anccon)
  
  # Create a new matrix for each set of simulated data
  QmatMat <- matrix(data=NA, nrow = nrow(results), ncol = (10+length(summary_row_labels)))
  colnames(QmatMat) <- c("Trial_#", "Unknown_Taxon", "True_4States", "TBL", "Tree_Height", 
                         "Std_TBL", "n00", "n01", "n10", "n11", summary_row_labels)
  
  # Pull the information relevant to the rates matrices and store it in a new matrix
  QmatMat[, "Trial_#"] <- results[, "Trial"]
  QmatMat[, "Unknown_Taxon"] <- results[, "rTaxon"]
  QmatMat[, "True_4States"] <- results[, "True_4States"]
  QmatMat[, "n00"] <- results[, "n00"]
  QmatMat[, "n01"] <- results[, "n01"]
  QmatMat[, "n10"] <- results[, "n10"]
  QmatMat[, "n11"] <- results[, "n11"]
  QmatMat[, "TBL"] <- AncCon[, "TBL"]
  QmatMat[, "Tree_Height"] <- AncCon[, "Tree_Height"]
  QmatMat[, "Std_TBL"] <- AncCon[, "Std_TBL"]
  
  for (i in 1:nrow(results)) {
    if (RJmodel == "MCMC" || RJmodel == "BOTH") {
      # Gather the independent q-matrix
      Qmat_IndMCMC <- estimated_qmat(type = type,
                                     i = i,
                                     IndependentRates = TRUE,
                                     ConstantRates = TRUE,
                                     RJ = FALSE,
                                     Predict = FALSE)
      
      # Record the rates for this trial
      QmatMat[i, "MCMC_Alpha_1"] <- Qmat_IndMCMC[1, 3]
      QmatMat[i, "MCMC_Alpha_2"] <- Qmat_IndMCMC[1, 2] 
      QmatMat[i, "MCMC_Beta_1"] <- Qmat_IndMCMC[3, 1]
      QmatMat[i, "MCMC_Beta_2"] <- Qmat_IndMCMC[2, 1]
      
      # Pull the log file and record the ESS
      logname <- paste0("ConstantRates/", type, "/Single/", type, ".Ind.MCMC.", i, ".Rates.Log.txt")
      logdata <- read.delim(logname, skip = 45)
      QmatMat[i, "MCMC_ESS_Ind"] <- ESS(logdata[, "Lh"])
      
      # Gather the dependent q-matrix
      Qmat_DepMCMC <- estimated_qmat(type = type,
                                     i = i,
                                     IndependentRates = FALSE,
                                     ConstantRates = TRUE,
                                     RJ = FALSE,
                                     Predict = FALSE)
      
      # Record the rates for this trial
      QmatMat[i, "MCMC_q12"] <- Qmat_DepMCMC[1, 2]
      QmatMat[i, "MCMC_q13"] <- Qmat_DepMCMC[1, 3] 
      QmatMat[i, "MCMC_q21"] <- Qmat_DepMCMC[2, 1]
      QmatMat[i, "MCMC_q24"] <- Qmat_DepMCMC[2, 4]
      QmatMat[i, "MCMC_q31"] <- Qmat_DepMCMC[3, 1]
      QmatMat[i, "MCMC_q34"] <- Qmat_DepMCMC[3, 4] 
      QmatMat[i, "MCMC_q42"] <- Qmat_DepMCMC[4, 2]
      QmatMat[i, "MCMC_q43"] <- Qmat_DepMCMC[4, 3]
      
      # Pull the log file and record the ESS
      logname <- paste0("ConstantRates/", type, "/Single/", type, ".Dep.MCMC.", i, ".Rates.Log.txt")
      logdata <- read.delim(logname, skip = 53)
      QmatMat[i, "MCMC_ESS_Dep"] <- ESS(logdata[, "Lh"])
    }
    
    if (RJmodel == "RJMCMC" || RJmodel == "BOTH") {
      # Gather the independent q-matrix
      Qmat_IndRJ <- estimated_qmat(type = type,
                                   i = i,
                                   IndependentRates = TRUE,
                                   ConstantRates = TRUE,
                                   RJ = TRUE,
                                   Predict = FALSE)
      
      # Record the rates for this trial
      QmatMat[i, "RJ_Alpha_1"] <- Qmat_IndRJ[1, 3]
      QmatMat[i, "RJ_Alpha_2"] <- Qmat_IndRJ[1, 2] 
      QmatMat[i, "RJ_Beta_1"] <- Qmat_IndRJ[3, 1]
      QmatMat[i, "RJ_Beta_2"] <- Qmat_IndRJ[2, 1]
      
      # Pull the log file and record the ESS
      logname <- paste0("ConstantRates/", type, "/Single/", type, ".Ind.RJMCMC.", i, ".Rates.Log.txt")
      logdata <- read.delim(logname, skip = 44)
      QmatMat[i, "RJ_ESS_Ind"] <- ESS(logdata[, "Lh"])
      
      # Gather the dependent q-matrix
      Qmat_DepRJ <- estimated_qmat(type = type,
                                   i = i,
                                   IndependentRates = FALSE,
                                   ConstantRates = TRUE,
                                   RJ = TRUE,
                                   Predict = FALSE)
      
      # Record the rates for this trial
      QmatMat[i, "RJ_q12"] <- Qmat_DepRJ[1, 2]
      QmatMat[i, "RJ_q13"] <- Qmat_DepRJ[1, 3] 
      QmatMat[i, "RJ_q21"] <- Qmat_DepRJ[2, 1]
      QmatMat[i, "RJ_q24"] <- Qmat_DepRJ[2, 4]
      QmatMat[i, "RJ_q31"] <- Qmat_DepRJ[3, 1]
      QmatMat[i, "RJ_q34"] <- Qmat_DepRJ[3, 4] 
      QmatMat[i, "RJ_q42"] <- Qmat_DepRJ[4, 2]
      QmatMat[i, "RJ_q43"] <- Qmat_DepRJ[4, 3]
      
      # Pull the log file and record the ESS
      logname <- paste0("ConstantRates/", type, "/Single/", type, ".Dep.RJMCMC.", i, ".Rates.Log.txt")
      logdata <- read.delim(logname, skip = 48)
      QmatMat[i, "RJ_ESS_Dep"] <- ESS(logdata[, "Lh"])
    }
  }
  
  # Now, summarize the results
  if (RJmodel == "MCMC" || RJmodel == "BOTH") {
    # Save the median independent rates
    summary["MCMC_Alpha_1", type] <- median(QmatMat[, "MCMC_Alpha_1"])
    summary["MCMC_Alpha_2", type] <- median(QmatMat[, "MCMC_Alpha_2"])
    summary["MCMC_Beta_1", type] <- median(QmatMat[, "MCMC_Beta_1"])
    summary["MCMC_Beta_2", type] <- median(QmatMat[, "MCMC_Beta_2"])
    # Save the median dependent rates
    summary["MCMC_q12", type] <- median(QmatMat[, "MCMC_q12"])
    summary["MCMC_q13", type] <- median(QmatMat[, "MCMC_q13"])
    summary["MCMC_q21", type] <- median(QmatMat[, "MCMC_q21"])
    summary["MCMC_q24", type] <- median(QmatMat[, "MCMC_q24"])
    summary["MCMC_q31", type] <- median(QmatMat[, "MCMC_q31"])
    summary["MCMC_q34", type] <- median(QmatMat[, "MCMC_q34"])
    summary["MCMC_q42", type] <- median(QmatMat[, "MCMC_q42"])
    summary["MCMC_q43", type] <- median(QmatMat[, "MCMC_q43"])
    # Save the median ESS's
    summary["MCMC_ESS_Ind", type] <- median(QmatMat[, "MCMC_ESS_Ind"])
    summary["MCMC_ESS_Dep", type] <- median(QmatMat[, "MCMC_ESS_Dep"])
  }
  
  if (RJmodel == "RJMCMC" || RJmodel == "BOTH") {
    # Save the median independent rates
    summary["RJ_Alpha_1", type] <- median(QmatMat[, "RJ_Alpha_1"])
    summary["RJ_Alpha_2", type] <- median(QmatMat[, "RJ_Alpha_2"])
    summary["RJ_Beta_1", type] <- median(QmatMat[, "RJ_Beta_1"])
    summary["RJ_Beta_2", type] <- median(QmatMat[, "RJ_Beta_2"])
    # Save the median dependent rates
    summary["RJ_q12", type] <- median(QmatMat[, "RJ_q12"])
    summary["RJ_q13", type] <- median(QmatMat[, "RJ_q13"])
    summary["RJ_q21", type] <- median(QmatMat[, "RJ_q21"])
    summary["RJ_q24", type] <- median(QmatMat[, "RJ_q24"])
    summary["RJ_q31", type] <- median(QmatMat[, "RJ_q31"])
    summary["RJ_q34", type] <- median(QmatMat[, "RJ_q34"])
    summary["RJ_q42", type] <- median(QmatMat[, "RJ_q42"])
    summary["RJ_q43", type] <- median(QmatMat[, "RJ_q43"])
    # Save the median ESS's
    summary["RJ_ESS_Ind", type] <- median(QmatMat[, "RJ_ESS_Ind"])
    summary["RJ_ESS_Dep", type] <- median(QmatMat[, "RJ_ESS_Dep"])
  }
  
  # Write the Qmat matrix
  QmatName <- paste0("Results/ConstantRates/Single/QmatInfo/", type, ".QmatInfo.txt")
  write.table(QmatMat, file = QmatName, quote = F, row.names = F, col.names = T, sep = "\t")
}


# Now, pull the random info
# Print a progress check
print("Beginning to summarize estimated rates matrices from tests with Random Rates.")

# Call the necessary paths
results_name <- "Results/Random/Random.Single.ResultsFull.txt"
anccon_name <- "Results/Random/AncestralStateReconstruction/Random.AncReconAndTreeHeight.txt"

# Call the results table  
colnames_results <- strsplit(readLines(results_name, n = 1), "\t")[[1]]
results <- read.table(results_name, skip = 1, sep = "\t", col.names = colnames_results)

# Call the Ancestral State Reconstruction Table
colnames_anccon <- strsplit(readLines(anccon_name, n = 1), "\t")[[1]]
AncCon <- read.table(anccon_name, skip = 1, sep = "\t", col.names = colnames_anccon)

for (i in 1:nrow(results)) {
  if (RJmodel == "MCMC" || RJmodel == "BOTH") {
    # Gather the independent q-matrix
    Qmat_IndMCMC <- estimated_qmat(type = type,
                                   i = i,
                                   IndependentRates = TRUE,
                                   ConstantRates = TRUE,
                                   RJ = FALSE,
                                   Predict = FALSE,
                                   Random = TRUE)
    
    # Record the rates for this trial
    QmatMat[i, "MCMC_Alpha_1"] <- Qmat_IndMCMC[1, 3]
    QmatMat[i, "MCMC_Alpha_2"] <- Qmat_IndMCMC[1, 2] 
    QmatMat[i, "MCMC_Beta_1"] <- Qmat_IndMCMC[3, 1]
    QmatMat[i, "MCMC_Beta_2"] <- Qmat_IndMCMC[2, 1]
    
    # Pull the log file and record the ESS
    logname <- paste0("Random/Single/Random.Ind.MCMC.", i, ".Rates.Log.txt")
    logdata <- read.delim(logname, skip = 45)
    QmatMat[i, "MCMC_ESS_Ind"] <- ESS(logdata[, "Lh"])
    
    # Gather the dependent q-matrix
    Qmat_DepMCMC <- estimated_qmat(type = type,
                                   i = i,
                                   IndependentRates = FALSE,
                                   ConstantRates = TRUE,
                                   RJ = FALSE,
                                   Predict = FALSE,
                                   Random = TRUE)
    
    # Record the rates for this trial
    QmatMat[i, "MCMC_q12"] <- Qmat_DepMCMC[1, 2]
    QmatMat[i, "MCMC_q13"] <- Qmat_DepMCMC[1, 3] 
    QmatMat[i, "MCMC_q21"] <- Qmat_DepMCMC[2, 1]
    QmatMat[i, "MCMC_q24"] <- Qmat_DepMCMC[2, 4]
    QmatMat[i, "MCMC_q31"] <- Qmat_DepMCMC[3, 1]
    QmatMat[i, "MCMC_q34"] <- Qmat_DepMCMC[3, 4] 
    QmatMat[i, "MCMC_q42"] <- Qmat_DepMCMC[4, 2]
    QmatMat[i, "MCMC_q43"] <- Qmat_DepMCMC[4, 3]
    
    # Pull the log file and record the ESS
    logname <- paste0("Random/Single/Random.Dep.MCMC.", i, ".Rates.Log.txt")
    logdata <- read.delim(logname, skip = 53)
    QmatMat[i, "MCMC_ESS_Dep"] <- ESS(logdata[, "Lh"])
  }
  
  if (RJmodel == "RJMCMC" || RJmodel == "BOTH") {
    # Gather the independent q-matrix
    Qmat_IndRJ <- estimated_qmat(type = type,
                                 i = i,
                                 IndependentRates = TRUE,
                                 ConstantRates = TRUE,
                                 RJ = TRUE,
                                 Predict = FALSE,
                                 Random = TRUE)
    
    # Record the rates for this trial
    QmatMat[i, "RJ_Alpha_1"] <- Qmat_IndRJ[1, 3]
    QmatMat[i, "RJ_Alpha_2"] <- Qmat_IndRJ[1, 2] 
    QmatMat[i, "RJ_Beta_1"] <- Qmat_IndRJ[3, 1]
    QmatMat[i, "RJ_Beta_2"] <- Qmat_IndRJ[2, 1]
    
    # Pull the log file and record the ESS
    logname <- paste0("Random/Single/Random.Ind.RJMCMC.", i, ".Rates.Log.txt")
    logdata <- read.delim(logname, skip = 44)
    QmatMat[i, "RJ_ESS_Ind"] <- ESS(logdata[, "Lh"])
    
    # Gather the dependent q-matrix
    Qmat_DepRJ <- estimated_qmat(type = type,
                                 i = i,
                                 IndependentRates = FALSE,
                                 ConstantRates = TRUE,
                                 RJ = TRUE,
                                 Predict = FALSE,
                                 Random = TRUE)
    
    # Record the rates for this trial
    QmatMat[i, "RJ_q12"] <- Qmat_DepRJ[1, 2]
    QmatMat[i, "RJ_q13"] <- Qmat_DepRJ[1, 3] 
    QmatMat[i, "RJ_q21"] <- Qmat_DepRJ[2, 1]
    QmatMat[i, "RJ_q24"] <- Qmat_DepRJ[2, 4]
    QmatMat[i, "RJ_q31"] <- Qmat_DepRJ[3, 1]
    QmatMat[i, "RJ_q34"] <- Qmat_DepRJ[3, 4] 
    QmatMat[i, "RJ_q42"] <- Qmat_DepRJ[4, 2]
    QmatMat[i, "RJ_q43"] <- Qmat_DepRJ[4, 3]
    
    # Pull the log file and record the ESS
    logname <- paste0("Random/Single/Random.Dep.RJMCMC.", i, ".Rates.Log.txt")
    logdata <- read.delim(logname, skip = 48)
    QmatMat[i, "RJ_ESS_Dep"] <- ESS(logdata[, "Lh"])
  }
}

# Now, summarize the results
if (RJmodel == "MCMC" || RJmodel == "BOTH") {
  # Save the median independent rates
  summary["MCMC_Alpha_1", "Random"] <- median(QmatMat[, "MCMC_Alpha_1"])
  summary["MCMC_Alpha_2", "Random"] <- median(QmatMat[, "MCMC_Alpha_2"])
  summary["MCMC_Beta_1", "Random"] <- median(QmatMat[, "MCMC_Beta_1"])
  summary["MCMC_Beta_2", "Random"] <- median(QmatMat[, "MCMC_Beta_2"])
  # Save the median dependent rates
  summary["MCMC_q12", "Random"] <- median(QmatMat[, "MCMC_q12"])
  summary["MCMC_q13", "Random"] <- median(QmatMat[, "MCMC_q13"])
  summary["MCMC_q21", "Random"] <- median(QmatMat[, "MCMC_q21"])
  summary["MCMC_q24", "Random"] <- median(QmatMat[, "MCMC_q24"])
  summary["MCMC_q31", "Random"] <- median(QmatMat[, "MCMC_q31"])
  summary["MCMC_q34", "Random"] <- median(QmatMat[, "MCMC_q34"])
  summary["MCMC_q42", "Random"] <- median(QmatMat[, "MCMC_q42"])
  summary["MCMC_q43", "Random"] <- median(QmatMat[, "MCMC_q43"])
  # Save the median ESS's
  summary["MCMC_ESS_Ind", "Random"] <- median(QmatMat[, "MCMC_ESS_Ind"])
  summary["MCMC_ESS_Dep", "Random"] <- median(QmatMat[, "MCMC_ESS_Dep"])
}

if (RJmodel == "RJMCMC" || RJmodel == "BOTH") {
  # Save the median independent rates
  summary["RJ_Alpha_1", "Random"] <- median(QmatMat[, "RJ_Alpha_1"])
  summary["RJ_Alpha_2", "Random"] <- median(QmatMat[, "RJ_Alpha_2"])
  summary["RJ_Beta_1", "Random"] <- median(QmatMat[, "RJ_Beta_1"])
  summary["RJ_Beta_2", "Random"] <- median(QmatMat[, "RJ_Beta_2"])
  # Save the median dependent rates
  summary["RJ_q12", "Random"] <- median(QmatMat[, "RJ_q12"])
  summary["RJ_q13", "Random"] <- median(QmatMat[, "RJ_q13"])
  summary["RJ_q21", "Random"] <- median(QmatMat[, "RJ_q21"])
  summary["RJ_q24", "Random"] <- median(QmatMat[, "RJ_q24"])
  summary["RJ_q31", "Random"] <- median(QmatMat[, "RJ_q31"])
  summary["RJ_q34", "Random"] <- median(QmatMat[, "RJ_q34"])
  summary["RJ_q42", "Random"] <- median(QmatMat[, "RJ_q42"])
  summary["RJ_q43", "Random"] <- median(QmatMat[, "RJ_q43"])
  # Save the median ESS's
  summary["RJ_ESS_Ind", "Random"] <- median(QmatMat[, "RJ_ESS_Ind"])
  summary["RJ_ESS_Dep", "Random"] <- median(QmatMat[, "RJ_ESS_Dep"])
}

# Write the Qmat matrix
QmatName <- "Results/Random/QmatInfo/Random.QmatInfo.txt"
write.table(QmatMat, file = QmatName, quote = F, row.names = F, col.names = T, sep = "\t")

# Write the summary matrix
SummaryName <- paste0("Results/ConstantRates/Single/QmatInfo/Summary.QmatInfo.txt")
write.table(summary, file = SummaryName, quote = F, row.names = T, col.names = T, sep = "\t")







# Repeat if variable rates were also tested
if (isTRUE(variable_rates)) {
  for (type in types) {
    # Print a progress check
    print(paste0("Beginning to summarize estimated rates matrices from ", type, " tests with Variable Rates."))
    
    # Call the necessary paths
    results_name <- paste0("Results/VariableRates/Single/", type, ".Single.ResultsFull.txt")
    anccon_name <- paste0("Results/VariableRates/Single/AncestralStateReconstruction/", type, ".AncReconAndTreeHeight.txt")
    
    # Call the results table  
    colnames_results <- strsplit(readLines(results_name, n = 1), "\t")[[1]]
    results <- read.table(results_name, skip = 1, sep = "\t", col.names = colnames_results)
    
    # Call the Ancestral State Reconstruction Table
    colnames_anccon <- strsplit(readLines(anccon_name, n = 1), "\t")[[1]]
    AncCon <- read.table(anccon_name, skip = 1, sep = "\t", col.names = colnames_anccon)
    
    # Create a new matrix for each set of simulated data
    QmatMat <- matrix(data=NA, nrow = nrow(results), ncol = (10+length(summary_row_labels)))
    colnames(QmatMat) <- c("Trial_#", "Unknown_Taxon", "True_4States", "TBL", "Tree_Height", 
                           "Std_TBL", "n00", "n01", "n10", "n11", summary_row_labels)
    
    # Pull the information relevant to the rates matrices and store it in a new matrix
    QmatMat[, "Trial_#"] <- results[, "Trial"]
    QmatMat[, "Unknown_Taxon"] <- results[, "rTaxon"]
    QmatMat[, "True_4States"] <- results[, "True_4States"]
    QmatMat[, "n00"] <- results[, "n00"]
    QmatMat[, "n01"] <- results[, "n01"]
    QmatMat[, "n10"] <- results[, "n10"]
    QmatMat[, "n11"] <- results[, "n11"]
    QmatMat[, "TBL"] <- AncCon[, "TBL"]
    QmatMat[, "Tree_Height"] <- AncCon[, "Tree_Height"]
    QmatMat[, "Std_TBL"] <- AncCon[, "Std_TBL"]
    
    for (i in 1:nrow(results)) {
      if (RJmodel == "MCMC" || RJmodel == "BOTH") {
        # Gather the independent q-matrix
        Qmat_IndMCMC <- estimated_qmat(type = type,
                                       i = i,
                                       IndependentRates = TRUE,
                                       ConstantRates = FALSE,
                                       RJ = FALSE,
                                       Predict = FALSE)
        
        # Record the rates for this trial
        QmatMat[i, "MCMC_Alpha_1"] <- Qmat_IndMCMC[1, 3]
        QmatMat[i, "MCMC_Alpha_2"] <- Qmat_IndMCMC[1, 2] 
        QmatMat[i, "MCMC_Beta_1"] <- Qmat_IndMCMC[3, 1]
        QmatMat[i, "MCMC_Beta_2"] <- Qmat_IndMCMC[2, 1]
        
        # Pull the log file and record the ESS
        logname <- paste0("VariableRates/", type, "/Single/", type, ".Ind.MCMC.", i, ".Rates.Log.txt")
        logdata <- read.delim(logname, skip = 45)
        QmatMat[i, "MCMC_ESS_Ind"] <- ESS(logdata[, "Lh"])
        
        # Gather the dependent q-matrix
        Qmat_DepMCMC <- estimated_qmat(type = type,
                                       i = i,
                                       IndependentRates = FALSE,
                                       ConstantRates = FALSE,
                                       RJ = FALSE,
                                       Predict = FALSE)
        
        # Record the rates for this trial
        QmatMat[i, "MCMC_q12"] <- Qmat_DepMCMC[1, 2]
        QmatMat[i, "MCMC_q13"] <- Qmat_DepMCMC[1, 3] 
        QmatMat[i, "MCMC_q21"] <- Qmat_DepMCMC[2, 1]
        QmatMat[i, "MCMC_q24"] <- Qmat_DepMCMC[2, 4]
        QmatMat[i, "MCMC_q31"] <- Qmat_DepMCMC[3, 1]
        QmatMat[i, "MCMC_q34"] <- Qmat_DepMCMC[3, 4] 
        QmatMat[i, "MCMC_q42"] <- Qmat_DepMCMC[4, 2]
        QmatMat[i, "MCMC_q43"] <- Qmat_DepMCMC[4, 3]
        
        # Pull the log file and record the ESS
        logname <- paste0("VariableRates/", type, "/Single/", type, ".Dep.MCMC.", i, ".Rates.Log.txt")
        logdata <- read.delim(logname, skip = 53)
        QmatMat[i, "MCMC_ESS_Dep"] <- ESS(logdata[, "Lh"])
      }
      
      if (RJmodel == "RJMCMC" || RJmodel == "BOTH") {
        # Gather the independent q-matrix
        Qmat_IndRJ <- estimated_qmat(type = type,
                                     i = i,
                                     IndependentRates = TRUE,
                                     ConstantRates = FALSE,
                                     RJ = TRUE,
                                     Predict = FALSE)
        
        # Record the rates for this trial
        QmatMat[i, "RJ_Alpha_1"] <- Qmat_IndRJ[1, 3]
        QmatMat[i, "RJ_Alpha_2"] <- Qmat_IndRJ[1, 2] 
        QmatMat[i, "RJ_Beta_1"] <- Qmat_IndRJ[3, 1]
        QmatMat[i, "RJ_Beta_2"] <- Qmat_IndRJ[2, 1]
        
        # Pull the log file and record the ESS
        logname <- paste0("VariableRates/", type, "/Single/", type, ".Ind.RJMCMC.", i, ".Rates.Log.txt")
        logdata <- read.delim(logname, skip = 44)
        QmatMat[i, "RJ_ESS_Ind"] <- ESS(logdata[, "Lh"])
        
        # Gather the dependent q-matrix
        Qmat_DepRJ <- estimated_qmat(type = type,
                                     i = i,
                                     IndependentRates = FALSE,
                                     ConstantRates = FALSE,
                                     RJ = TRUE,
                                     Predict = FALSE)
        
        # Record the rates for this trial
        QmatMat[i, "RJ_q12"] <- Qmat_DepRJ[1, 2]
        QmatMat[i, "RJ_q13"] <- Qmat_DepRJ[1, 3] 
        QmatMat[i, "RJ_q21"] <- Qmat_DepRJ[2, 1]
        QmatMat[i, "RJ_q24"] <- Qmat_DepRJ[2, 4]
        QmatMat[i, "RJ_q31"] <- Qmat_DepRJ[3, 1]
        QmatMat[i, "RJ_q34"] <- Qmat_DepRJ[3, 4] 
        QmatMat[i, "RJ_q42"] <- Qmat_DepRJ[4, 2]
        QmatMat[i, "RJ_q43"] <- Qmat_DepRJ[4, 3]
        
        # Pull the log file and record the ESS
        logname <- paste0("VariableRates/", type, "/Single/", type, ".Dep.RJMCMC.", i, ".Rates.Log.txt")
        logdata <- read.delim(logname, skip = 48)
        QmatMat[i, "RJ_ESS_Dep"] <- ESS(logdata[, "Lh"])
      }
    }
    
    # Now, summarize the results
    if (RJmodel == "MCMC" || RJmodel == "BOTH") {
      # Save the median independent rates
      summary["MCMC_Alpha_1", type] <- median(QmatMat[, "MCMC_Alpha_1"])
      summary["MCMC_Alpha_2", type] <- median(QmatMat[, "MCMC_Alpha_2"])
      summary["MCMC_Beta_1", type] <- median(QmatMat[, "MCMC_Beta_1"])
      summary["MCMC_Beta_2", type] <- median(QmatMat[, "MCMC_Beta_2"])
      # Save the median dependent rates
      summary["MCMC_q12", type] <- median(QmatMat[, "MCMC_q12"])
      summary["MCMC_q13", type] <- median(QmatMat[, "MCMC_q13"])
      summary["MCMC_q21", type] <- median(QmatMat[, "MCMC_q21"])
      summary["MCMC_q24", type] <- median(QmatMat[, "MCMC_q24"])
      summary["MCMC_q31", type] <- median(QmatMat[, "MCMC_q31"])
      summary["MCMC_q34", type] <- median(QmatMat[, "MCMC_q34"])
      summary["MCMC_q42", type] <- median(QmatMat[, "MCMC_q42"])
      summary["MCMC_q43", type] <- median(QmatMat[, "MCMC_q43"])
      # Save the median ESS's
      summary["MCMC_ESS_Ind", type] <- median(QmatMat[, "MCMC_ESS_Ind"])
      summary["MCMC_ESS_Dep", type] <- median(QmatMat[, "MCMC_ESS_Dep"])
    }
    
    if (RJmodel == "RJMCMC" || RJmodel == "BOTH") {
      # Save the median independent rates
      summary["RJ_Alpha_1", type] <- median(QmatMat[, "RJ_Alpha_1"])
      summary["RJ_Alpha_2", type] <- median(QmatMat[, "RJ_Alpha_2"])
      summary["RJ_Beta_1", type] <- median(QmatMat[, "RJ_Beta_1"])
      summary["RJ_Beta_2", type] <- median(QmatMat[, "RJ_Beta_2"])
      # Save the median dependent rates
      summary["RJ_q12", type] <- median(QmatMat[, "RJ_q12"])
      summary["RJ_q13", type] <- median(QmatMat[, "RJ_q13"])
      summary["RJ_q21", type] <- median(QmatMat[, "RJ_q21"])
      summary["RJ_q24", type] <- median(QmatMat[, "RJ_q24"])
      summary["RJ_q31", type] <- median(QmatMat[, "RJ_q31"])
      summary["RJ_q34", type] <- median(QmatMat[, "RJ_q34"])
      summary["RJ_q42", type] <- median(QmatMat[, "RJ_q42"])
      summary["RJ_q43", type] <- median(QmatMat[, "RJ_q43"])
      # Save the median ESS's
      summary["RJ_ESS_Ind", type] <- median(QmatMat[, "RJ_ESS_Ind"])
      summary["RJ_ESS_Dep", type] <- median(QmatMat[, "RJ_ESS_Dep"])
    }
    
    # Write the Qmat matrix
    QmatName <- paste0("Results/VariableRates/Single/QmatInfo/", type, ".QmatInfo.txt")
    write.table(QmatMat, file = QmatName, quote = F, row.names = F, col.names = T, sep = "\t")
  }
  
  
  # Now, pull the random info
  # Print a progress check
  print("Beginning to summarize estimated rates matrices from tests with Random Rates.")
  
  # Get the data
  CR_info <- read.delim(file = "Results/ConstantRates/Single/QmatInfo/Summary.QmatInfo.txt", header = TRUE)
  
  
  # Move the random data into the VR table
  summary[, "Random"] <- CR_info[, "Random"]
  
  
  # Write the summary matrix
  SummaryName <- paste0("Results/VariableRates/Single/QmatInfo/Summary.QmatInfo.txt")
  write.table(summary, file = SummaryName, quote = F, row.names = T, col.names = T, sep = "\t")
}