### This script was written to be able to compare corHMM to the rest of the simulation

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





### This is where the real code will start
### For loop to record all of this information for Constant Rates
for (type in types) {
  # Print a progress check
  print(paste0("Beginning the ancestral reconstruction for ", type, " with Constant Rates."))
  
  # Call the necessary paths
  data_path <- paste0("ConstantRates/", type, "/Data/", type)
  results_name <- paste0("Results/ConstantRates/Single/", type, ".Single.ResultsFull.txt")
  AncRecon_name <- paste0("Results/ConstantRates/Single/AncestralStateReconstruction/", type, ".AncReconAndTreeHeight.txt")
  
  # Call the Results and Ancestral State Reconstruction tables 
  colnames_results <- strsplit(readLines(results_name, n = 1), "\t")[[1]]
  results <- read.table(results_name, skip = 1, sep = "\t", col.names = colnames_results)
  
  colnames_AncRecon <- strsplit(readLines(AncRecon_name, n = 1), "\t")[[1]]
  AncConMat <- read.table(AncRecon_name, skip = 1, sep = "\t", col.names = colnames_AncRecon)
  
  # Create a matrix to store the information
  matrix_name <- paste0("Results/ConstantRates/Single/corHMM/", type, ".corHmmResults.txt")
  
  corMat <- matrix(data=NA, nrow = trial_amount, ncol = 28)
  colnames(corMat) <- c("Trial_#", "Unknown_Taxon", "True_A", "True_B", "True_4States", "corHMM_Pred_4States", "corHMM_PredProb",
                        "corHMM_acc", "corHMM_LL", "True_Is_Max?", "PredProb_of_True", "corHMM_Lh_00", "corHMM_Lh_01", "corHMM_Lh_10", "corHMM_Lh_11",
                        "Est_q12", "Est_q13", "Est_q21", "Est_q24", "Est_q31", "Est_q34", "Est_q42", "Est_q43", "terminal.std", "n00", "n01", "n10", "n11")


  # Start the for-loop that will run each prediction and fill it in
  for (i in 1:trial_amount) {
    # First, call rTaxon by subsetting the rows where column 1 matches i
    matching_rows <- which(results[, "Trial"] == i)
    corMat[matching_rows, "Trial_#"] <- i
    corMat[matching_rows, "Unknown_Taxon"] <- rTaxon <- results[matching_rows, "rTaxon"]
    corMat[matching_rows, "True_A"] <- True_A <- results[matching_rows, "Trait_A"]
    corMat[matching_rows, "True_B"] <- True_B <- results[matching_rows, "Trait_B"]
    corMat[matching_rows, "True_4States"] <- True_4 <- results[matching_rows, "True_4States"]
    
    # Call the data file for this iteration
    data_name <- paste0(data_path, ".", i, ".Full_data.txt")
    data <- read.table(data_name, skip = 1, sep = "\t")
    
    # Set up the data for prediction
    predict_data_AB <- data[, 1:3]
    predict_data_AB[, 2] <- as.numeric(predict_data_AB[, 2]) 
    predict_data_AB[, 3] <- as.numeric(predict_data_AB[, 3])
    predict_data_AB[rTaxon, 3] <- "?"
  
  
    # Pull the tree to read the terminal branch length
    treename <- paste0("Trees/Full_tree.", i, ".tre")
    tree_complete <- read.nexus(treename)  
    
    
    
    # The corHMM function breaks if n10 and/or n11 is 0 & we're predicting for "0_"
    ### So this if statement adjusts the code as necessary to prevent breakage
    if (True_A == 0 && results[matching_rows, "n10"] == 0 && results[matching_rows, "n11"] == 0) {
      # Set up a matrix for when n10 = 0
      blank_qmat <- matrix(c(
        NA,1,
        2,NA),
        2,2,byrow=TRUE)
      
      # This is the full run of the corHMM software
      corhmm_full <- corHMM(tree_complete, predict_data_AB, rate.mat = blank_qmat,
                            rate.cat = 1, node.states = "marginal", get.tip.states = TRUE)
      
      # Adjust the calculated Lhs for 4 states instead of 3
      corHMM_Lhs <- round(corhmm_full$tip.states[rTaxon, ], digits = 6)
      corHMM_Lhs[4] <- corHMM_Lhs[3] <- 0
      
      
      # Pull the prediction info from the run
      corMat[matching_rows, c("corHMM_Lh_00", "corHMM_Lh_01", "corHMM_Lh_10", "corHMM_Lh_11")] <- corHMM_Lhs
      
      # Save the estimated q-matrix information
      corMat[matching_rows, "Est_q12"] <- corhmm_full$solution[1, 2]
      corMat[matching_rows, "Est_q21"] <- corhmm_full$solution[2, 1]
      corMat[matching_rows, c("Est_q13", "Est_q24", "Est_q31","Est_q34", 
                              "Est_q42", "Est_q43")] <- 0
    
      
      ### Repeat the check, but for trials with no "10" tips only 
    } else if (True_A == 0 && results[matching_rows, "n10"] == 0) {
      
      # Set up a matrix for when n10 = 0
      blank_qmat <- matrix(c(
        NA,1,0,
        2,NA,3,
        0,4,NA),
        3,3,byrow=TRUE)
      
      # This is the full run of the corHMM software
      corhmm_full <- corHMM(tree_complete, predict_data_AB, rate.mat = blank_qmat,
                            rate.cat = 1, node.states = "marginal", get.tip.states = TRUE)
      
      # Adjust the calculated Lhs for 4 states instead of 3
      corHMM_Lhs <- round(corhmm_full$tip.states[rTaxon, ], digits = 6)
      corHMM_Lhs[4] <- corHMM_Lhs[3]
      corHMM_Lhs[3] <- 0
      
      # Pull the prediction info from the run
      corMat[matching_rows, c("corHMM_Lh_00", "corHMM_Lh_01", "corHMM_Lh_10", "corHMM_Lh_11")] <- corHMM_Lhs
      
      # Save the estimated q-matrix information
      corMat[matching_rows, "Est_q12"] <- corhmm_full$solution[1, 2]
      corMat[matching_rows, c("Est_q21", "Est_q24")] <- corhmm_full$solution[2, c(1,3)]
      corMat[matching_rows, "Est_q42"] <- corhmm_full$solution[3, 2]
      corMat[matching_rows, c("Est_q13","Est_q31", "Est_q34", "Est_q43")] <- 0
      
      
      ### Repeat the check, but for trials with no "11" tips only
    } else if (True_A == 0 && results[matching_rows, "n11"] == 0) {
      
      # Set up a matrix for when n11 = 0
      blank_qmat <- matrix(c(
        NA,1,2,
        3,NA,0,
        4,0,NA),
        3,3,byrow=TRUE)
      
      # This is the full run of the corHMM software
      corhmm_full <- corHMM(tree_complete, predict_data_AB, rate.mat = blank_qmat,
                            rate.cat = 1, node.states = "marginal", get.tip.states = TRUE)
      
      # Adjust the calculated Lhs for 4 states instead of 3
      corHMM_Lhs <- round(corhmm_full$tip.states[rTaxon, ], digits = 6)
      corHMM_Lhs[4] <- 0
      
      # Pull the prediction info from the run
      corMat[matching_rows, c("corHMM_Lh_00", "corHMM_Lh_01", "corHMM_Lh_10", "corHMM_Lh_11")] <- corHMM_Lhs
      
      # Save the estimated q-matrix information
      corMat[matching_rows, c("Est_q12", "Est_q13")] <- corhmm_full$solution[1, 2:3]
      corMat[matching_rows, "Est_q21"] <- corhmm_full$solution[2, 1]
      corMat[matching_rows, "Est_q31"] <- corhmm_full$solution[3, 1]
      corMat[matching_rows, c("Est_q24", "Est_q34", "Est_q42", "Est_q43")] <- 0
      
      
      
      ### This is the standard process for all other trials
    } else {
      
      # Set up the proper q-matrix
      blank_qmat <- matrix(c(
        NA,1,2,0,
        3,NA,0,4,
        5,0,NA,6,
        0,7,8,NA),
        4,4,byrow=TRUE)
      
      ### This is the full run of the corHMM software
      corhmm_full <- corHMM(tree_complete, predict_data_AB, rate.mat = blank_qmat,
                            rate.cat = 1, node.states = "marginal", get.tip.states = TRUE)
      
      # This pulls the prediction info from the run
      corMat[matching_rows, c("corHMM_Lh_00", "corHMM_Lh_01", "corHMM_Lh_10", "corHMM_Lh_11")] <- 
        corHMM_Lhs <- round(corhmm_full$tip.states[rTaxon, ], digits = 6)
      corMat[matching_rows, "corHMM_Pred_4States"] <- corHMM_pred <- which.max(corhmm_full$tip.states[rTaxon, ])
      corMat[matching_rows, "corHMM_PredProb"] <- corHMM_pred_prob <- corHMM_Lhs[corHMM_pred]
      
      # Save the estimated q-matrix information
      corMat[matching_rows, c("Est_q12", "Est_q13")] <- corhmm_full$solution[1, 2:3]
      corMat[matching_rows, c("Est_q21", "Est_q24")] <- corhmm_full$solution[2, c(1,4)]
      corMat[matching_rows, c("Est_q31", "Est_q34")] <- corhmm_full$solution[3, c(1,4)]
      corMat[matching_rows, c("Est_q42", "Est_q43")] <- corhmm_full$solution[4, 2:3]
    }
    
    
    
    ### Pull the prediction info from the run
    corMat[matching_rows, "corHMM_Pred_4States"] <- corHMM_pred <- which.max(corhmm_full$tip.states[rTaxon, ])
    corMat[matching_rows, "corHMM_PredProb"] <- corHMM_pred_prob <- corHMM_Lhs[corHMM_pred]
    
    # Adjust predicted probability to the probability of a 1
    if(corHMM_pred == 1 || corHMM_pred == 3) {corHMM_pred_prob <- 1 - corHMM_pred_prob}
    
    # Calculate some quality metrics
    corMat[matching_rows, "corHMM_acc"] <- calculate_accuracy(True_B, corHMM_pred_prob) 
    corMat[matching_rows, "corHMM_LL"] <- LogLoss(True_B, corHMM_pred_prob)
    
    # Determine how the model did for the true value
    corMat[matching_rows, "True_Is_Max?"] <- ifelse(corHMM_pred == corMat[matching_rows, "True_4States"], 1, 0)
    corMat[matching_rows, "PredProb_of_True"] <- corHMM_Lhs[corMat[matching_rows, "True_4States"]]
    
    # Save the tree and population information
    corMat[matching_rows, "terminal.std"] <- AncConMat[matching_rows, "Std_TBL"]
    corMat[matching_rows, "n00"] <- results[matching_rows, "n00"]
    corMat[matching_rows, "n01"] <- results[matching_rows, "n01"]
    corMat[matching_rows, "n10"] <- results[matching_rows, "n10"]
    corMat[matching_rows, "n11"] <- results[matching_rows, "n11"]
  }

  # Save this table
  write.table(corMat, file = matrix_name, quote = F, row.names = F, col.names = T, sep = "\t")


  # Finally, send a progress signal
  print(paste0("Finished with the ancestral reconstructions for ", type, " with Constant Rates"))
} 

###
# Random Trials
#####
### Repeat for Randomly generated Data
# Print a progress check
print(paste0("Beginning the ancestral reconstruction for Random with Constant Rates."))

# Call the necessary paths
data_path <- "Random/Data/Random"
results_name <- "Results/Random/Random.Single.ResultsFull.txt"

# Call the results table  
colnames_results <- strsplit(readLines(results_name, n = 1), "\t")[[1]]
results <- read.table(results_name, skip = 1, sep = "\t", col.names = colnames_results)

# Create a matrix to store the information
matrix_name <- paste0("Results/Random/corHMM/Random.corHmmResults.txt")

corMat <- matrix(data=NA, nrow = trial_amount, ncol = 23)
# corMat <- matrix(data=NA, nrow = nrow(results), ncol = 23)
colnames(corMat) <- c("Trial_#", "Unknown_Taxon", "True_A", "True_B", "True_4States", "corHMM_Pred_4States", "corHMM_PredProb",
                      "corHMM_acc", "corHMM_LL", "True_Is_Max?", "PredProb_of_True", "corHMM_Lh_00", "corHMM_Lh_01", "corHMM_Lh_10", "corHMM_Lh_11",
                      "Est_q12", "Est_q13", "Est_q21", "Est_q24", "Est_q31", "Est_q34", "Est_q42", "Est_q43")


# Start the for-loop that will run each prediction and fill it in
for (i in 1:trial_amount) {
  # First, call rTaxon by subsetting the rows where column 1 matches i
  matching_rows <- which(results[, "Trial"] == i)
  corMat[matching_rows, "Trial_#"] <- i
  corMat[matching_rows, "Unknown_Taxon"] <- rTaxon <- results[matching_rows, "rTaxon"]
  corMat[matching_rows, "True_A"] <- True_A <- results[matching_rows, "Trait_A"]
  corMat[matching_rows, "True_B"] <- True_B <- results[matching_rows, "Trait_B"]
  corMat[matching_rows, "True_4States"] <- True_4 <- results[matching_rows, "True_4States"]
  
  # Call the data file for this iteration
  data_name <- paste0(data_path, ".", i, ".Full_data.txt")
  data <- read.table(data_name, skip = 1, sep = "\t")
  
  # Set up the data for prediction
  predict_data_AB <- data[, 1:3]
  predict_data_AB[, 2] <- as.numeric(predict_data_AB[, 2]) 
  predict_data_AB[, 3] <- as.numeric(predict_data_AB[, 3])
  predict_data_AB[rTaxon, 3] <- "?"
  
  
  # Pull the tree to read the terminal branch length
  treename <- paste0("Trees/Full_tree.", i, ".tre")
  tree_complete <- read.nexus(treename)  
  
  
  
  # Set up the proper q-matrix
  blank_qmat <- matrix(c(
    NA,1,2,0,
    3,NA,0,4,
    5,0,NA,6,
    0,7,8,NA),
    4,4,byrow=TRUE)
  
  ### This is the full run of the corHMM software
  corhmm_full <- corHMM(tree_complete, predict_data_AB, rate.mat = blank_qmat,
                        rate.cat = 1, node.states = "marginal", get.tip.states = TRUE)
  
  # This pulls the prediction info from the run
  corMat[matching_rows, c("corHMM_Lh_00", "corHMM_Lh_01", "corHMM_Lh_10", "corHMM_Lh_11")] <- 
    corHMM_Lhs <- round(corhmm_full$tip.states[rTaxon, ], digits = 6)
  corMat[matching_rows, "corHMM_Pred_4States"] <- corHMM_pred <- which.max(corhmm_full$tip.states[rTaxon, ])
  corMat[matching_rows, "corHMM_PredProb"] <- corHMM_pred_prob <- corHMM_Lhs[corHMM_pred]
  
  # Save the estimated q-matrix information
  corMat[matching_rows, c("Est_q12", "Est_q13")] <- corhmm_full$solution[1, 2:3]
  corMat[matching_rows, c("Est_q21", "Est_q24")] <- corhmm_full$solution[2, c(1,4)]
  corMat[matching_rows, c("Est_q31", "Est_q34")] <- corhmm_full$solution[3, c(1,4)]
  corMat[matching_rows, c("Est_q42", "Est_q43")] <- corhmm_full$solution[4, 2:3]
  
  
  ### Pull the prediction info from the run
  corMat[matching_rows, "corHMM_Pred_4States"] <- corHMM_pred <- which.max(corhmm_full$tip.states[rTaxon, ])
  corMat[matching_rows, "corHMM_PredProb"] <- corHMM_pred_prob <- corHMM_Lhs[corHMM_pred]
  
  # Adjust predicted probability to the probability of a 1
  if(corHMM_pred == 1 || corHMM_pred == 3) {corHMM_pred_prob <- 1 - corHMM_pred_prob}
  
  # Calculate some quality metrics
  corMat[matching_rows, "corHMM_acc"] <- calculate_accuracy(True_B, corHMM_pred_prob) 
  corMat[matching_rows, "corHMM_LL"] <- LogLoss(True_B, corHMM_pred_prob)
  
  # Determine how the model did for the true value
  corMat[matching_rows, "True_Is_Max?"] <- ifelse(corHMM_pred == corMat[matching_rows, "True_4States"], 1, 0)
  corMat[matching_rows, "PredProb_of_True"] <- corHMM_Lhs[corMat[matching_rows, "True_4States"]]
    
  # Save the tree and population information
  corMat[matching_rows, "terminal.std"] <- AncConMat[matching_rows, "Std_TBL"]
  corMat[matching_rows, "n00"] <- results[matching_rows, "n00"]
  corMat[matching_rows, "n01"] <- results[matching_rows, "n01"]
  corMat[matching_rows, "n10"] <- results[matching_rows, "n10"]
  corMat[matching_rows, "n11"] <- results[matching_rows, "n11"]
}

# Save this table
write.table(corMat, file = matrix_name, quote = F, row.names = F, col.names = T, sep = "\t")


# Finally, send a progress signal
print(paste0("Finished with the ancestral reconstructions for Random Rates"))
