# This script is for creating the full datatable for each run



### First, we need to pull the necessary information
# Call the bash file
bash_file <- Sys.glob("Discrete_Simulation.V*")
version_str <- sub(".*(V[0-9]+).*", "\\1", bash_file)
source(paste0("Scripts/DiscreteFunctions.", version_str, ".R"))

# Get the types and trials arrays
types <- grep("^types=", readLines(bash_file), value = TRUE)
types <- gsub("types=\\(|\\)", "", types)
types <- gsub("\"", "", types)
types <- strsplit(types, " ")[[1]]

# Get the special cases you may run
variable_rates <- read_logical("variable_rates", bash_file)
multiple_prediction <- read_logical("multiple_prediction", bash_file)
clade_prediction <- read_logical("clade_prediction", bash_file)
multistate_prediction <- read_logical("multistate_prediction", bash_file)

# Get the MCMC type
RJmodel <- grep("^RJmodel=", readLines(bash_file), value = TRUE)
RJmodel <- gsub(".*=(.*)", "\\1", RJmodel)

# Number of iterations when running as batch job/ on HPC
trial_amount <- grep("^num_iterations=", readLines(bash_file), value = TRUE)
trial_amount <- as.numeric(gsub("[^0-9.-]", "", trial_amount))



### Second, we need to know how big the results table will be
# We start with what we know (rTaxon info will be added later)
tests <- c("Beta_Binom", "Naive_Bayes")
test_size <- 2

# Expand for the BayesTraits trials
if (RJmodel == "MCMC") {
  BTruns <- c("MCMC_Ind", "MCMC_Dep")
  BTsize <- 2
  ACR_labels <- c("Sister_IndMCMC_MaxLh", "Sister_IndMCMC_Prediction", "Sister_DepMCMC_MaxLh", "Sister_DepMCMC_Prediction",
               "Ancestor_IndMCMC_MaxLh", "Ancestor_IndMCMC_Prediction","Ancestor_DepMCMC_MaxLh", "Ancestor_DepMCMC_Prediction",
               "TrueState_Sister_IndMCMC_Lh", "TrueState_Sister_DepMCMC_Lh", "TrueState_Ancestor_IndMCMC_Lh", "TrueState_Ancestor_DepMCMC_Lh")
  ACRsize <- 12
  if (isTRUE(multistate_prediction)) {
    BTruns <- c("MCMC_Multi", BTruns)
    BTsize <- BTsize + 1
  }
}

if (RJmodel == "RJMCMC") {
  BTruns <- c("RJ_Ind", "RJ_Dep")
  BTsize <- 2
  ACR_labels <- c("Sister_IndRJ_MaxLh", "Sister_IndRJ_Prediction", "Sister_DepRJ_MaxLh", "Sister_DepRJ_Prediction",
               "Ancestor_IndRJ_MaxLh", "Ancestor_IndRJ_Prediction","Ancestor_DepRJ_MaxLh", "Ancestor_DepRJ_Prediction",
               "TrueState_Sister_IndRJ_Lh", "TrueState_Sister_DepRJ_Lh", "TrueState_Ancestor_IndRJ_Lh", "TrueState_Ancestor_DepRJ_Lh")
  ACRsize <- 12
  if (isTRUE(multistate_prediction)) {
    BTruns <- c("RJ_Multi", BTruns)
    BTsize <- BTsize + 1
  }
}

if (RJmodel == "BOTH") {
  BTruns <- c("MCMC_Ind", "MCMC_Dep", "RJ_Ind", "RJ_Dep")
  BTsize <- 4
  ACR_labels <- c("Sister_IndMCMC_MaxLh", "Sister_IndMCMC_Prediction", "Sister_DepMCMC_MaxLh", "Sister_DepMCMC_Prediction",
               "Sister_IndRJ_MaxLh", "Sister_IndRJ_Prediction", "Sister_DepRJ_MaxLh", "Sister_DepRJ_Prediction",
               "Ancestor_IndMCMC_MaxLh", "Ancestor_IndMCMC_Prediction","Ancestor_DepMCMC_MaxLh", "Ancestor_DepMCMC_Prediction",
               "Ancestor_IndRJ_MaxLh", "Ancestor_IndRJ_Prediction","Ancestor_DepRJ_MaxLh", "Ancestor_DepRJ_Prediction",
               "TrueState_Sister_IndMCMC_Lh", "TrueState_Sister_DepMCMC_Lh", "TrueState_Sister_IndRJ_Lh", "TrueState_Sister_DepRJ_Lh",
               "TrueState_Ancestor_IndMCMC_Lh", "TrueState_Ancestor_DepMCMC_Lh", "TrueState_Ancestor_IndRJ_Lh", "TrueState_Ancestor_DepRJ_Lh")
  ACRsize <- 24
  if (isTRUE(multistate_prediction)) {
    BTruns <- c("MCMC_Multi", "MCMC_Ind", "MCMC_Dep", "RJ_Multi", "RJ_Ind", "RJ_Dep")
    BTsize <- BTsize + 2
  }
}

# Determine the total number of models being tested and their names
tests <- c(tests, BTruns)
test_size <- test_size + BTsize

# Now, we separate out each into three columns
prob_labels <- paste0(tests, "_Prob")
acc_labels <- paste0(tests, "_acc")
LL_labels <- paste0(tests, "_LL")

# Lastly, we want to include the state frequency counts
state_labels <- c("True_4States","n00", "n01", "n10", "n11", "Avg_Sister_A", "Avg_Sister_B")

# We readjust our new column names and column count to match
new_col_names <- c(state_labels, prob_labels, acc_labels, LL_labels, ACR_labels)
new_col_number <- (test_size * 3) + 7 + ACRsize



#############################
### Now we create the new matrices for the data
# First we go through the settings we know have run
for (type in types) {
  # First, pull the unknown matrix
  unknown_labels <- c("Trial", "rTaxon", "Trait_A", "Trait_B", "Terminal_Branch_Length")
  unknown_name <- paste0("ConstantRates/", type, "/Single/", type, ".Unknown_info.txt")
  unknown_info <- read.table(unknown_name, skip = 1, col.names = unknown_labels)
  row_num <- nrow(unknown_info)
  full_data_path <- paste0("ConstantRates/", type, "/Data/", type)

  # Next, create the new section of the results matrix
  matrix_b <- matrix(data=NA, nrow = row_num, ncol = new_col_number)
  colnames(matrix_b) <- new_col_names

  # Combine the unknown information with the new space for results
  results <- cbind(unknown_info, matrix_b)

  # Fill in the counts data that has already by calculated
  for (i in 1:trial_amount) {
    # Find the rows that match the trial number
    trial_i <- which(results[, 1] == i)
    matching_rows <- which(results[, "Trial"] == i)

    # Pull the true state data
    full_data_name <- paste0(full_data_path, ".", i, ".Full_data.txt")
    trait_data <- read.table(full_data_name, skip = 1, sep = "\t")
    trait_data[,4] <- match(trait_data[,4], c(0,1,10,11))

    # Put the 4-state character into the results matrix
    rTaxon <- results[matching_rows, "rTaxon"]
    for (j in seq_along(rTaxon)) {results[matching_rows[j], "True_4States"] <- trait_data[rTaxon[j], 4]}

    # Call the count data and sort it into a vector
    count_name <- paste0("ConstantRates/", type, "/Data/", type, ".", i, ".Counts.txt")
    counts <- read.table(count_name, skip = 1)
    count_vector <- c(counts[1, 2], counts[2, 2], counts[1, 3], counts[2, 3])

    # Remove the unknown taxa and store it for later
    truestate_vector <- results[matching_rows, "True_4States"]
    for (j in seq_along(rTaxon)) {count_vector[truestate_vector[j]] <- count_vector[truestate_vector[j]] - 1}
    results[trial_i, c("n00", "n01", "n10", "n11")] <- count_vector
  }

  # Finally, write the matrix somewhere it can be accessed later
  results_name <- paste0("Results/ConstantRates/Single/", type, ".Single.ResultsFull.txt")
  write.table(results, file = results_name, quote = F, sep = "\t", row.names = F, col.names = T)
}



### Now we do this for the Random Data
# First, pull the unknown matrix
unknown_labels <- c("Trial", "rTaxon", "Trait_A", "Trait_B", "Terminal_Branch_Length")
unknown_name <- "Random/Single/Random.Unknown_info.txt"
unknown_info <- read.table(unknown_name, skip = 1, col.names = unknown_labels)
row_num <- nrow(unknown_info)
full_data_path <- "Random/Data/Random"

# Next, create the new section of the results matrix
matrix_b <- matrix(data=NA, nrow = row_num, ncol = new_col_number)
colnames(matrix_b) <- new_col_names

# Combine the unknown information with the new space for results
results <- cbind(unknown_info, matrix_b)

# Fill in the counts data that has already by calculated
for (i in 1:trial_amount) {
  # Find the rows that match the trial number
  trial_i <- which(results[, 1] == i)
  matching_rows <- which(results[, "Trial"] == i)

  # Pull the true state data
  full_data_name <- paste0(full_data_path, ".", i, ".Full_data.txt")
  trait_data <- read.table(full_data_name, skip = 1, sep = "\t")
  trait_data[,4] <- match(trait_data[,4], c(0,1,10,11))

  # Put the 4-state character into the results matrix
  rTaxon <- results[matching_rows, "rTaxon"]
  for (j in seq_along(rTaxon)) {results[matching_rows[j], "True_4States"] <- trait_data[rTaxon[j], 4]}

  # Call the count data and sort it into a vector
  count_name <- paste0("Random/Data/Random.", i, ".Counts.txt")
  counts <- read.table(count_name, skip = 1)
  count_vector <- c(counts[1, 2], counts[2, 2], counts[1, 3], counts[2, 3])

  # Remove the unknown taxa and store it for later
  truestate_vector <- results[matching_rows, "True_4States"]
  for (j in seq_along(rTaxon)) {count_vector[truestate_vector[j]] <- count_vector[truestate_vector[j]] - 1}
  results[trial_i, c("n00", "n01", "n10", "n11")] <- count_vector
}

# Finally, write the matrix somewhere it can be accessed later
results_name <- "Results/Random/Random.Single.ResultsFull.txt"
write.table(results, file = results_name, quote = F, sep = "\t", row.names = F, col.names = T)



#############################
### Below are the reruns for all of the possible runs
# First of these is single prediction, variable rates
if (variable_rates == TRUE) {
  for (type in types) {
    # First, pull the unknown matrix
    unknown_labels <- c("Trial", "rTaxon", "Trait_A", "Trait_B", "Terminal_Branch_Length")
    unknown_name <- paste0("VariableRates/", type, "/Single/", type, ".Unknown_info.txt")
    unknown_info <- read.table(unknown_name, skip = 1, col.names = unknown_labels)
    row_num <- nrow(unknown_info)
    full_data_path <- paste0("VariableRates/", type, "/Data/", type)

    # Next, create the new section of the results matrix
    matrix_b <- matrix(data=NA, nrow = row_num, ncol = new_col_number)
    colnames(matrix_b) <- new_col_names

    # Combine the unknown information with the new space for results
    results <- cbind(unknown_info, matrix_b)

    # Fill in the counts data that has already by calculated
    for (i in 1:trial_amount) {
      # Find the rows that match the trial number
      trial_i <- which(results[, 1] == i)
      matching_rows <- which(results[, "Trial"] == i)

      # Pull the true state data
      full_data_name <- paste0(full_data_path, ".", i, ".Full_data.txt")
      trait_data <- read.table(full_data_name, skip = 1, sep = "\t")
      trait_data[,4] <- match(trait_data[,4], c(0,1,10,11))

      # Put the 4-state character into the results matrix
      rTaxon <- results[matching_rows, "rTaxon"]
      for (j in seq_along(rTaxon)) {results[matching_rows[j], "True_4States"] <- trait_data[rTaxon[j], 4]}

      # Call the count data and sort it into a vector
      count_name <- paste0("VariableRates/", type, "/Data/", type, ".", i, ".Counts.txt")
      counts <- read.table(count_name, skip = 1)
      count_vector <- c(counts[1, 2], counts[2, 2], counts[1, 3], counts[2, 3])

      # Remove the unknown taxa and store it for later
      truestate_vector <- results[matching_rows, "True_4States"]
      for (j in seq_along(rTaxon)) {count_vector[truestate_vector[j]] <- count_vector[truestate_vector[j]] - 1}
      results[trial_i, c("n00", "n01", "n10", "n11")] <- count_vector
    }

    # Finally, write the matrix somewhere it can be accessed later
    results_name <- paste0("Results/VariableRates/Single/", type, ".Single.ResultsFull.txt")
    write.table(results, file = results_name, quote = F, sep = "\t", row.names = F, col.names = T)
  }
}


#############################
### Now, check for multiple prediction for constant rates
if (multiple_prediction == TRUE) {
  for (type in types) {
    # First, pull the unknown matrix
    unknown_labels <- c("Trial", "rTaxon", "Trait_A", "Trait_B", "Terminal_Branch_Length")
    unknown_name <- paste0("ConstantRates/", type, "/Multiple/", type, ".Unknown_info.txt")
    unknown_info <- read.table(unknown_name, skip = 1, col.names = unknown_labels)
    row_num <- nrow(unknown_info)
    full_data_path <- paste0("ConstantRates/", type, "/Data/", type)

    # Next, create the new section of the results matrix
    matrix_b <- matrix(data=NA, nrow = row_num, ncol = new_col_number)
    colnames(matrix_b) <- new_col_names

    # Combine the unknown information with the new space for results
    results <- cbind(unknown_info, matrix_b)

    # Fill in the counts data that has already by calculated
    for (i in 1:trial_amount) {
      # Find the rows that match the trial number
      trial_i <- which(results[, 1] == i)
      matching_rows <- which(results[, "Trial"] == i)

      # Pull the true state data
      full_data_name <- paste0(full_data_path, ".", i, ".Full_data.txt")
      trait_data <- read.table(full_data_name, skip = 1, sep = "\t")
      trait_data[,4] <- match(trait_data[,4], c(0,1,10,11))

      # Put the 4-state character into the results matrix
      rTaxon <- results[matching_rows, "rTaxon"]
      for (j in seq_along(rTaxon)) {results[matching_rows[j], "True_4States"] <- trait_data[rTaxon[j], 4]}

      # Call the count data and sort it into a vector
      count_name <- paste0("ConstantRates/", type, "/Data/", type, ".", i, ".Counts.txt")
      counts <- read.table(count_name, skip = 1)
      count_vector <- c(counts[1, 2], counts[2, 2], counts[1, 3], counts[2, 3])

      # Remove the unknown taxa and store it for later
      truestate_vector <- results[matching_rows, "True_4States"]
      for (j in seq_along(rTaxon)) {count_vector[truestate_vector[j]] <- count_vector[truestate_vector[j]] - 1}
      results[trial_i, c("n00", "n01", "n10", "n11")] <- count_vector
    }

    # Finally, write the matrix somewhere it can be accessed later
    results_name <- paste0("Results/ConstantRates/Multiple/", type, ".Multiple.ResultsFull.txt")
    write.table(results, file = results_name, quote = F, sep = "\t", row.names = F, col.names = T)
  }
}


### Check multiple prediction, variable rates
if (variable_rates == TRUE && multiple_prediction == TRUE) {
  for (type in types) {
    # First, pull the unknown matrix
    unknown_labels <- c("Trial", "rTaxon", "Trait_A", "Trait_B", "Terminal_Branch_Length")
    unknown_name <- paste0("VariableRates/", type, "/Multiple/", type, ".Unknown_info.txt")
    unknown_info <- read.table(unknown_name, skip = 1, col.names = unknown_labels)
    row_num <- nrow(unknown_info)
    full_data_path <- paste0("VariableRates/", type, "/Data/", type)

    # Next, create the new section of the results matrix
    matrix_b <- matrix(data=NA, nrow = row_num, ncol = new_col_number)
    colnames(matrix_b) <- new_col_names

    # Combine the unknown information with the new space for results
    results <- cbind(unknown_info, matrix_b)

    # Fill in the counts data that has already by calculated
    for (i in 1:trial_amount) {
      # Find the rows that match the trial number
      trial_i <- which(results[, 1] == i)
      matching_rows <- which(results[, "Trial"] == i)

      # Pull the true state data
      full_data_name <- paste0(full_data_path, ".", i, ".Full_data.txt")
      trait_data <- read.table(full_data_name, skip = 1, sep = "\t")
      trait_data[,4] <- match(trait_data[,4], c(0,1,10,11))

      # Put the 4-state character into the results matrix
      rTaxon <- results[matching_rows, "rTaxon"]
      for (j in seq_along(rTaxon)) {results[matching_rows[j], "True_4States"] <- trait_data[rTaxon[j], 4]}

      # Call the count data and sort it into a vector
      count_name <- paste0("VariableRates/", type, "/Data/", type, ".", i, ".Counts.txt")
      counts <- read.table(count_name, skip = 1)
      count_vector <- c(counts[1, 2], counts[2, 2], counts[1, 3], counts[2, 3])

      # Remove the unknown taxa and store it for later
      truestate_vector <- results[matching_rows, "True_4States"]
      for (j in seq_along(rTaxon)) {count_vector[truestate_vector[j]] <- count_vector[truestate_vector[j]] - 1}
      results[trial_i, c("n00", "n01", "n10", "n11")] <- count_vector
    }

    # Finally, write the matrix somewhere it can be accessed later
    results_name <- paste0("Results/VariableRates/Multiple/", type, ".Multiple.ResultsFull.txt")
    write.table(results, file = results_name, quote = F, sep = "\t", row.names = F, col.names = T)
  }
}



#############################
### Next, check for clade predictions for constant rates
if (clade_prediction == TRUE) {
  for (type in types) {
    # First, pull the unknown matrix
    unknown_labels <- c("Trial", "rTaxon", "Trait_A", "Trait_B", "Terminal_Branch_Length")
    unknown_name <- paste0("ConstantRates/", type, "/Clade/", type, ".Unknown_info.txt")
    unknown_info <- read.table(unknown_name, skip = 1, col.names = unknown_labels)
    row_num <- nrow(unknown_info)
    full_data_path <- paste0("ConstantRates/", type, "/Data/", type)

    # Next, create the new section of the results matrix
    matrix_b <- matrix(data=NA, nrow = row_num, ncol = new_col_number)
    colnames(matrix_b) <- new_col_names

    # Combine the unknown information with the new space for results
    results <- cbind(unknown_info, matrix_b)

    # Fill in the counts data that has already by calculated
    for (i in 1:trial_amount) {
      # Find the rows that match the trial number
      trial_i <- which(results[, 1] == i)
      matching_rows <- which(results[, "Trial"] == i)

      # Pull the true state data
      full_data_name <- paste0(full_data_path, ".", i, ".Full_data.txt")
      trait_data <- read.table(full_data_name, skip = 1, sep = "\t")
      trait_data[,4] <- match(trait_data[,4], c(0,1,10,11))

      # Put the 4-state character into the results matrix
      rTaxon <- results[matching_rows, "rTaxon"]
      for (j in seq_along(rTaxon)) {results[matching_rows[j], "True_4States"] <- trait_data[rTaxon[j], 4]}

      # Call the count data and sort it into a vector
      count_name <- paste0("ConstantRates/", type, "/Data/", type, ".", i, ".Counts.txt")
      counts <- read.table(count_name, skip = 1)
      count_vector <- c(counts[1, 2], counts[2, 2], counts[1, 3], counts[2, 3])

      # Remove the unknown taxa and store it for later
      truestate_vector <- results[matching_rows, "True_4States"]
      for (j in seq_along(rTaxon)) {count_vector[truestate_vector[j]] <- count_vector[truestate_vector[j]] - 1}
      results[trial_i, c("n00", "n01", "n10", "n11")] <- count_vector
    }

    # Finally, write the matrix somewhere it can be accessed later
    results_name <- paste0("Results/ConstantRates/Clade/", type, ".Clade.ResultsFull.txt")
    write.table(results, file = results_name, quote = F, sep = "\t", row.names = F, col.names = T)
  }
}


### Check clade prediction, variable rates
if (variable_rates == TRUE && clade_prediction == TRUE) {
  for (type in types) {
    # First, pull the unknown matrix
    unknown_labels <- c("Trial", "rTaxon", "Trait_A", "Trait_B", "Terminal_Branch_Length")
    unknown_name <- paste0("VariableRates/", type, "/Clade/", type, ".Unknown_info.txt")
    unknown_info <- read.table(unknown_name, skip = 1, col.names = unknown_labels)
    row_num <- nrow(unknown_info)
    full_data_path <- paste0("VariableRates/", type, "/Data/", type)

    # Next, create the new section of the results matrix
    matrix_b <- matrix(data=NA, nrow = row_num, ncol = new_col_number)
    colnames(matrix_b) <- new_col_names

    # Combine the unknown information with the new space for results
    results <- cbind(unknown_info, matrix_b)

    # Fill in the counts data that has already by calculated
    for (i in 1:trial_amount) {
      # Find the rows that match the trial number
      trial_i <- which(results[, 1] == i)
      matching_rows <- which(results[, "Trial"] == i)

      # Pull the true state data
      full_data_name <- paste0(full_data_path, ".", i, ".Full_data.txt")
      trait_data <- read.table(full_data_name, skip = 1, sep = "\t")
      trait_data[,4] <- match(trait_data[,4], c(0,1,10,11))

      # Put the 4-state character into the results matrix
      rTaxon <- results[matching_rows, "rTaxon"]
      for (j in seq_along(rTaxon)) {results[matching_rows[j], "True_4States"] <- trait_data[rTaxon[j], 4]}

      # Call the count data and sort it into a vector
      count_name <- paste0("VariableRates/", type, "/Data/", type, ".", i, ".Counts.txt")
      counts <- read.table(count_name, skip = 1)
      count_vector <- c(counts[1, 2], counts[2, 2], counts[1, 3], counts[2, 3])

      # Remove the unknown taxa and store it for later
      truestate_vector <- results[matching_rows, "True_4States"]
      for (j in seq_along(rTaxon)) {count_vector[truestate_vector[j]] <- count_vector[truestate_vector[j]] - 1}
      results[trial_i, c("n00", "n01", "n10", "n11")] <- count_vector
    }

    # Finally, write the matrix somewhere it can be accessed later
    results_name <- paste0("Results/VariableRates/Clade/", type, ".Clade.ResultsFull.txt")
    write.table(results, file = results_name, quote = F, sep = "\t", row.names = F, col.names = T)
  }
}



#############################
### Next, check if we are running randomly generated multiple prediction
if (multiple_prediction == TRUE) {
  unknown_labels <- c("Trial", "rTaxon", "Trait_A", "Trait_B", "Terminal_Branch_Length")
  unknown_name <- "Random/Multiple/Random.Unknown_info.txt"
  unknown_info <- read.table(unknown_name, skip = 1, col.names = unknown_labels)
  row_num <- nrow(unknown_info)
  full_data_path <- "Random/Data/Random"

  # Next, create the new section of the results matrix
  matrix_b <- matrix(data=NA, nrow = row_num, ncol = new_col_number)
  colnames(matrix_b) <- new_col_names

  # Combine the unknown information with the new space for results
  results <- cbind(unknown_info, matrix_b)

  # Fill in the counts data that has already by calculated
  for (i in 1:trial_amount) {
    # Find the rows that match the trial number
    trial_i <- which(results[, 1] == i)
    matching_rows <- which(results[, "Trial"] == i)

    # Pull the true state data
    full_data_name <- paste0(full_data_path, ".", i, ".Full_data.txt")
    trait_data <- read.table(full_data_name, skip = 1, sep = "\t")
    trait_data[,4] <- match(trait_data[,4], c(0,1,10,11))

    # Put the 4-state character into the results matrix
    rTaxon <- results[matching_rows, "rTaxon"]
    for (j in seq_along(rTaxon)) {results[matching_rows[j], "True_4States"] <- trait_data[rTaxon[j], 4]}

    # Call the count data and sort it into a vector
    count_name <- paste0("Random/Data/Random.", i, ".Counts.txt")
    counts <- read.table(count_name, skip = 1)
    count_vector <- c(counts[1, 2], counts[2, 2], counts[1, 3], counts[2, 3])

    # Remove the unknown taxa and store it for later
    truestate_vector <- results[matching_rows, "True_4States"]
    for (j in seq_along(rTaxon)) {count_vector[truestate_vector[j]] <- count_vector[truestate_vector[j]] - 1}
    results[trial_i, c("n00", "n01", "n10", "n11")] <- count_vector
  }

  # Finally, write the matrix somewhere it can be accessed later
  results_name <- "Results/Random/Random.Multiple.ResultsFull.txt"
  write.table(results, file = results_name, quote = F, sep = "\t", row.names = F, col.names = T)
}



### Finally, check if we are running randomly generated clade prediction
if (clade_prediction == TRUE) {
  unknown_labels <- c("Trial", "rTaxon", "Trait_A", "Trait_B", "Terminal_Branch_Length")
  unknown_name <- "Random/Clade/Random.Unknown_info.txt"
  unknown_info <- read.table(unknown_name, skip = 1, col.names = unknown_labels)
  row_num <- nrow(unknown_info)
  full_data_path <- "Random/Data/Random"

  # Next, create the new section of the results matrix
  matrix_b <- matrix(data=NA, nrow = row_num, ncol = new_col_number)
  colnames(matrix_b) <- new_col_names

  # Combine the unknown information with the new space for results
  results <- cbind(unknown_info, matrix_b)

  # Fill in the counts data that has already by calculated
  for (i in 1:trial_amount) {
    # Find the rows that match the trial number
    trial_i <- which(results[, 1] == i)
    matching_rows <- which(results[, "Trial"] == i)

    # Pull the true state data
    full_data_name <- paste0(full_data_path, ".", i, ".Full_data.txt")
    trait_data <- read.table(full_data_name, skip = 1, sep = "\t")
    trait_data[,4] <- match(trait_data[,4], c(0,1,10,11))

    # Put the 4-state character into the results matrix
    rTaxon <- results[matching_rows, "rTaxon"]
    for (j in seq_along(rTaxon)) {results[matching_rows[j], "True_4States"] <- trait_data[rTaxon[j], 4]}

    # Call the count data and sort it into a vector
    count_name <- paste0("Random/Data/Random.", i, ".Counts.txt")
    counts <- read.table(count_name, skip = 1)
    count_vector <- c(counts[1, 2], counts[2, 2], counts[1, 3], counts[2, 3])

    # Remove the unknown taxa and store it for later
    truestate_vector <- results[matching_rows, "True_4States"]
    for (j in seq_along(rTaxon)) {count_vector[truestate_vector[j]] <- count_vector[truestate_vector[j]] - 1}
    results[trial_i, c("n00", "n01", "n10", "n11")] <- count_vector
  }

  # Finally, write the matrix somewhere it can be accessed later
  results_name <- "Results/Random/Random.Clade.ResultsFull.txt"
  write.table(results, file = results_name, quote = F, sep = "\t", row.names = F, col.names = T)
}

print("Finished generating Results tables. They can be found in the Results folder.")
print("Results tables are currently mostly blank. Following scripts will fill them in.")