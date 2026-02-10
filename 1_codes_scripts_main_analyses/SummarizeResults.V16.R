# This script is for summarizing the results sections


#####
# Start by establishing the necessary parameters
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



#####
# ### To begin, we need to know how big the results table will be
# We start with what we know (rTaxon info will be added later)
tests <- c("Beta_Binom", "Naive_Bayes")
test_size <- 2
# Expand for the BayesTraits trials
if (RJmodel == "MCMC") {
  BTruns <- c("MCMC_Ind", "MCMC_Dep")
  BTsize <- 2
  if (isTRUE(multistate_prediction)) {
    BTruns <- c("MCMC_Multi", BTruns)
    BTsize <- BTsize + 1
  }
}
if (RJmodel == "RJMCMC") {
  BTruns <- c("RJ_Ind", "RJ_Dep")
  BTsize <- 2
  if (isTRUE(multistate_prediction)) {
    BTruns <- c("RJ_Multi", BTruns)
    BTsize <- BTsize + 1
  }
}
if (RJmodel == "BOTH") {
  BTruns <- c("MCMC_Ind", "MCMC_Dep", "RJ_Ind", "RJ_Dep")
  BTsize <- 4
  if (isTRUE(multistate_prediction)) {
    BTruns <- c("MCMC_Multi", "MCMC_Ind", "MCMC_Dep", "RJ_Multi", "RJ_Ind", "RJ_Dep")
    BTsize <- BTsize + 2
  }
}

# Determine the total number of models being tested and their names
tests <- c(tests, BTruns)
test_size <- test_size + BTsize

# Now, we separate out each into three columns
acc_labels <- paste0(tests, "_acc")
LL_labels <- paste0(tests, "_LL")

# Add the medians for each state
n_cols <- c("n00", "n01", "n10", "n11")
n_labels <- paste0("Median_", n_cols)

# We readjust our new column names and column count to match
col_names <- c("Matrix_Type", acc_labels, LL_labels, n_labels)
col_number <- (test_size * 2) + 1 + 4
row_number <- length(types) + 1

# Establish the framework for the final tables
summary_blank <- matrix(data = NA, nrow = row_number, ncol = col_number)
colnames(summary_blank) <- col_names


#####
# Begin summarizing with Constant Rates, Single Prediction
CR.Tip <- summary_blank

# Get the results path
path <- "Results/ConstantRates/Single/"

for (i in 1:length(types)) {
  type <- types[i]
  data_name <- paste0(path, type, ".Single.ResultsFull.txt")
  colnames_data <- strsplit(readLines(data_name, n = 1), "\t")[[1]]
  data <- read.table(data_name, skip = 1, sep = "\t", col.names = colnames_data)

  # Identify columns ending in "_acc" or "_LL"
  target_cols <- grep("(_acc|_LL)$", colnames(data), value = TRUE)

  # Take the mean of each of these columns
  col_means <- colMeans(data[, target_cols, drop = FALSE], na.rm = TRUE)

  # Take the median of the state distributions as well
  CR.Tip[i, n_labels[1]] <- median(data[, n_cols[1]], na.rm = TRUE)
  CR.Tip[i, n_labels[2]] <- median(data[, n_cols[2]], na.rm = TRUE)
  CR.Tip[i, n_labels[3]] <- median(data[, n_cols[3]], na.rm = TRUE)
  CR.Tip[i, n_labels[4]] <- median(data[, n_cols[4]], na.rm = TRUE)

  # Append the column names and data to include the matrix
  CR.Tip[i, 1] <- type
  CR.Tip[i, target_cols] <- col_means
}
# Add in the random results as a control
data_name <- paste0("Results/Random/Random.Single.ResultsFull.txt")
colnames_data <- strsplit(readLines(data_name, n = 1), "\t")[[1]]
data <- read.table(data_name, skip = 1, sep = "\t", col.names = colnames_data)

# Identify columns ending in "_acc" or "_LL"
target_cols <- grep("(_acc|_LL)$", colnames(data), value = TRUE)

# Take the mean of each of these columns
col_means <- colMeans(data[, target_cols, drop = FALSE], na.rm = TRUE)

# Append the column names and data to include the matrix
random_index <- nrow(CR.Tip)
CR.Tip[random_index, 1] <- "Random"
CR.Tip[random_index, target_cols] <- col_means
CR.Tip[random_index, n_labels[1]] <- median(data[, n_cols[1]], na.rm = TRUE)
CR.Tip[random_index, n_labels[2]] <- median(data[, n_cols[2]], na.rm = TRUE)
CR.Tip[random_index, n_labels[3]] <- median(data[, n_cols[3]], na.rm = TRUE)
CR.Tip[random_index, n_labels[4]] <- median(data[, n_cols[4]], na.rm = TRUE)

# Save the table
write.table(CR.Tip, file = "Results/SummarizedResults.CR.Tip.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)



#####
# Now Variable Rates, if run
if (isTRUE(variable_rates)) {
  VR.Tip <- summary_blank

  # Get the results path
  path <- "Results/VariableRates/Single/"

  for (i in 1:length(types)) {
    type <- types[i]
    data_name <- paste0(path, type, ".Single.ResultsFull.txt")
    colnames_data <- strsplit(readLines(data_name, n = 1), "\t")[[1]]
    data <- read.table(data_name, skip = 1, sep = "\t", col.names = colnames_data)

    # Identify columns ending in "_acc" or "_LL"
    target_cols <- grep("(_acc|_LL)$", colnames(data), value = TRUE)

    # Take the mean of each of these columns
    col_means <- colMeans(data[, target_cols, drop = FALSE], na.rm = TRUE)

    # Take the median of the state distributions as well
    VR.Tip[i, n_labels[1]] <- median(data[, n_cols[1]], na.rm = TRUE)
    VR.Tip[i, n_labels[2]] <- median(data[, n_cols[2]], na.rm = TRUE)
    VR.Tip[i, n_labels[3]] <- median(data[, n_cols[3]], na.rm = TRUE)
    VR.Tip[i, n_labels[4]] <- median(data[, n_cols[4]], na.rm = TRUE)

    # Append the column names and data to include the matrix
    VR.Tip[i, 1] <- type
    VR.Tip[i, target_cols] <- col_means
  }
  # Add in the random results as a control
  data_name <- paste0("Results/Random/Random.Single.ResultsFull.txt")
  colnames_data <- strsplit(readLines(data_name, n = 1), "\t")[[1]]
  data <- read.table(data_name, skip = 1, sep = "\t", col.names = colnames_data)

  # Identify columns ending in "_acc" or "_LL"
  target_cols <- grep("(_acc|_LL)$", colnames(data), value = TRUE)

  # Take the mean of each of these columns
  col_means <- colMeans(data[, target_cols, drop = FALSE], na.rm = TRUE)

  # Append the column names and data to include the matrix
  random_index <- nrow(VR.Tip)
  VR.Tip[random_index, 1] <- "Random"
  VR.Tip[random_index, target_cols] <- col_means
  VR.Tip[random_index, n_labels[1]] <- median(data[, n_cols[1]], na.rm = TRUE)
  VR.Tip[random_index, n_labels[2]] <- median(data[, n_cols[2]], na.rm = TRUE)
  VR.Tip[random_index, n_labels[3]] <- median(data[, n_cols[3]], na.rm = TRUE)
  VR.Tip[random_index, n_labels[4]] <- median(data[, n_cols[4]], na.rm = TRUE)

  # Save the table
  write.table(VR.Tip, file = "Results/SummarizedResults.VR.Tip.txt", sep = "\t",
              row.names = F, col.names = T, quote = F)
}



#####
# Now loop back in case of multiple prediction
if (isTRUE(multiple_prediction)) {
  # Begin summarizing with Constant Rates, Multiple Prediction
  CR.Multi <- summary_blank

  # Get the results path
  path <- "Results/ConstantRates/Multiple/"

  for (i in 1:length(types)) {
    type <- types[i]
    data_name <- paste0(path, type, ".Multiple.ResultsFull.txt")
    colnames_data <- strsplit(readLines(data_name, n = 1), "\t")[[1]]
    data <- read.table(data_name, skip = 1, sep = "\t", col.names = colnames_data)

    # Identify columns ending in "_acc" or "_LL"
    target_cols <- grep("(_acc|_LL)$", colnames(data), value = TRUE)

    # Take the mean of each of these columns
    col_means <- colMeans(data[, target_cols, drop = FALSE], na.rm = TRUE)

    # Take the median of the state distributions as well
    CR.Multi[i, n_labels[1]] <- median(data[, n_cols[1]], na.rm = TRUE)
    CR.Multi[i, n_labels[2]] <- median(data[, n_cols[2]], na.rm = TRUE)
    CR.Multi[i, n_labels[3]] <- median(data[, n_cols[3]], na.rm = TRUE)
    CR.Multi[i, n_labels[4]] <- median(data[, n_cols[4]], na.rm = TRUE)

    # Append the column names and data to include the matrix
    CR.Multi[i, 1] <- type
    CR.Multi[i, target_cols] <- col_means
  }
  # Add in the random results as a control
  data_name <- paste0("Results/Random/Random.Multiple.ResultsFull.txt")
  colnames_data <- strsplit(readLines(data_name, n = 1), "\t")[[1]]
  data <- read.table(data_name, skip = 1, sep = "\t", col.names = colnames_data)

  # Identify columns ending in "_acc" or "_LL"
  target_cols <- grep("(_acc|_LL)$", colnames(data), value = TRUE)

  # Take the mean of each of these columns
  col_means <- colMeans(data[, target_cols, drop = FALSE], na.rm = TRUE)

  # Append the column names and data to include the matrix
  random_index <- nrow(CR.Multi)
  CR.Multi[random_index, 1] <- "Random"
  CR.Multi[random_index, target_cols] <- col_means
  CR.Multi[random_index, n_labels[1]] <- median(data[, n_cols[1]], na.rm = TRUE)
  CR.Multi[random_index, n_labels[2]] <- median(data[, n_cols[2]], na.rm = TRUE)
  CR.Multi[random_index, n_labels[3]] <- median(data[, n_cols[3]], na.rm = TRUE)
  CR.Multi[random_index, n_labels[4]] <- median(data[, n_cols[4]], na.rm = TRUE)

  # Save the table
  write.table(CR.Multi, file = "Results/SummarizedResults.CR.Multiple.txt", sep = "\t",
              row.names = F, col.names = T, quote = F)



#####
  # Now Variable Rates, if run
  if (isTRUE(variable_rates)) {
    VR.Multi <- summary_blank

    # Get the results path
    path <- "Results/VariableRates/Multiple/"

    for (i in 1:length(types)) {
      type <- types[i]
      data_name <- paste0(path, type, ".Multiple.ResultsFull.txt")
      colnames_data <- strsplit(readLines(data_name, n = 1), "\t")[[1]]
      data <- read.table(data_name, skip = 1, sep = "\t", col.names = colnames_data)

      # Identify columns ending in "_acc" or "_LL"
      target_cols <- grep("(_acc|_LL)$", colnames(data), value = TRUE)

      # Take the mean of each of these columns
      col_means <- colMeans(data[, target_cols, drop = FALSE], na.rm = TRUE)

      # Take the median of the state distributions as well
      VR.Multi[i, n_labels[1]] <- median(data[, n_cols[1]], na.rm = TRUE)
      VR.Multi[i, n_labels[2]] <- median(data[, n_cols[2]], na.rm = TRUE)
      VR.Multi[i, n_labels[3]] <- median(data[, n_cols[3]], na.rm = TRUE)
      VR.Multi[i, n_labels[4]] <- median(data[, n_cols[4]], na.rm = TRUE)

      # Append the column names and data to include the matrix
      VR.Multi[i, 1] <- type
      VR.Multi[i, target_cols] <- col_means
    }
    # Add in the random results as a control
    data_name <- paste0("Results/Random/Random.Multiple.ResultsFull.txt")
    colnames_data <- strsplit(readLines(data_name, n = 1), "\t")[[1]]
    data <- read.table(data_name, skip = 1, sep = "\t", col.names = colnames_data)

    # Identify columns ending in "_acc" or "_LL"
    target_cols <- grep("(_acc|_LL)$", colnames(data), value = TRUE)

    # Take the mean of each of these columns
    col_means <- colMeans(data[, target_cols, drop = FALSE], na.rm = TRUE)

    # Append the column names and data to include the matrix
    random_index <- nrow(VR.Multi)
    VR.Multi[random_index, 1] <- "Random"
    VR.Multi[random_index, target_cols] <- col_means
    VR.Multi[random_index, n_labels[1]] <- median(data[, n_cols[1]], na.rm = TRUE)
    VR.Multi[random_index, n_labels[2]] <- median(data[, n_cols[2]], na.rm = TRUE)
    VR.Multi[random_index, n_labels[3]] <- median(data[, n_cols[3]], na.rm = TRUE)
    VR.Multi[random_index, n_labels[4]] <- median(data[, n_cols[4]], na.rm = TRUE)

    # Save the table
    write.table(VR.Multi, file = "Results/SummarizedResults.VR.Multiple.txt", sep = "\t",
                row.names = F, col.names = T, quote = F)
  }
}



#####
# Now loop back in case of clade prediction
if (isTRUE(clade_prediction)) {
  # Begin summarizing with Constant Rates
  CR.Clade <- summary_blank

  # Get the results path
  path <- "Results/ConstantRates/Clade/"

  for (i in 1:length(types)) {
    type <- types[i]
    data_name <- paste0(path, type, ".Clade.ResultsFull.txt")
    colnames_data <- strsplit(readLines(data_name, n = 1), "\t")[[1]]
    data <- read.table(data_name, skip = 1, sep = "\t", col.names = colnames_data)

    # Identify columns ending in "_acc" or "_LL"
    target_cols <- grep("(_acc|_LL)$", colnames(data), value = TRUE)

    # Take the mean of each of these columns
    col_means <- colMeans(data[, target_cols, drop = FALSE], na.rm = TRUE)

    # Take the median of the state distributions as well
    CR.Clade[i, n_labels[1]] <- median(data[, n_cols[1]], na.rm = TRUE)
    CR.Clade[i, n_labels[2]] <- median(data[, n_cols[2]], na.rm = TRUE)
    CR.Clade[i, n_labels[3]] <- median(data[, n_cols[3]], na.rm = TRUE)
    CR.Clade[i, n_labels[4]] <- median(data[, n_cols[4]], na.rm = TRUE)

    # Append the column names and data to include the matrix
    CR.Clade[i, 1] <- type
    CR.Clade[i, target_cols] <- col_means
  }
  # Add in the random results as a control
  data_name <- paste0("Results/Random/Random.Clade.ResultsFull.txt")
  colnames_data <- strsplit(readLines(data_name, n = 1), "\t")[[1]]
  data <- read.table(data_name, skip = 1, sep = "\t", col.names = colnames_data)

  # Identify columns ending in "_acc" or "_LL"
  target_cols <- grep("(_acc|_LL)$", colnames(data), value = TRUE)

  # Take the mean of each of these columns
  col_means <- colMeans(data[, target_cols, drop = FALSE], na.rm = TRUE)

  # Append the column names and data to include the matrix
  random_index <- nrow(CR.Clade)
  CR.Clade[random_index, 1] <- "Random"
  CR.Clade[random_index, target_cols] <- col_means
  CR.Clade[random_index, n_labels[1]] <- median(data[, n_cols[1]], na.rm = TRUE)
  CR.Clade[random_index, n_labels[2]] <- median(data[, n_cols[2]], na.rm = TRUE)
  CR.Clade[random_index, n_labels[3]] <- median(data[, n_cols[3]], na.rm = TRUE)
  CR.Clade[random_index, n_labels[4]] <- median(data[, n_cols[4]], na.rm = TRUE)

  # Save the table
  write.table(CR.Clade, file = "Results/SummarizedResults.CR.Clade.txt", sep = "\t",
              row.names = F, col.names = T, quote = F)



  #####
  # Now Variable Rates, if run
  if (isTRUE(variable_rates)) {
    VR.Clade <- summary_blank

    # Get the results path
    path <- "Results/VariableRates/Clade/"

    for (i in 1:length(types)) {
      type <- types[i]
      data_name <- paste0(path, type, ".Clade.ResultsFull.txt")
      colnames_data <- strsplit(readLines(data_name, n = 1), "\t")[[1]]
      data <- read.table(data_name, skip = 1, sep = "\t", col.names = colnames_data)

      # Identify columns ending in "_acc" or "_LL"
      target_cols <- grep("(_acc|_LL)$", colnames(data), value = TRUE)

      # Take the mean of each of these columns
      col_means <- colMeans(data[, target_cols, drop = FALSE], na.rm = TRUE)

      # Take the median of the state distributions as well
      VR.Clade[i, n_labels[1]] <- median(data[, n_cols[1]], na.rm = TRUE)
      VR.Clade[i, n_labels[2]] <- median(data[, n_cols[2]], na.rm = TRUE)
      VR.Clade[i, n_labels[3]] <- median(data[, n_cols[3]], na.rm = TRUE)
      VR.Clade[i, n_labels[4]] <- median(data[, n_cols[4]], na.rm = TRUE)

      # Append the column names and data to include the matrix
      VR.Clade[i, 1] <- type
      VR.Clade[i, target_cols] <- col_means
    }
    # Add in the random results as a control
    data_name <- paste0("Results/Random/Random.Clade.ResultsFull.txt")
    colnames_data <- strsplit(readLines(data_name, n = 1), "\t")[[1]]
    data <- read.table(data_name, skip = 1, sep = "\t", col.names = colnames_data)

    # Identify columns ending in "_acc" or "_LL"
    target_cols <- grep("(_acc|_LL)$", colnames(data), value = TRUE)

    # Take the mean of each of these columns
    col_means <- colMeans(data[, target_cols, drop = FALSE], na.rm = TRUE)

    # Append the column names and data to include the matrix
    random_index <- nrow(VR.Clade)
    VR.Clade[random_index, 1] <- "Random"
    VR.Clade[random_index, target_cols] <- col_means
    VR.Clade[random_index, n_labels[1]] <- median(data[, n_cols[1]], na.rm = TRUE)
    VR.Clade[random_index, n_labels[2]] <- median(data[, n_cols[2]], na.rm = TRUE)
    VR.Clade[random_index, n_labels[3]] <- median(data[, n_cols[3]], na.rm = TRUE)
    VR.Clade[random_index, n_labels[4]] <- median(data[, n_cols[4]], na.rm = TRUE)

    # Save the table
    write.table(VR.Clade, file = "Results/SummarizedResults.VR.Clade.txt", sep = "\t",
                row.names = F, col.names = T, quote = F)
  }
}
print("Finished Summarizing all available results.")
