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





### Install and call the necessary
# First, a list of the packages we need
required_packages <- c("ape", "phytools", "castor")

# Install any packages that aren't already installed
installed <- required_packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(required_packages[!installed])
}

# Load the libraries
lapply(required_packages, library, character.only = TRUE)





### For loop to record all of this information for Constant Rates
for (type in types) {
  # Print a progress check
  print(paste0("Beginning the ancestral reconstruction for ", type, " with Constant Rates."))

  # Call the necessary paths
  data_path <- paste0("ConstantRates/", type, "/Data/", type)
  results_name <- paste0("Results/ConstantRates/Single/", type, ".Single.ResultsFull.txt")

  # Call the results table
  colnames_results <- strsplit(readLines(results_name, n = 1), "\t")[[1]]
  results <- read.table(results_name, skip = 1, sep = "\t", col.names = colnames_results)

  # Create a matrix to store the information based on tests run
  matrix_name <- paste0("Results/ConstantRates/Single/AncestralStateReconstruction/", type, ".AncReconAndTreeHeight.txt")

  if (RJmodel == "MCMC") {
    AncConMat <- matrix(data=NA, nrow = nrow(results), ncol = 42)
    colnames(AncConMat) <- c("Trial_#", "Unknown_Taxon", "True_A", "True_B", "True_4States", "Avg_Sis_A", "Avg_Sis_B",
                             "Sis_IndMCMC_00", "Sis_IndMCMC_01", "Sis_IndMCMC_10", "Sis_IndMCMC_11",
                             "Sis_DepMCMC_00", "Sis_DepMCMC_01", "Sis_DepMCMC_10", "Sis_DepMCMC_11",
                             "Sis_IndMCMC_Pred", "Sis_IndMCMC_Lh", "Sis_DepMCMC_Pred", "Sis_DepMCMC_Lh",
                             "Anc_IndMCMC_00", "Anc_IndMCMC_01", "Anc_IndMCMC_10", "Anc_IndMCMC_11",
                             "Anc_DepMCMC_00", "Anc_DepMCMC_01", "Anc_DepMCMC_10", "Anc_DepMCMC_11",
                             "Anc_IndMCMC_Pred", "Anc_IndMCMC_Lh", "Anc_DepMCMC_Pred", "Anc_DepMCMC_Lh",
                             "True_IsMax_Sis_IndMCMC", "Max_Sis_IndMCMC_Lh", "True_IsMax_Sis_DepMCMC", "Max_Sis_DepMCMC_Lh",
                             "True_IsMax_Anc_IndMCMC", "Max_Anc_IndMCMC_Lh", "True_IsMax_Anc_DepMCMC", "Max_Anc_DepMCMC_Lh",
                             "TBL", "Tree_Height", "Std_TBL")
  } else if (RJmodel == "RJMCMC") {
    AncConMat <- matrix(data=NA, nrow = nrow(results), ncol = 42)
    colnames(AncConMat) <- c("Trial_#", "Unknown_Taxon", "True_A", "True_B", "True_4States", "Avg_Sis_A", "Avg_Sis_B",
                             "Sis_IndRJ_00", "Sis_IndRJ_01", "Sis_IndRJ_10", "Sis_IndRJ_11",
                             "Sis_DepRJ_00", "Sis_DepRJ_01", "Sis_DepRJ_10", "Sis_DepRJ_11",
                             "Sis_IndRJ_Pred", "Sis_IndRJ_Lh", "Sis_DepRJ_Pred", "Sis_DepRJ_Lh",
                             "Anc_IndRJ_00", "Anc_IndRJ_01", "Anc_IndRJ_10", "Anc_IndRJ_11",
                             "Anc_DepRJ_00", "Anc_DepRJ_01", "Anc_DepRJ_10", "Anc_DepRJ_11",
                             "Anc_IndRJ_Pred", "Anc_IndRJ_Lh", "Anc_DepRJ_Pred", "Anc_DepRJ_Lh",
                             "True_IsMax_Sis_IndRJ", "Max_Sis_IndRJ_Lh", "True_IsMax_Sis_DepRJ", "Max_Sis_DepRJ_Lh",
                             "True_IsMax_Anc_IndRJ", "Max_Anc_IndRJ_Lh", "True_IsMax_Anc_DepRJ", "Max_Anc_DepRJ_Lh",
                             "TBL", "Tree_Height", "Std_TBL")
  } else if (RJmodel == "BOTH") {
    AncConMat <- matrix(data=NA, nrow = nrow(results), ncol = 74)
    colnames(AncConMat) <- c("Trial_#", "Unknown_Taxon", "True_A", "True_B", "True_4States", "Avg_Sis_A", "Avg_Sis_B",
                             "Sis_IndMCMC_00", "Sis_IndMCMC_01", "Sis_IndMCMC_10", "Sis_IndMCMC_11",
                             "Sis_DepMCMC_00", "Sis_DepMCMC_01", "Sis_DepMCMC_10", "Sis_DepMCMC_11",
                             "Sis_IndRJ_00", "Sis_IndRJ_01", "Sis_IndRJ_10", "Sis_IndRJ_11",
                             "Sis_DepRJ_00", "Sis_DepRJ_01", "Sis_DepRJ_10", "Sis_DepRJ_11",
                             "Sis_IndMCMC_Pred", "Sis_IndMCMC_Lh", "Sis_DepMCMC_Pred", "Sis_DepMCMC_Lh",
                             "Sis_IndRJ_Pred", "Sis_IndRJ_Lh", "Sis_DepRJ_Pred", "Sis_DepRJ_Lh",
                             "Anc_IndMCMC_00", "Anc_IndMCMC_01", "Anc_IndMCMC_10", "Anc_IndMCMC_11",
                             "Anc_DepMCMC_00", "Anc_DepMCMC_01", "Anc_DepMCMC_10", "Anc_DepMCMC_11",
                             "Anc_IndRJ_00", "Anc_IndRJ_01", "Anc_IndRJ_10", "Anc_IndRJ_11",
                             "Anc_DepRJ_00", "Anc_DepRJ_01", "Anc_DepRJ_10", "Anc_DepRJ_11",
                             "Anc_IndMCMC_Pred", "Anc_IndMCMC_Lh", "Anc_DepMCMC_Pred", "Anc_DepMCMC_Lh",
                             "Anc_IndRJ_Pred", "Anc_IndRJ_Lh", "Anc_DepRJ_Pred", "Anc_DepRJ_Lh",
                             "True_IsMax_Sis_IndMCMC", "Max_Sis_IndMCMC_Lh", "True_IsMax_Sis_DepMCMC", "Max_Sis_DepMCMC_Lh",
                             "True_IsMax_Sis_IndRJ", "Max_Sis_IndRJ_Lh", "True_IsMax_Sis_DepRJ", "Max_Sis_DepRJ_Lh",
                             "True_IsMax_Anc_IndMCMC", "Max_Anc_IndMCMC_Lh", "True_IsMax_Anc_DepMCMC", "Max_Anc_DepMCMC_Lh",
                             "True_IsMax_Anc_IndRJ", "Max_Anc_IndRJ_Lh", "True_IsMax_Anc_DepRJ", "Max_Anc_DepRJ_Lh",
                             "TBL", "Tree_Height", "Std_TBL")
  }



  # Start the for-loop that will run each prediction and fill it in
  for (i in 1:trial_amount) {
    # First, call rTaxon by subsetting the rows where column 1 matches i
    matching_rows <- which(results[, "Trial"] == i)
    rTaxon <- results[matching_rows, "rTaxon"]

    # Call the data file for this iteration
    data_name <- paste0(data_path, ".", i, ".Full_data.txt")
    data <- read.table(data_name, skip = 1, sep = "\t")
    data[, 4] <- match(data[,4], c(0,1,10,11))


    # Pull the tree to read the terminal branch length
    treename <- paste0("Trees/Full_tree.", i, ".tre")
    tree_complete <- read.nexus(treename)

    # Loop through every taxon for which we are predicting
    for (j in 1:length(rTaxon)) {
      # This collects the numerical average value of the sister taxon's trait information
      AncConMat[matching_rows[j], "Avg_Sis_A"] <- sisterData(tree_complete, rTaxon[j], data, "A")
      AncConMat[matching_rows[j], "Avg_Sis_B"] <- sisterData(tree_complete, rTaxon[j], data, "B")

      # Reset the sister node number and find the other important nodes
      sisters <- getSisters(tree_complete, rTaxon[j], mode = "number")
      sister_tips <- getAllDescendantTips(tree_complete, sisters)
      grandparent <- getParent(tree_complete, getMRCA(tree_complete, c(rTaxon, sister_tips)))
      cousins <- getAllDescendantTips(tree_complete, grandparent)
      cousins <- cousins[cousins != rTaxon[j]]

      # Find the node of interest in a tree where you don't know the unknown's information
      edited_tree <- drop.tip(tree_complete, rTaxon[j])
      edited_data <- data[-rTaxon,]
      edited_cousins <- match(cousins, edited_tree$tip.label)
      grandparent_node_asr <- getMRCA(edited_tree, edited_cousins)



      ### Ancestral State Reconstructin with both traits
      if (as.numeric(sisters) > Ntip(tree_complete) && length(unique(data[, 4])) > 1) {
        # Adjust node number to row number of the output
        sister_asr <- getMRCA(edited_tree, sister_tips)
        sister_asr <- sister_asr - Ntip(edited_tree)

        ### First pull the Independent MCMC model, and perform ML Ancestral State Reconstruction on the nearest sister taxon
        if (RJmodel == "MCMC" || RJmodel == "BOTH") {
          IndMCMC_qmat <- estimated_qmat(type, i, IndependentRates = TRUE, ConstantRates = TRUE, RJ = FALSE, Predict = FALSE)
          ace_IndMCMC <- asr_mk_model(tree = edited_tree, tip_states = edited_data[, 4], Nstates = 4, transition_matrix = IndMCMC_qmat, reroot = FALSE, root_prior = "flat")
          # Pull the likelihood
          Sis_IndMCMC_Lhs <- round(ace_IndMCMC$ancestral_likelihoods[sister_asr, ], digits = 6)
          # Store them in the matrix
          AncConMat[matching_rows[j], c("Sis_IndMCMC_00", "Sis_IndMCMC_01", "Sis_IndMCMC_10", "Sis_IndMCMC_11")] <- Sis_IndMCMC_Lhs

        ### Next Dependent MCMC models
          DepMCMC_qmat <- estimated_qmat(type, i, IndependentRates = FALSE, ConstantRates = T, RJ = FALSE, Predict = FALSE)
          ace_DepMCMC <- asr_mk_model(tree = edited_tree, tip_states = edited_data[, 4], Nstates = 4, transition_matrix = DepMCMC_qmat, reroot = FALSE, root_prior = "flat")
          # Pull the likelihood
          Sis_DepMCMC_Lhs <- round(ace_DepMCMC$ancestral_likelihoods[sister_asr, ], digits = 6)
          # Store them in the matrix
          AncConMat[matching_rows[j], c("Sis_DepMCMC_00", "Sis_DepMCMC_01", "Sis_DepMCMC_10", "Sis_DepMCMC_11")] <- Sis_DepMCMC_Lhs
        }


        ### Next Ind RJ models
        if (RJmodel == "RJMCMC" || RJmodel == "BOTH") {
          IndRJ_qmat <- estimated_qmat(type, i, IndependentRates = TRUE, ConstantRates = T, RJ = T, Predict = FALSE)
          ace_IndRJ <- asr_mk_model(tree = edited_tree, tip_states = edited_data[, 4], Nstates = 4, transition_matrix = IndRJ_qmat, reroot = FALSE, root_prior = "flat")
          # Pull the likelihood
          Sis_IndRJ_Lhs <- round(ace_IndRJ$ancestral_likelihoods[sister_asr, ], digits = 6)
          # Store them in the matrix
          AncConMat[matching_rows[j], c("Sis_IndRJ_00", "Sis_IndRJ_01", "Sis_IndRJ_10", "Sis_IndRJ_11")] <- Sis_IndRJ_Lhs

        ### Finally, the Dep RJ models
          DepRJ_qmat <- estimated_qmat(type, i, IndependentRates = TRUE, ConstantRates = T, RJ = FALSE, Predict = FALSE)
          ace_DepRJ <- asr_mk_model(tree = edited_tree, tip_states = edited_data[, 4], Nstates = 4, transition_matrix = DepRJ_qmat, reroot = FALSE, root_prior = "flat")
          # Pull the likelihood
          Sis_DepRJ_Lhs <- round(ace_DepRJ$ancestral_likelihoods[sister_asr, ], digits = 6)
          # Store them in the matrix
          AncConMat[matching_rows[j], c("Sis_DepRJ_00", "Sis_DepRJ_01", "Sis_DepRJ_10", "Sis_DepRJ_11")] <- Sis_DepRJ_Lhs
        }


        ### Line accounting for non-variation in the tree or single tip prediction
      } else {

        # In case of non-variation in the tree, make this adjustment for sister-node reconstruction
        # This only impacts cases of no trait variation across the tree when the sister node is internal
        # Although wonky, it covers edge cases without impacting sister-leaf reconstruction
        if (as.numeric(sisters) > Ntip(tree_complete)) {sisters <- sisters - Ntip(tree_complete)}

        # First cover the MCMC models
        if (RJmodel == "MCMC" || RJmodel == "BOTH") {
          Sis_IndMCMC_Lhs <- NULL
          Sis_IndMCMC_Lhs[1] <- Sis_IndMCMC_Lhs[2] <- Sis_IndMCMC_Lhs[3] <- Sis_IndMCMC_Lhs[4] <- 0
          Sis_IndMCMC_Lhs[ data[sisters,4] ] <- 1
          AncConMat[matching_rows[j], c("Sis_IndMCMC_00", "Sis_IndMCMC_01", "Sis_IndMCMC_10", "Sis_IndMCMC_11")] <- Sis_IndMCMC_Lhs

          Sis_DepMCMC_Lhs <- NULL
          Sis_DepMCMC_Lhs[1] <- Sis_DepMCMC_Lhs[2] <- Sis_DepMCMC_Lhs[3] <- Sis_DepMCMC_Lhs[4] <- 0
          Sis_DepMCMC_Lhs[ data[sisters,4] ] <- 1
          AncConMat[matching_rows[j], c("Sis_DepMCMC_00", "Sis_DepMCMC_01", "Sis_DepMCMC_10", "Sis_DepMCMC_11")] <- Sis_DepMCMC_Lhs
        }

        # Next cover the RJMCMC models
        if (RJmodel == "RJMCMC" || RJmodel == "BOTH") {
          Sis_IndRJ_Lhs <- NULL
          Sis_IndRJ_Lhs[1] <- Sis_IndRJ_Lhs[2] <- Sis_IndRJ_Lhs[3] <- Sis_IndRJ_Lhs[4] <- 0
          Sis_IndRJ_Lhs[ data[sisters,4] ] <- 1
          AncConMat[matching_rows[j], c("Sis_IndRJ_00", "Sis_IndRJ_01", "Sis_IndRJ_10", "Sis_IndRJ_11")] <- Sis_IndRJ_Lhs

          Sis_DepRJ_Lhs <- NULL
          Sis_DepRJ_Lhs[1] <- Sis_DepRJ_Lhs[2] <- Sis_DepRJ_Lhs[3] <- Sis_DepRJ_Lhs[4] <- 0
          Sis_DepRJ_Lhs[ data[sisters,4] ] <- 1
          AncConMat[matching_rows[j], c("Sis_DepRJ_00", "Sis_DepRJ_01", "Sis_DepRJ_10", "Sis_DepRJ_11")] <- Sis_DepRJ_Lhs
        }
      }


      ### Now pull the information for the ancestral node (the grandpa)
      # First, adjust to match the internal node number
      if (!is.null(grandparent)) {
        grandparent_node_asr <- grandparent_node_asr - Ntip(edited_tree)
      }


      ### Perform ML Ancestral State Reconstruction on the parent node of the Phylogenetic Bracket
      # First the MCMC models
      if ((RJmodel == "MCMC" || RJmodel == "BOTH") && !is.null(grandparent)) {
        IndMCMC_qmat <- estimated_qmat(type, i, IndependentRates = TRUE, ConstantRates = TRUE, RJ = FALSE, Predict = FALSE)
        Anc_ace_IndMCMC <- asr_mk_model(tree = edited_tree, tip_states = edited_data[, 4], Nstates = 4, transition_matrix = IndMCMC_qmat, reroot = FALSE, root_prior = "flat")
        # Pull the likelihood
        Anc_IndMCMC_Lhs <- round(Anc_ace_IndMCMC$ancestral_likelihoods[grandparent_node_asr, ], digits = 6)
        # Store them in the matrix
        AncConMat[matching_rows[j], c("Anc_IndMCMC_00", "Anc_IndMCMC_01", "Anc_IndMCMC_10", "Anc_IndMCMC_11")] <- Anc_IndMCMC_Lhs

        ### Next Dependent MCMC models
        DepMCMC_qmat <- estimated_qmat(type, i, IndependentRates = FALSE, ConstantRates = T, RJ = FALSE, Predict = FALSE)
        Anc_ace_DepMCMC <- asr_mk_model(tree = edited_tree, tip_states = edited_data[, 4], Nstates = 4, transition_matrix = DepMCMC_qmat, reroot = FALSE, root_prior = "flat")
        # Pull the likelihood
        Anc_DepMCMC_Lhs <- round(Anc_ace_DepMCMC$ancestral_likelihoods[grandparent_node_asr, ], digits = 6)
        # Store them in the matrix
        AncConMat[matching_rows[j], c("Anc_DepMCMC_00", "Anc_DepMCMC_01", "Anc_DepMCMC_10", "Anc_DepMCMC_11")] <- Anc_DepMCMC_Lhs
      }

      ### Next RJ models, starting with the Ind model
      if ((RJmodel == "RJMCMC" || RJmodel == "BOTH") && !is.null(grandparent)) {
        IndRJ_qmat <- estimated_qmat(type, i, IndependentRates = TRUE, ConstantRates = T, RJ = T, Predict = FALSE)
        Anc_ace_IndRJ <- asr_mk_model(tree = edited_tree, tip_states = edited_data[, 4], Nstates = 4, transition_matrix = IndRJ_qmat, reroot = FALSE, root_prior = "flat")
        # Pull the likelihood
        Anc_IndRJ_Lhs <- round(Anc_ace_IndRJ$ancestral_likelihoods[grandparent_node_asr, ], digits = 6)
        # Store them in the matrix
        AncConMat[matching_rows[j], c("Anc_IndRJ_00", "Anc_IndRJ_01", "Anc_IndRJ_10", "Anc_IndRJ_11")] <- Anc_IndRJ_Lhs

        ### Finally, the Dep RJ models
        DepRJ_qmat <- estimated_qmat(type, i, IndependentRates = F, ConstantRates = T, RJ = T, Predict = FALSE)
        Anc_ace_DepRJ <- asr_mk_model(tree = edited_tree, tip_states = edited_data[, 4], Nstates = 4, transition_matrix = DepRJ_qmat, reroot = FALSE, root_prior = "flat")
        # Pull the likelihood
        Anc_DepRJ_Lhs <- round(Anc_ace_DepRJ$ancestral_likelihoods[grandparent_node_asr, ], digits = 6)
        # Store them in the matrix
        AncConMat[matching_rows[j], c("Anc_DepRJ_00", "Anc_DepRJ_01", "Anc_DepRJ_10", "Anc_DepRJ_11")] <- Anc_DepRJ_Lhs
      }


      ### Adjust the code if the parent node is the root
      # First the MCMC models
      if ((RJmodel == "MCMC" || RJmodel == "BOTH") && is.null(grandparent)) {
        # Reset the likelihoods to clarify this trial does not have a valid "grandparent"
        Anc_IndMCMC_Lhs <- c(NA, NA, NA, NA)
        # Store them in the matrix
        AncConMat[matching_rows[j], c("Anc_IndMCMC_00", "Anc_IndMCMC_01", "Anc_IndMCMC_10", "Anc_IndMCMC_11")] <- Anc_IndMCMC_Lhs

        ### Next Dependent MCMC models
        # Reset the likelihoods to clarify this trial does not have a valid "grandparent"
        Anc_DepMCMC_Lhs <- c(NA, NA, NA, NA)
        # Store them in the matrix
        AncConMat[matching_rows[j], c("Anc_DepMCMC_00", "Anc_DepMCMC_01", "Anc_DepMCMC_10", "Anc_DepMCMC_11")] <- Anc_DepMCMC_Lhs
      }

      ### Next RJ models, starting with the Ind model
      if ((RJmodel == "RJMCMC" || RJmodel == "BOTH") && is.null(grandparent)) {
        # Reset the likelihoods to clarify this trial does not have a valid "grandparent"
        Anc_IndRJ_Lhs <- c(NA, NA, NA, NA)
        # Store them in the matrix
        AncConMat[matching_rows[j], c("Anc_IndRJ_00", "Anc_IndRJ_01", "Anc_IndRJ_10", "Anc_IndRJ_11")] <- Anc_IndRJ_Lhs

        ### Finally, the Dep RJ models
        # Reset the likelihoods to clarify this trial does not have a valid "grandparent"
        Anc_DepRJ_Lhs <- c(NA, NA, NA, NA)
        # Store them in the matrix
        AncConMat[matching_rows[j], c("Anc_DepRJ_00", "Anc_DepRJ_01", "Anc_DepRJ_10", "Anc_DepRJ_11")] <- Anc_DepRJ_Lhs
      }


      ### Pull the predicted states and likelihoods for the tests performed
      if ((RJmodel == "MCMC" || RJmodel == "BOTH") && !is.null(grandparent)) {
        # First for the sister
        AncConMat[matching_rows[j], "Sis_IndMCMC_Pred"] <- which.max(Sis_IndMCMC_Lhs)
        AncConMat[matching_rows[j], "Sis_IndMCMC_Lh"] <- max(Sis_IndMCMC_Lhs)


        AncConMat[matching_rows[j], "Sis_DepMCMC_Pred"] <- which.max(Sis_DepMCMC_Lhs)
        AncConMat[matching_rows[j], "Sis_DepMCMC_Lh"] <- max(Sis_DepMCMC_Lhs)

        # Now for the grandparent
        AncConMat[matching_rows[j], "Anc_IndMCMC_Pred"] <- which.max(Anc_IndMCMC_Lhs)
        AncConMat[matching_rows[j], "Anc_IndMCMC_Lh"] <- max(Anc_IndMCMC_Lhs)

        AncConMat[matching_rows[j], "Anc_DepMCMC_Pred"] <- which.max(Anc_DepMCMC_Lhs)
        AncConMat[matching_rows[j], "Anc_DepMCMC_Lh"] <- max(Anc_DepMCMC_Lhs)

        # Separate out the likelihood of the "true" state of rTaxon
        results[matching_rows[j], "TrueState_Sister_IndMCMC_Lh"] <- Sis_IndMCMC_Lhs[results[matching_rows[j], "True_4States"]]
        results[matching_rows[j], "TrueState_Sister_DepMCMC_Lh"] <- Sis_DepMCMC_Lhs[results[matching_rows[j], "True_4States"]]

        results[matching_rows[j], "TrueState_Ancestor_IndMCMC_Lh"] <- Anc_IndMCMC_Lhs[results[matching_rows[j], "True_4States"]]
        results[matching_rows[j], "TrueState_Ancestor_DepMCMC_Lh"] <- Anc_DepMCMC_Lhs[results[matching_rows[j], "True_4States"]]
      }

      if ((RJmodel == "RJMCMC" || RJmodel == "BOTH") && !is.null(grandparent)) {
        # First for the sister
        AncConMat[matching_rows[j], "Sis_IndRJ_Pred"] <- which.max(Sis_IndRJ_Lhs)
        AncConMat[matching_rows[j], "Sis_IndRJ_Lh"] <- max(Sis_IndRJ_Lhs)

        AncConMat[matching_rows[j], "Sis_DepRJ_Pred"] <- which.max(Sis_DepRJ_Lhs)
        AncConMat[matching_rows[j], "Sis_DepRJ_Lh"] <- max(Sis_DepRJ_Lhs)

        # Now for the grandparent
        AncConMat[matching_rows[j], "Anc_IndRJ_Pred"] <- which.max(Anc_IndRJ_Lhs)
        AncConMat[matching_rows[j], "Anc_IndRJ_Lh"] <- max(Anc_IndRJ_Lhs)

        AncConMat[matching_rows[j], "Anc_DepRJ_Pred"] <- which.max(Anc_DepRJ_Lhs)
        AncConMat[matching_rows[j], "Anc_DepRJ_Lh"] <- max(Anc_DepRJ_Lhs)

        # Separate out the likelihood of the "true" state of rTaxon
        results[matching_rows[j], "TrueState_Sister_IndRJ_Lh"] <- Sis_IndMCMC_Lhs[results[matching_rows[j], "True_4States"]]
        results[matching_rows[j], "TrueState_Sister_DepRJ_Lh"] <- Sis_DepMCMC_Lhs[results[matching_rows[j], "True_4States"]]

        results[matching_rows[j], "TrueState_Ancestor_IndRJ_Lh"] <- Anc_IndMCMC_Lhs[results[matching_rows[j], "True_4States"]]
        results[matching_rows[j], "TrueState_Ancestor_DepRJ_Lh"] <- Anc_DepMCMC_Lhs[results[matching_rows[j], "True_4States"]]
      }


      ### A fix for when there is no grandparent node (parent is the root)
      # These are again separated by tests run
      if ((RJmodel == "MCMC" || RJmodel == "BOTH") && is.null(grandparent)) {
        # Pull the predicted states and likelihoods for the sister
        AncConMat[matching_rows[j], "Sis_IndMCMC_Pred"] <- which.max(Sis_IndMCMC_Lhs)
        AncConMat[matching_rows[j], "Sis_IndMCMC_Lh"] <- max(Sis_IndMCMC_Lhs)

        AncConMat[matching_rows[j], "Sis_DepMCMC_Pred"] <- which.max(Sis_DepMCMC_Lhs)
        AncConMat[matching_rows[j], "Sis_DepMCMC_Lh"] <- max(Sis_DepMCMC_Lhs)

        # Separate out the likelihood of the "true" state of rTaxon
        results[matching_rows[j], "TrueState_Sister_IndMCMC_Lh"] <- Sis_IndMCMC_Lhs[results[matching_rows[j], "True_4States"]]
        results[matching_rows[j], "TrueState_Sister_DepMCMC_Lh"] <- Sis_DepMCMC_Lhs[results[matching_rows[j], "True_4States"]]

        # However, none of the "ancestor" state reconstruction can be done without a known ancestor.
        AncConMat[matching_rows[j], "Anc_IndMCMC_Pred"] <- NA
        AncConMat[matching_rows[j], "Anc_IndMCMC_Lh"] <- NA

        AncConMat[matching_rows[j], "Anc_DepMCMC_Pred"] <- NA
        AncConMat[matching_rows[j], "Anc_DepMCMC_Lh"] <- NA

        results[matching_rows[j], "TrueState_Ancestor_IndMCMC_Lh"] <- NA
        results[matching_rows[j], "TrueState_Ancestor_DepMCMC_Lh"] <- NA
      }

      if ((RJmodel == "RJMCMC" || RJmodel == "BOTH") && is.null(grandparent)) {
        # Pull the predicted states and likelihoods for the sister
        AncConMat[matching_rows[j], "Sis_IndRJ_Pred"] <- which.max(Sis_IndRJ_Lhs)
        AncConMat[matching_rows[j], "Sis_IndRJ_Lh"] <- max(Sis_IndRJ_Lhs)

        AncConMat[matching_rows[j], "Sis_DepRJ_Pred"] <- which.max(Sis_DepRJ_Lhs)
        AncConMat[matching_rows[j], "Sis_DepRJ_Lh"] <- max(Sis_DepRJ_Lhs)

        # Separate out the likelihood of the "true" state of rTaxon
        results[matching_rows[j], "TrueState_Sister_IndRJ_Lh"] <- Sis_IndRJ_Lhs[results[matching_rows[j], "True_4States"]]
        results[matching_rows[j], "TrueState_Sister_DepRJ_Lh"] <- Sis_DepRJ_Lhs[results[matching_rows[j], "True_4States"]]

        # However, none of the "ancestor" state reconstruction can be done without a known ancestor.
        AncConMat[matching_rows[j], "Anc_IndRJ_Pred"] <- NA
        AncConMat[matching_rows[j], "Anc_IndRJ_Lh"] <- NA

        AncConMat[matching_rows[j], "Anc_DepRJ_Pred"] <- NA
        AncConMat[matching_rows[j], "Anc_DepRJ_Lh"] <- NA

        results[matching_rows[j], "TrueState_Ancestor_IndRJ_Lh"] <- NA
        results[matching_rows[j], "TrueState_Ancestor_DepRJ_Lh"] <- NA
      }


      ### These lines add to the results whether or not the ACR would agree with our predictions or the true, unknown state
      # First the MCMC models
      if (RJmodel == "MCMC" || RJmodel == "BOTH") {
        # First determine if the maximum likelihood from the ASR matches the unknown taxon's character state
        AncConMat[matching_rows[j], "True_IsMax_Sis_IndMCMC"] <- ifelse(round(Sis_IndMCMC_Lhs[results[matching_rows[j], "True_4States"]], digits = 6) == max(round(Sis_IndMCMC_Lhs, digits = 6)), 1, 0)
        AncConMat[matching_rows[j], "True_IsMax_Sis_DepMCMC"] <- ifelse(round(Sis_DepMCMC_Lhs[results[matching_rows[j], "True_4States"]], digits = 6) == max(round(Sis_DepMCMC_Lhs, digits = 6)), 1, 0)

        AncConMat[matching_rows[j], "True_IsMax_Anc_IndMCMC"] <- ifelse(!is.null(grandparent), ifelse(round(Anc_IndMCMC_Lhs[results[matching_rows[j], "True_4States"]], digits = 6) == max(round(Anc_IndMCMC_Lhs, digits = 6)), 1, 0), NA)
        AncConMat[matching_rows[j], "True_IsMax_Anc_DepMCMC"] <- ifelse(!is.null(grandparent), ifelse(round(Anc_DepMCMC_Lhs[results[matching_rows[j], "True_4States"]], digits = 6) == max(round(Anc_DepMCMC_Lhs, digits = 6)), 1, 0), NA)

        # Now save whichever is the maximum likelihood
        AncConMat[matching_rows[j], "Max_Sis_IndMCMC_Lh"] <- max(round(Sis_IndMCMC_Lhs, digits = 6))
        AncConMat[matching_rows[j], "Max_Sis_DepMCMC_Lh"] <- max(round(Sis_DepMCMC_Lhs, digits = 6))

        AncConMat[matching_rows[j], "Max_Anc_IndMCMC_Lh"] <- ifelse(!is.null(grandparent), max(round(Anc_IndMCMC_Lhs, digits = 6)), NA)
        AncConMat[matching_rows[j], "Max_Anc_DepMCMC_Lh"] <- ifelse(!is.null(grandparent), max(round(Anc_DepMCMC_Lhs, digits = 6)), NA)
      }

      # Next RJ models, starting with the Ind model
      if (RJmodel == "RJMCMC" || RJmodel == "BOTH") {
        # First determine if the maximum likelihood from the ASR matches the unknown taxon's character state
        AncConMat[matching_rows[j], "True_IsMax_Sis_IndRJ"] <- ifelse(round(Sis_IndRJ_Lhs[results[matching_rows[j], "True_4States"]], digits = 6) == max(round(Sis_IndRJ_Lhs, digits = 6)), 1, 0)
        AncConMat[matching_rows[j], "True_IsMax_Sis_DepRJ"] <- ifelse(round(Sis_DepRJ_Lhs[results[matching_rows[j], "True_4States"]], digits = 6) == max(round(Sis_DepRJ_Lhs, digits = 6)), 1, 0)

        AncConMat[matching_rows[j], "True_IsMax_Anc_IndRJ"] <- ifelse(!is.null(grandparent), ifelse(round(Anc_IndRJ_Lhs[results[matching_rows[j], "True_4States"]], digits = 6) == max(round(Anc_IndRJ_Lhs, digits = 6)), 1, 0), NA)
        AncConMat[matching_rows[j], "True_IsMax_Anc_DepRJ"] <- ifelse(!is.null(grandparent), ifelse(round(Anc_DepRJ_Lhs[results[matching_rows[j], "True_4States"]], digits = 6) == max(round(Anc_DepRJ_Lhs, digits = 6)), 1, 0), NA)

        # Now save whichever is the maximum likelihood
        AncConMat[matching_rows[j], "Max_Sis_IndRJ_Lh"] <- max(round(Sis_IndRJ_Lhs, digits = 6))
        AncConMat[matching_rows[j], "Max_Sis_DepRJ_Lh"] <- max(round(Sis_DepRJ_Lhs, digits = 6))

        AncConMat[matching_rows[j], "Max_Anc_IndRJ_Lh"] <- ifelse(!is.null(grandparent), max(round(Anc_IndRJ_Lhs, digits = 6)), NA)
        AncConMat[matching_rows[j], "Max_Anc_DepRJ_Lh"] <- ifelse(!is.null(grandparent), max(round(Anc_DepRJ_Lhs, digits = 6)), NA)
      }


      ### Pull Terminal Branch Length and Tree Height
      AncConMat[matching_rows[j], "TBL"] <- TBL <- unknown_branch_length(rTaxon[j], tree_complete)
      AncConMat[matching_rows[j], "Tree_Height"] <- tree_height <- max(node.depth.edgelength(tree_complete))
      AncConMat[matching_rows[j], "Std_TBL"] <- TBL/tree_height

      ### Finally, add the trial information to this table
      AncConMat[matching_rows[j], "Trial_#"] <- i
      AncConMat[matching_rows[j], "Unknown_Taxon"] <- rTaxon[j]
      AncConMat[matching_rows[j], "True_A"] <- data[rTaxon[j], 2]
      AncConMat[matching_rows[j], "True_B"] <- data[rTaxon[j], 3]
      AncConMat[matching_rows[j], "True_4States"] <- data[rTaxon[j], 4]
    }
  }
  ### Save this table
  write.table(AncConMat, file = matrix_name, quote = F, row.names = F, col.names = T, sep = "\t")


  ### Add the important information into the main results folder, starting with Averaged Sister Trait Values
  results[, "Avg_Sister_A"] <- AncConMat[, "Avg_Sis_A"]
  results[, "Avg_Sister_B"] <- AncConMat[, "Avg_Sis_B"]

  ### Get the important information from the ancestral reconstructions, and add it to the main results file
  if (RJmodel == "MCMC" || RJmodel == "BOTH") {
    results[, "Sister_IndMCMC_MaxLh"] <- AncConMat[, "Sis_IndMCMC_Lh"]
    results[, "Sister_IndMCMC_Prediction"] <- AncConMat[, "Sis_IndMCMC_Pred"]
    results[, "Sister_DepMCMC_MaxLh"] <- AncConMat[, "Sis_DepMCMC_Lh"]
    results[, "Sister_DepMCMC_Prediction"] <- AncConMat[, "Sis_DepMCMC_Pred"]

    results[, "Ancestor_IndMCMC_MaxLh"] <- AncConMat[, "Anc_IndMCMC_Lh"]
    results[, "Ancestor_IndMCMC_Prediction"] <- AncConMat[, "Anc_IndMCMC_Pred"]
    results[, "Ancestor_DepMCMC_MaxLh"] <- AncConMat[, "Anc_DepMCMC_Lh"]
    results[, "Ancestor_DepMCMC_Prediction"] <- AncConMat[, "Anc_DepMCMC_Pred"]
  }
  if (RJmodel == "RJMCMC" || RJmodel == "BOTH") {
    results[, "Sister_IndRJ_MaxLh"] <- AncConMat[, "Sis_IndRJ_Lh"]
    results[, "Sister_IndRJ_Prediction"] <- AncConMat[, "Sis_IndRJ_Pred"]
    results[, "Sister_DepRJ_MaxLh"] <- AncConMat[, "Sis_DepRJ_Lh"]
    results[, "Sister_DepRJ_Prediction"] <- AncConMat[, "Sis_DepRJ_Pred"]

    results[, "Ancestor_IndRJ_MaxLh"] <- AncConMat[, "Anc_IndRJ_Lh"]
    results[, "Ancestor_IndRJ_Prediction"] <- AncConMat[, "Anc_IndRJ_Pred"]
    results[, "Ancestor_DepRJ_MaxLh"] <- AncConMat[, "Anc_DepRJ_Lh"]
    results[, "Ancestor_DepRJ_Prediction"] <- AncConMat[, "Anc_DepRJ_Pred"]
  }

  ### Now save the changes to the results matrices
  write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t")

  ### Finally, send a progress signal
  print(paste0("Finished with the ancestral reconstructions for ", type, " with Constant Rates"))
}




#####
### Now repeat for variable rates
if (variable_rates == TRUE) {
  for (type in types) {
    # Print a progress check
    print(paste0("Beginning the ancestral reconstruction for ", type, " with Variable Rates."))

    # Call the necessary paths
    data_path <- paste0("VariableRates/", type, "/Data/", type)
    results_name <- paste0("Results/VariableRates/Single/", type, ".Single.ResultsFull.txt")

    # Call the results table
    colnames_results <- strsplit(readLines(results_name, n = 1), "\t")[[1]]
    results <- read.table(results_name, skip = 1, sep = "\t", col.names = colnames_results)

    # Create a matrix to store the information
    matrix_name <- paste0("Results/VariableRates/Single/AncestralStateReconstruction/", type, ".AncReconAndTreeHeight.txt")

    if (RJmodel == "MCMC") {
      AncConMat <- matrix(data=NA, nrow = nrow(results), ncol = 42)
      colnames(AncConMat) <- c("Trial_#", "Unknown_Taxon", "True_A", "True_B", "True_4States", "Avg_Sis_A", "Avg_Sis_B",
                               "Sis_IndMCMC_00", "Sis_IndMCMC_01", "Sis_IndMCMC_10", "Sis_IndMCMC_11",
                               "Sis_DepMCMC_00", "Sis_DepMCMC_01", "Sis_DepMCMC_10", "Sis_DepMCMC_11",
                               "Sis_IndMCMC_Pred", "Sis_IndMCMC_Lh", "Sis_DepMCMC_Pred", "Sis_DepMCMC_Lh",
                               "Anc_IndMCMC_00", "Anc_IndMCMC_01", "Anc_IndMCMC_10", "Anc_IndMCMC_11",
                               "Anc_DepMCMC_00", "Anc_DepMCMC_01", "Anc_DepMCMC_10", "Anc_DepMCMC_11",
                               "Anc_IndMCMC_Pred", "Anc_IndMCMC_Lh", "Anc_DepMCMC_Pred", "Anc_DepMCMC_Lh",
                               "True_IsMax_Sis_IndMCMC", "Max_Sis_IndMCMC_Lh", "True_IsMax_Sis_DepMCMC", "Max_Sis_DepMCMC_Lh",
                               "True_IsMax_Anc_IndMCMC", "Max_Anc_IndMCMC_Lh", "True_IsMax_Anc_DepMCMC", "Max_Anc_DepMCMC_Lh",
                               "TBL", "Tree_Height", "Std_TBL")
    } else if (RJmodel == "RJMCMC") {
      AncConMat <- matrix(data=NA, nrow = nrow(results), ncol = 42)
      colnames(AncConMat) <- c("Trial_#", "Unknown_Taxon", "True_A", "True_B", "True_4States", "Avg_Sis_A", "Avg_Sis_B",
                               "Sis_IndRJ_00", "Sis_IndRJ_01", "Sis_IndRJ_10", "Sis_IndRJ_11",
                               "Sis_DepRJ_00", "Sis_DepRJ_01", "Sis_DepRJ_10", "Sis_DepRJ_11",
                               "Sis_IndRJ_Pred", "Sis_IndRJ_Lh", "Sis_DepRJ_Pred", "Sis_DepRJ_Lh",
                               "Anc_IndRJ_00", "Anc_IndRJ_01", "Anc_IndRJ_10", "Anc_IndRJ_11",
                               "Anc_DepRJ_00", "Anc_DepRJ_01", "Anc_DepRJ_10", "Anc_DepRJ_11",
                               "Anc_IndRJ_Pred", "Anc_IndRJ_Lh", "Anc_DepRJ_Pred", "Anc_DepRJ_Lh",
                               "True_IsMax_Sis_IndRJ", "Max_Sis_IndRJ_Lh", "True_IsMax_Sis_DepRJ", "Max_Sis_DepRJ_Lh",
                               "True_IsMax_Anc_IndRJ", "Max_Anc_IndRJ_Lh", "True_IsMax_Anc_DepRJ", "Max_Anc_DepRJ_Lh",
                               "TBL", "Tree_Height", "Std_TBL")
    } else if (RJmodel == "BOTH") {
      AncConMat <- matrix(data=NA, nrow = nrow(results), ncol = 74)
      colnames(AncConMat) <- c("Trial_#", "Unknown_Taxon", "True_A", "True_B", "True_4States", "Avg_Sis_A", "Avg_Sis_B",
                               "Sis_IndMCMC_00", "Sis_IndMCMC_01", "Sis_IndMCMC_10", "Sis_IndMCMC_11",
                               "Sis_DepMCMC_00", "Sis_DepMCMC_01", "Sis_DepMCMC_10", "Sis_DepMCMC_11",
                               "Sis_IndRJ_00", "Sis_IndRJ_01", "Sis_IndRJ_10", "Sis_IndRJ_11",
                               "Sis_DepRJ_00", "Sis_DepRJ_01", "Sis_DepRJ_10", "Sis_DepRJ_11",
                               "Sis_IndMCMC_Pred", "Sis_IndMCMC_Lh", "Sis_DepMCMC_Pred", "Sis_DepMCMC_Lh",
                               "Sis_IndRJ_Pred", "Sis_IndRJ_Lh", "Sis_DepRJ_Pred", "Sis_DepRJ_Lh",
                               "Anc_IndMCMC_00", "Anc_IndMCMC_01", "Anc_IndMCMC_10", "Anc_IndMCMC_11",
                               "Anc_DepMCMC_00", "Anc_DepMCMC_01", "Anc_DepMCMC_10", "Anc_DepMCMC_11",
                               "Anc_IndRJ_00", "Anc_IndRJ_01", "Anc_IndRJ_10", "Anc_IndRJ_11",
                               "Anc_DepRJ_00", "Anc_DepRJ_01", "Anc_DepRJ_10", "Anc_DepRJ_11",
                               "Anc_IndMCMC_Pred", "Anc_IndMCMC_Lh", "Anc_DepMCMC_Pred", "Anc_DepMCMC_Lh",
                               "Anc_IndRJ_Pred", "Anc_IndRJ_Lh", "Anc_DepRJ_Pred", "Anc_DepRJ_Lh",
                               "True_IsMax_Sis_IndMCMC", "Max_Sis_IndMCMC_Lh", "True_IsMax_Sis_DepMCMC", "Max_Sis_DepMCMC_Lh",
                               "True_IsMax_Sis_IndRJ", "Max_Sis_IndRJ_Lh", "True_IsMax_Sis_DepRJ", "Max_Sis_DepRJ_Lh",
                               "True_IsMax_Anc_IndMCMC", "Max_Anc_IndMCMC_Lh", "True_IsMax_Anc_DepMCMC", "Max_Anc_DepMCMC_Lh",
                               "True_IsMax_Anc_IndRJ", "Max_Anc_IndRJ_Lh", "True_IsMax_Anc_DepRJ", "Max_Anc_DepRJ_Lh",
                               "TBL", "Tree_Height", "Std_TBL")
    }



    # Start the for-loop that will run each prediction and fill it in
    for (i in 1:trial_amount) {
      # First, call rTaxon by subsetting the rows where column 1 matches i
      matching_rows <- which(results[, "Trial"] == i)
      rTaxon <- results[matching_rows, "rTaxon"]

      # Call the data file for this iteration
      data_name <- paste0(data_path, ".", i, ".Full_data.txt")
      data <- read.table(data_name, skip = 1, sep = "\t")
      data[, 4] <- match(data[,4], c(0,1,10,11))


      # Pull the tree to read the terminal branch length
      treename <- paste0("Trees/Full_tree.", i, ".tre")
      tree_complete <- read.nexus(treename)

      # Loop through every taxon for which we are predicting
      for (j in 1:length(rTaxon)) {
        # This collects the numerical average value of the sister taxon's trait information
        AncConMat[matching_rows[j], "Avg_Sis_A"] <- sisterData(tree_complete, rTaxon[j], data, "A")
        AncConMat[matching_rows[j], "Avg_Sis_B"] <- sisterData(tree_complete, rTaxon[j], data, "B")

        # Reset the sister node number and find the other important nodes
        sisters <- getSisters(tree_complete, rTaxon[j], mode = "number")
        sister_tips <- getAllDescendantTips(tree_complete, sisters)
        grandparent <- getParent(tree_complete, getMRCA(tree_complete, c(rTaxon, sister_tips)))
        cousins <- getAllDescendantTips(tree_complete, grandparent)
        cousins <- cousins[cousins != rTaxon[j]]

        # Find the node of interest in a tree where you don't know the unknown's information
        edited_tree <- drop.tip(tree_complete, rTaxon[j])
        edited_data <- data[-rTaxon,]
        edited_cousins <- match(cousins, edited_tree$tip.label)
        grandparent_node_asr <- getMRCA(edited_tree, edited_cousins)



        ### Ancestral State Reconstructin with both traits
        if (as.numeric(sisters) > Ntip(tree_complete) && length(unique(data[, 4])) > 1) {
          # Adjust node number to row number of the output
          sister_asr <- getMRCA(edited_tree, sister_tips)
          sister_asr <- sister_asr - Ntip(edited_tree)

          ### First pull the Independent MCMC model, and perform ML Ancestral State Reconstruction on the nearest sister taxon
          if (RJmodel == "MCMC" || RJmodel == "BOTH") {
            IndMCMC_qmat <- estimated_qmat(type, i, IndependentRates = TRUE, ConstantRates = TRUE, RJ = FALSE, Predict = FALSE)
            ace_IndMCMC <- asr_mk_model(tree = edited_tree, tip_states = edited_data[, 4], Nstates = 4, transition_matrix = IndMCMC_qmat, reroot = FALSE, root_prior = "flat")
            # Pull the likelihood
            Sis_IndMCMC_Lhs <- round(ace_IndMCMC$ancestral_likelihoods[sister_asr, ], digits = 6)
            # Store them in the matrix
            AncConMat[matching_rows[j], c("Sis_IndMCMC_00", "Sis_IndMCMC_01", "Sis_IndMCMC_10", "Sis_IndMCMC_11")] <- Sis_IndMCMC_Lhs

          ### Next Dependent MCMC models
            DepMCMC_qmat <- estimated_qmat(type, i, IndependentRates = FALSE, ConstantRates = T, RJ = FALSE, Predict = FALSE)
            ace_DepMCMC <- asr_mk_model(tree = edited_tree, tip_states = edited_data[, 4], Nstates = 4, transition_matrix = DepMCMC_qmat, reroot = FALSE, root_prior = "flat")
            # Pull the likelihood
            Sis_DepMCMC_Lhs <- round(ace_DepMCMC$ancestral_likelihoods[sister_asr, ], digits = 6)
            # Store them in the matrix
            AncConMat[matching_rows[j], c("Sis_DepMCMC_00", "Sis_DepMCMC_01", "Sis_DepMCMC_10", "Sis_DepMCMC_11")] <- Sis_DepMCMC_Lhs
          }


          ### Next Ind RJ models
          if (RJmodel == "RJMCMC" || RJmodel == "BOTH") {
            IndRJ_qmat <- estimated_qmat(type, i, IndependentRates = TRUE, ConstantRates = T, RJ = T, Predict = FALSE)
            ace_IndRJ <- asr_mk_model(tree = edited_tree, tip_states = edited_data[, 4], Nstates = 4, transition_matrix = IndRJ_qmat, reroot = FALSE, root_prior = "flat")
            # Pull the likelihood
            Sis_IndRJ_Lhs <- round(ace_IndRJ$ancestral_likelihoods[sister_asr, ], digits = 6)
            # Store them in the matrix
            AncConMat[matching_rows[j], c("Sis_IndRJ_00", "Sis_IndRJ_01", "Sis_IndRJ_10", "Sis_IndRJ_11")] <- Sis_IndRJ_Lhs

          ### Finally, the Dep RJ models
            DepRJ_qmat <- estimated_qmat(type, i, IndependentRates = TRUE, ConstantRates = T, RJ = FALSE, Predict = FALSE)
            ace_DepRJ <- asr_mk_model(tree = edited_tree, tip_states = edited_data[, 4], Nstates = 4, transition_matrix = DepRJ_qmat, reroot = FALSE, root_prior = "flat")
            # Pull the likelihood
            Sis_DepRJ_Lhs <- round(ace_DepRJ$ancestral_likelihoods[sister_asr, ], digits = 6)
            # Store them in the matrix
            AncConMat[matching_rows[j], c("Sis_DepRJ_00", "Sis_DepRJ_01", "Sis_DepRJ_10", "Sis_DepRJ_11")] <- Sis_DepRJ_Lhs
          }


          ### Line accounting for non-variation in the tree or single tip prediction
        } else {

          # In case of non-variation in the tree, make this adjustment for sister-node reconstruction
          # This only impacts cases of no trait variation across the tree when the sister node is internal
          # Although wonky, it covers edge cases without impacting sister-leaf reconstruction
          if (as.numeric(sisters) > Ntip(tree_complete)) {sisters <- sisters - Ntip(tree_complete)}

          # First cover the MCMC models
          if (RJmodel == "MCMC" || RJmodel == "BOTH") {
            Sis_IndMCMC_Lhs <- NULL
            Sis_IndMCMC_Lhs[1] <- Sis_IndMCMC_Lhs[2] <- Sis_IndMCMC_Lhs[3] <- Sis_IndMCMC_Lhs[4] <- 0
            Sis_IndMCMC_Lhs[ data[sisters,4] ] <- 1
            AncConMat[matching_rows[j], c("Sis_IndMCMC_00", "Sis_IndMCMC_01", "Sis_IndMCMC_10", "Sis_IndMCMC_11")] <- Sis_IndMCMC_Lhs

            Sis_DepMCMC_Lhs <- NULL
            Sis_DepMCMC_Lhs[1] <- Sis_DepMCMC_Lhs[2] <- Sis_DepMCMC_Lhs[3] <- Sis_DepMCMC_Lhs[4] <- 0
            Sis_DepMCMC_Lhs[ data[sisters,4] ] <- 1
            AncConMat[matching_rows[j], c("Sis_DepMCMC_00", "Sis_DepMCMC_01", "Sis_DepMCMC_10", "Sis_DepMCMC_11")] <- Sis_DepMCMC_Lhs
          }

          # Next cover the RJMCMC models
          if (RJmodel == "RJMCMC" || RJmodel == "BOTH") {
            Sis_IndRJ_Lhs <- NULL
            Sis_IndRJ_Lhs[1] <- Sis_IndRJ_Lhs[2] <- Sis_IndRJ_Lhs[3] <- Sis_IndRJ_Lhs[4] <- 0
            Sis_IndRJ_Lhs[ data[sisters,4] ] <- 1
            AncConMat[matching_rows[j], c("Sis_IndRJ_00", "Sis_IndRJ_01", "Sis_IndRJ_10", "Sis_IndRJ_11")] <- Sis_IndRJ_Lhs

            Sis_DepRJ_Lhs <- NULL
            Sis_DepRJ_Lhs[1] <- Sis_DepRJ_Lhs[2] <- Sis_DepRJ_Lhs[3] <- Sis_DepRJ_Lhs[4] <- 0
            Sis_DepRJ_Lhs[ data[sisters,4] ] <- 1
            AncConMat[matching_rows[j], c("Sis_DepRJ_00", "Sis_DepRJ_01", "Sis_DepRJ_10", "Sis_DepRJ_11")] <- Sis_DepRJ_Lhs
          }
        }


        ### Now pull the information for the ancestral node (the grandpa)
        # First, adjust to match the internal node number
        if (!is.null(grandparent)) {
          grandparent_node_asr <- grandparent_node_asr - Ntip(edited_tree)
        }


        ### Perform ML Ancestral State Reconstruction on the parent node of the Phylogenetic Bracket
        # First the MCMC models
        if ((RJmodel == "MCMC" || RJmodel == "BOTH") && !is.null(grandparent)) {
          IndMCMC_qmat <- estimated_qmat(type, i, IndependentRates = TRUE, ConstantRates = TRUE, RJ = FALSE, Predict = FALSE)
          Anc_ace_IndMCMC <- asr_mk_model(tree = edited_tree, tip_states = edited_data[, 4], Nstates = 4, transition_matrix = IndMCMC_qmat, reroot = FALSE, root_prior = "flat")
          # Pull the likelihood
          Anc_IndMCMC_Lhs <- round(Anc_ace_IndMCMC$ancestral_likelihoods[grandparent_node_asr, ], digits = 6)
          # Store them in the matrix
          AncConMat[matching_rows[j], c("Anc_IndMCMC_00", "Anc_IndMCMC_01", "Anc_IndMCMC_10", "Anc_IndMCMC_11")] <- Anc_IndMCMC_Lhs

          ### Next Dependent MCMC models
          DepMCMC_qmat <- estimated_qmat(type, i, IndependentRates = FALSE, ConstantRates = T, RJ = FALSE, Predict = FALSE)
          Anc_ace_DepMCMC <- asr_mk_model(tree = edited_tree, tip_states = edited_data[, 4], Nstates = 4, transition_matrix = DepMCMC_qmat, reroot = FALSE, root_prior = "flat")
          # Pull the likelihood
          Anc_DepMCMC_Lhs <- round(Anc_ace_DepMCMC$ancestral_likelihoods[grandparent_node_asr, ], digits = 6)
          # Store them in the matrix
          AncConMat[matching_rows[j], c("Anc_DepMCMC_00", "Anc_DepMCMC_01", "Anc_DepMCMC_10", "Anc_DepMCMC_11")] <- Anc_DepMCMC_Lhs
        }

        ### Next RJ models, starting with the Ind model
        if ((RJmodel == "RJMCMC" || RJmodel == "BOTH") && !is.null(grandparent)) {
          IndRJ_qmat <- estimated_qmat(type, i, IndependentRates = TRUE, ConstantRates = T, RJ = T, Predict = FALSE)
          Anc_ace_IndRJ <- asr_mk_model(tree = edited_tree, tip_states = edited_data[, 4], Nstates = 4, transition_matrix = IndRJ_qmat, reroot = FALSE, root_prior = "flat")
          # Pull the likelihood
          Anc_IndRJ_Lhs <- round(Anc_ace_IndRJ$ancestral_likelihoods[grandparent_node_asr, ], digits = 6)
          # Store them in the matrix
          AncConMat[matching_rows[j], c("Anc_IndRJ_00", "Anc_IndRJ_01", "Anc_IndRJ_10", "Anc_IndRJ_11")] <- Anc_IndRJ_Lhs

          ### Finally, the Dep RJ models
          DepRJ_qmat <- estimated_qmat(type, i, IndependentRates = F, ConstantRates = T, RJ = T, Predict = FALSE)
          Anc_ace_DepRJ <- asr_mk_model(tree = edited_tree, tip_states = edited_data[, 4], Nstates = 4, transition_matrix = DepRJ_qmat, reroot = FALSE, root_prior = "flat")
          # Pull the likelihood
          Anc_DepRJ_Lhs <- round(Anc_ace_DepRJ$ancestral_likelihoods[grandparent_node_asr, ], digits = 6)
          # Store them in the matrix
          AncConMat[matching_rows[j], c("Anc_DepRJ_00", "Anc_DepRJ_01", "Anc_DepRJ_10", "Anc_DepRJ_11")] <- Anc_DepRJ_Lhs
        }


        ### Adjust the code if the parent node is the root
        # First the MCMC models
        if ((RJmodel == "MCMC" || RJmodel == "BOTH") && is.null(grandparent)) {
          # Reset the likelihoods to clarify this trial does not have a valid "grandparent"
          Anc_IndMCMC_Lhs <- c(NA, NA, NA, NA)
          # Store them in the matrix
          AncConMat[matching_rows[j], c("Anc_IndMCMC_00", "Anc_IndMCMC_01", "Anc_IndMCMC_10", "Anc_IndMCMC_11")] <- Anc_IndMCMC_Lhs

          ### Next Dependent MCMC models
          # Reset the likelihoods to clarify this trial does not have a valid "grandparent"
          Anc_DepMCMC_Lhs <- c(NA, NA, NA, NA)
          # Store them in the matrix
          AncConMat[matching_rows[j], c("Anc_DepMCMC_00", "Anc_DepMCMC_01", "Anc_DepMCMC_10", "Anc_DepMCMC_11")] <- Anc_DepMCMC_Lhs
        }

        ### Next RJ models, starting with the Ind model
        if ((RJmodel == "RJMCMC" || RJmodel == "BOTH") && is.null(grandparent)) {
          # Reset the likelihoods to clarify this trial does not have a valid "grandparent"
          Anc_IndRJ_Lhs <- c(NA, NA, NA, NA)
          # Store them in the matrix
          AncConMat[matching_rows[j], c("Anc_IndRJ_00", "Anc_IndRJ_01", "Anc_IndRJ_10", "Anc_IndRJ_11")] <- Anc_IndRJ_Lhs

          ### Finally, the Dep RJ models
          # Reset the likelihoods to clarify this trial does not have a valid "grandparent"
          Anc_DepRJ_Lhs <- c(NA, NA, NA, NA)
          # Store them in the matrix
          AncConMat[matching_rows[j], c("Anc_DepRJ_00", "Anc_DepRJ_01", "Anc_DepRJ_10", "Anc_DepRJ_11")] <- Anc_DepRJ_Lhs
        }


        ### Pull the predicted states and likelihoods for the tests performed
        if ((RJmodel == "MCMC" || RJmodel == "BOTH") && !is.null(grandparent)) {
          # First for the sister
          AncConMat[matching_rows[j], "Sis_IndMCMC_Pred"] <- which.max(Sis_IndMCMC_Lhs)
          AncConMat[matching_rows[j], "Sis_IndMCMC_Lh"] <- max(Sis_IndMCMC_Lhs)


          AncConMat[matching_rows[j], "Sis_DepMCMC_Pred"] <- which.max(Sis_DepMCMC_Lhs)
          AncConMat[matching_rows[j], "Sis_DepMCMC_Lh"] <- max(Sis_DepMCMC_Lhs)

          # Now for the grandparent
          AncConMat[matching_rows[j], "Anc_IndMCMC_Pred"] <- which.max(Anc_IndMCMC_Lhs)
          AncConMat[matching_rows[j], "Anc_IndMCMC_Lh"] <- max(Anc_IndMCMC_Lhs)

          AncConMat[matching_rows[j], "Anc_DepMCMC_Pred"] <- which.max(Anc_DepMCMC_Lhs)
          AncConMat[matching_rows[j], "Anc_DepMCMC_Lh"] <- max(Anc_DepMCMC_Lhs)

          # Separate out the likelihood of the "true" state of rTaxon
          results[matching_rows[j], "TrueState_Sister_IndMCMC_Lh"] <- Sis_IndMCMC_Lhs[results[matching_rows[j], "True_4States"]]
          results[matching_rows[j], "TrueState_Sister_DepMCMC_Lh"] <- Sis_DepMCMC_Lhs[results[matching_rows[j], "True_4States"]]

          results[matching_rows[j], "TrueState_Ancestor_IndMCMC_Lh"] <- Anc_IndMCMC_Lhs[results[matching_rows[j], "True_4States"]]
          results[matching_rows[j], "TrueState_Ancestor_DepMCMC_Lh"] <- Anc_DepMCMC_Lhs[results[matching_rows[j], "True_4States"]]
        }

        if ((RJmodel == "RJMCMC" || RJmodel == "BOTH") && !is.null(grandparent)) {
          # First for the sister
          AncConMat[matching_rows[j], "Sis_IndRJ_Pred"] <- which.max(Sis_IndRJ_Lhs)
          AncConMat[matching_rows[j], "Sis_IndRJ_Lh"] <- max(Sis_IndRJ_Lhs)

          AncConMat[matching_rows[j], "Sis_DepRJ_Pred"] <- which.max(Sis_DepRJ_Lhs)
          AncConMat[matching_rows[j], "Sis_DepRJ_Lh"] <- max(Sis_DepRJ_Lhs)

          # Now for the grandparent
          AncConMat[matching_rows[j], "Anc_IndRJ_Pred"] <- which.max(Anc_IndRJ_Lhs)
          AncConMat[matching_rows[j], "Anc_IndRJ_Lh"] <- max(Anc_IndRJ_Lhs)

          AncConMat[matching_rows[j], "Anc_DepRJ_Pred"] <- which.max(Anc_DepRJ_Lhs)
          AncConMat[matching_rows[j], "Anc_DepRJ_Lh"] <- max(Anc_DepRJ_Lhs)

          # Separate out the likelihood of the "true" state of rTaxon
          results[matching_rows[j], "TrueState_Sister_IndRJ_Lh"] <- Sis_IndMCMC_Lhs[results[matching_rows[j], "True_4States"]]
          results[matching_rows[j], "TrueState_Sister_DepRJ_Lh"] <- Sis_DepMCMC_Lhs[results[matching_rows[j], "True_4States"]]

          results[matching_rows[j], "TrueState_Ancestor_IndRJ_Lh"] <- Anc_IndMCMC_Lhs[results[matching_rows[j], "True_4States"]]
          results[matching_rows[j], "TrueState_Ancestor_DepRJ_Lh"] <- Anc_DepMCMC_Lhs[results[matching_rows[j], "True_4States"]]
        }


        ### A fix for when there is no grandparent node (parent is the root)
        # These are again separated by tests run
        if ((RJmodel == "MCMC" || RJmodel == "BOTH") && is.null(grandparent)) {
          # Pull the predicted states and likelihoods for the sister
          AncConMat[matching_rows[j], "Sis_IndMCMC_Pred"] <- which.max(Sis_IndMCMC_Lhs)
          AncConMat[matching_rows[j], "Sis_IndMCMC_Lh"] <- max(Sis_IndMCMC_Lhs)

          AncConMat[matching_rows[j], "Sis_DepMCMC_Pred"] <- which.max(Sis_DepMCMC_Lhs)
          AncConMat[matching_rows[j], "Sis_DepMCMC_Lh"] <- max(Sis_DepMCMC_Lhs)

          # Separate out the likelihood of the "true" state of rTaxon
          results[matching_rows[j], "TrueState_Sister_IndMCMC_Lh"] <- Sis_IndMCMC_Lhs[results[matching_rows[j], "True_4States"]]
          results[matching_rows[j], "TrueState_Sister_DepMCMC_Lh"] <- Sis_DepMCMC_Lhs[results[matching_rows[j], "True_4States"]]

          # However, none of the "ancestor" state reconstruction can be done without a known ancestor.
          AncConMat[matching_rows[j], "Anc_IndMCMC_Pred"] <- NA
          AncConMat[matching_rows[j], "Anc_IndMCMC_Lh"] <- NA

          AncConMat[matching_rows[j], "Anc_DepMCMC_Pred"] <- NA
          AncConMat[matching_rows[j], "Anc_DepMCMC_Lh"] <- NA

          results[matching_rows[j], "TrueState_Ancestor_IndMCMC_Lh"] <- NA
          results[matching_rows[j], "TrueState_Ancestor_DepMCMC_Lh"] <- NA
        }

        if ((RJmodel == "RJMCMC" || RJmodel == "BOTH") && is.null(grandparent)) {
          # Pull the predicted states and likelihoods for the sister
          AncConMat[matching_rows[j], "Sis_IndRJ_Pred"] <- which.max(Sis_IndRJ_Lhs)
          AncConMat[matching_rows[j], "Sis_IndRJ_Lh"] <- max(Sis_IndRJ_Lhs)

          AncConMat[matching_rows[j], "Sis_DepRJ_Pred"] <- which.max(Sis_DepRJ_Lhs)
          AncConMat[matching_rows[j], "Sis_DepRJ_Lh"] <- max(Sis_DepRJ_Lhs)

          # Separate out the likelihood of the "true" state of rTaxon
          results[matching_rows[j], "TrueState_Sister_IndRJ_Lh"] <- Sis_IndRJ_Lhs[results[matching_rows[j], "True_4States"]]
          results[matching_rows[j], "TrueState_Sister_DepRJ_Lh"] <- Sis_DepRJ_Lhs[results[matching_rows[j], "True_4States"]]

          # However, none of the "ancestor" state reconstruction can be done without a known ancestor.
          AncConMat[matching_rows[j], "Anc_IndRJ_Pred"] <- NA
          AncConMat[matching_rows[j], "Anc_IndRJ_Lh"] <- NA

          AncConMat[matching_rows[j], "Anc_DepRJ_Pred"] <- NA
          AncConMat[matching_rows[j], "Anc_DepRJ_Lh"] <- NA

          results[matching_rows[j], "TrueState_Ancestor_IndRJ_Lh"] <- NA
          results[matching_rows[j], "TrueState_Ancestor_DepRJ_Lh"] <- NA
        }


        ### These lines add to the results whether or not the ACR would agree with our predictions or the true, unknown state
        # First the MCMC models
        if (RJmodel == "MCMC" || RJmodel == "BOTH") {
          # First determine if the maximum likelihood from the ASR matches the unknown taxon's character state
          AncConMat[matching_rows[j], "True_IsMax_Sis_IndMCMC"] <- ifelse(Sis_IndMCMC_Lhs[results[matching_rows[j], "True_4States"]] == max(round(Sis_IndMCMC_Lhs, digits = 6)), 1, 0)
          AncConMat[matching_rows[j], "True_IsMax_Sis_DepMCMC"] <- ifelse(Sis_DepMCMC_Lhs[results[matching_rows[j], "True_4States"]] == max(round(Sis_DepMCMC_Lhs, digits = 6)), 1, 0)

          AncConMat[matching_rows[j], "True_IsMax_Anc_IndMCMC"] <- ifelse(!is.null(grandparent), ifelse(Anc_IndMCMC_Lhs[results[matching_rows[j], "True_4States"]] == max(round(Anc_IndMCMC_Lhs, digits = 6)), 1, 0), NA)
          AncConMat[matching_rows[j], "True_IsMax_Anc_DepMCMC"] <- ifelse(!is.null(grandparent), ifelse(Anc_DepMCMC_Lhs[results[matching_rows[j], "True_4States"]] == max(round(Anc_DepMCMC_Lhs, digits = 6)), 1, 0), NA)

          # Now save whichever is the maximum likelihood
          AncConMat[matching_rows[j], "Max_Sis_IndMCMC_Lh"] <- max(round(Sis_IndMCMC_Lhs, digits = 6))
          AncConMat[matching_rows[j], "Max_Sis_DepMCMC_Lh"] <- max(round(Sis_DepMCMC_Lhs, digits = 6))

          AncConMat[matching_rows[j], "Max_Anc_IndMCMC_Lh"] <- ifelse(!is.null(grandparent), max(round(Anc_IndMCMC_Lhs, digits = 6)), NA)
          AncConMat[matching_rows[j], "Max_Anc_DepMCMC_Lh"] <- ifelse(!is.null(grandparent), max(round(Anc_DepMCMC_Lhs, digits = 6)), NA)
        }

        # Next RJ models, starting with the Ind model
        if (RJmodel == "RJMCMC" || RJmodel == "BOTH") {
          # First determine if the maximum likelihood from the ASR matches the unknown taxon's character state
          AncConMat[matching_rows[j], "True_IsMax_Sis_IndRJ"] <- ifelse(Sis_IndRJ_Lhs[results[matching_rows[j], "True_4States"]] == max(round(Sis_IndRJ_Lhs, digits = 6)), 1, 0)
          AncConMat[matching_rows[j], "True_IsMax_Sis_DepRJ"] <- ifelse(Sis_DepRJ_Lhs[results[matching_rows[j], "True_4States"]] == max(round(Sis_DepRJ_Lhs, digits = 6)), 1, 0)

          AncConMat[matching_rows[j], "True_IsMax_Anc_IndRJ"] <- ifelse(!is.null(grandparent), ifelse(Anc_IndRJ_Lhs[results[matching_rows[j], "True_4States"]] == max(round(Anc_IndRJ_Lhs, digits = 6)), 1, 0), NA)
          AncConMat[matching_rows[j], "True_IsMax_Anc_DepRJ"] <- ifelse(!is.null(grandparent), ifelse(Anc_DepRJ_Lhs[results[matching_rows[j], "True_4States"]] == max(round(Anc_DepRJ_Lhs, digits = 6)), 1, 0), NA)

          # Now save whichever is the maximum likelihood
          AncConMat[matching_rows[j], "Max_Sis_IndRJ_Lh"] <- max(round(Sis_IndRJ_Lhs, digits = 6))
          AncConMat[matching_rows[j], "Max_Sis_DepRJ_Lh"] <- max(round(Sis_DepRJ_Lhs, digits = 6))

          AncConMat[matching_rows[j], "Max_Anc_IndRJ_Lh"] <- ifelse(!is.null(grandparent), max(round(Anc_IndRJ_Lhs, digits = 6)), NA)
          AncConMat[matching_rows[j], "Max_Anc_DepRJ_Lh"] <- ifelse(!is.null(grandparent), max(round(Anc_DepRJ_Lhs, digits = 6)), NA)
        }


        ### Pull Terminal Branch Length and Tree Height
        AncConMat[matching_rows[j], "TBL"] <- TBL <- unknown_branch_length(rTaxon[j], tree_complete)
        AncConMat[matching_rows[j], "Tree_Height"] <- tree_height <- max(node.depth.edgelength(tree_complete))
        AncConMat[matching_rows[j], "Std_TBL"] <- TBL/tree_height

        ### Finally, add the trial information to this table
        AncConMat[matching_rows[j], "Trial_#"] <- i
        AncConMat[matching_rows[j], "Unknown_Taxon"] <- rTaxon[j]
        AncConMat[matching_rows[j], "True_A"] <- data[rTaxon[j], 2]
        AncConMat[matching_rows[j], "True_B"] <- data[rTaxon[j], 3]
        AncConMat[matching_rows[j], "True_4States"] <- data[rTaxon[j], 4]
      }
    }
    ### Save this table
    write.table(AncConMat, file = matrix_name, quote = F, row.names = F, col.names = T, sep = "\t")


    ### Add the important information into the main results folder, starting with Averaged Sister Trait Values
    results[, "Avg_Sister_A"] <- AncConMat[, "Avg_Sis_A"]
    results[, "Avg_Sister_B"] <- AncConMat[, "Avg_Sis_B"]

    ### Get the important information from the ancestral reconstructions, and add it to the main results file
    if (RJmodel == "MCMC" || RJmodel == "BOTH") {
      results[, "Sister_IndMCMC_MaxLh"] <- AncConMat[, "Sis_IndMCMC_Lh"]
      results[, "Sister_IndMCMC_Prediction"] <- AncConMat[, "Sis_IndMCMC_Pred"]
      results[, "Sister_DepMCMC_MaxLh"] <- AncConMat[, "Sis_DepMCMC_Lh"]
      results[, "Sister_DepMCMC_Prediction"] <- AncConMat[, "Sis_DepMCMC_Pred"]

      results[, "Ancestor_IndMCMC_MaxLh"] <- AncConMat[, "Anc_IndMCMC_Lh"]
      results[, "Ancestor_IndMCMC_Prediction"] <- AncConMat[, "Anc_IndMCMC_Pred"]
      results[, "Ancestor_DepMCMC_MaxLh"] <- AncConMat[, "Anc_DepMCMC_Lh"]
      results[, "Ancestor_DepMCMC_Prediction"] <- AncConMat[, "Anc_DepMCMC_Pred"]
    }
    if (RJmodel == "RJMCMC" || RJmodel == "BOTH") {
      results[, "Sister_IndRJ_MaxLh"] <- AncConMat[, "Sis_IndRJ_Lh"]
      results[, "Sister_IndRJ_Prediction"] <- AncConMat[, "Sis_IndRJ_Pred"]
      results[, "Sister_DepRJ_MaxLh"] <- AncConMat[, "Sis_DepRJ_Lh"]
      results[, "Sister_DepRJ_Prediction"] <- AncConMat[, "Sis_DepRJ_Pred"]

      results[, "Ancestor_IndRJ_MaxLh"] <- AncConMat[, "Anc_IndRJ_Lh"]
      results[, "Ancestor_IndRJ_Prediction"] <- AncConMat[, "Anc_IndRJ_Pred"]
      results[, "Ancestor_DepRJ_MaxLh"] <- AncConMat[, "Anc_DepRJ_Lh"]
      results[, "Ancestor_DepRJ_Prediction"] <- AncConMat[, "Anc_DepRJ_Pred"]
    }

    # Now save the changes to the results matrices
    write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t")

    # Finally, send a progress signal
    print(paste0("Finished with the ancestral reconstructions for ", type, " with Variable Rates"))
  }
}




#####
### Finally, repeat this for Random Data
# Print a progress check
print(paste0("Beginning the ancestral reconstruction for Random Rates."))

# Call the necessary paths
data_path <- "Random/Data/Random"
results_name <- "Results/Random/Random.Single.ResultsFull.txt"

# Call the results table
colnames_results <- strsplit(readLines(results_name, n = 1), "\t")[[1]]
results <- read.table(results_name, skip = 1, sep = "\t", col.names = colnames_results)

# Create a matrix to store the information
matrix_name <- "Results/Random/AncestralStateReconstruction/Random.AncReconAndTreeHeight.txt"

if (RJmodel == "MCMC") {
  AncConMat <- matrix(data=NA, nrow = nrow(results), ncol = 42)
  colnames(AncConMat) <- c("Trial_#", "Unknown_Taxon", "True_A", "True_B", "True_4States", "Avg_Sis_A", "Avg_Sis_B",
                           "Sis_IndMCMC_00", "Sis_IndMCMC_01", "Sis_IndMCMC_10", "Sis_IndMCMC_11",
                           "Sis_DepMCMC_00", "Sis_DepMCMC_01", "Sis_DepMCMC_10", "Sis_DepMCMC_11",
                           "Sis_IndMCMC_Pred", "Sis_IndMCMC_Lh", "Sis_DepMCMC_Pred", "Sis_DepMCMC_Lh",
                           "Anc_IndMCMC_00", "Anc_IndMCMC_01", "Anc_IndMCMC_10", "Anc_IndMCMC_11",
                           "Anc_DepMCMC_00", "Anc_DepMCMC_01", "Anc_DepMCMC_10", "Anc_DepMCMC_11",
                           "Anc_IndMCMC_Pred", "Anc_IndMCMC_Lh", "Anc_DepMCMC_Pred", "Anc_DepMCMC_Lh",
                           "True_IsMax_Sis_IndMCMC", "Max_Sis_IndMCMC_Lh", "True_IsMax_Sis_DepMCMC", "Max_Sis_DepMCMC_Lh",
                           "True_IsMax_Anc_IndMCMC", "Max_Anc_IndMCMC_Lh", "True_IsMax_Anc_DepMCMC", "Max_Anc_DepMCMC_Lh",
                           "TBL", "Tree_Height", "Std_TBL")
} else if (RJmodel == "RJMCMC") {
  AncConMat <- matrix(data=NA, nrow = nrow(results), ncol = 42)
  colnames(AncConMat) <- c("Trial_#", "Unknown_Taxon", "True_A", "True_B", "True_4States", "Avg_Sis_A", "Avg_Sis_B",
                           "Sis_IndRJ_00", "Sis_IndRJ_01", "Sis_IndRJ_10", "Sis_IndRJ_11",
                           "Sis_DepRJ_00", "Sis_DepRJ_01", "Sis_DepRJ_10", "Sis_DepRJ_11",
                           "Sis_IndRJ_Pred", "Sis_IndRJ_Lh", "Sis_DepRJ_Pred", "Sis_DepRJ_Lh",
                           "Anc_IndRJ_00", "Anc_IndRJ_01", "Anc_IndRJ_10", "Anc_IndRJ_11",
                           "Anc_DepRJ_00", "Anc_DepRJ_01", "Anc_DepRJ_10", "Anc_DepRJ_11",
                           "Anc_IndRJ_Pred", "Anc_IndRJ_Lh", "Anc_DepRJ_Pred", "Anc_DepRJ_Lh",
                           "True_IsMax_Sis_IndRJ", "Max_Sis_IndRJ_Lh", "True_IsMax_Sis_DepRJ", "Max_Sis_DepRJ_Lh",
                           "True_IsMax_Anc_IndRJ", "Max_Anc_IndRJ_Lh", "True_IsMax_Anc_DepRJ", "Max_Anc_DepRJ_Lh",
                           "TBL", "Tree_Height", "Std_TBL")
} else if (RJmodel == "BOTH") {
  AncConMat <- matrix(data=NA, nrow = nrow(results), ncol = 74)
  colnames(AncConMat) <- c("Trial_#", "Unknown_Taxon", "True_A", "True_B", "True_4States", "Avg_Sis_A", "Avg_Sis_B",
                           "Sis_IndMCMC_00", "Sis_IndMCMC_01", "Sis_IndMCMC_10", "Sis_IndMCMC_11",
                           "Sis_DepMCMC_00", "Sis_DepMCMC_01", "Sis_DepMCMC_10", "Sis_DepMCMC_11",
                           "Sis_IndRJ_00", "Sis_IndRJ_01", "Sis_IndRJ_10", "Sis_IndRJ_11",
                           "Sis_DepRJ_00", "Sis_DepRJ_01", "Sis_DepRJ_10", "Sis_DepRJ_11",
                           "Sis_IndMCMC_Pred", "Sis_IndMCMC_Lh", "Sis_DepMCMC_Pred", "Sis_DepMCMC_Lh",
                           "Sis_IndRJ_Pred", "Sis_IndRJ_Lh", "Sis_DepRJ_Pred", "Sis_DepRJ_Lh",
                           "Anc_IndMCMC_00", "Anc_IndMCMC_01", "Anc_IndMCMC_10", "Anc_IndMCMC_11",
                           "Anc_DepMCMC_00", "Anc_DepMCMC_01", "Anc_DepMCMC_10", "Anc_DepMCMC_11",
                           "Anc_IndRJ_00", "Anc_IndRJ_01", "Anc_IndRJ_10", "Anc_IndRJ_11",
                           "Anc_DepRJ_00", "Anc_DepRJ_01", "Anc_DepRJ_10", "Anc_DepRJ_11",
                           "Anc_IndMCMC_Pred", "Anc_IndMCMC_Lh", "Anc_DepMCMC_Pred", "Anc_DepMCMC_Lh",
                           "Anc_IndRJ_Pred", "Anc_IndRJ_Lh", "Anc_DepRJ_Pred", "Anc_DepRJ_Lh",
                           "True_IsMax_Sis_IndMCMC", "Max_Sis_IndMCMC_Lh", "True_IsMax_Sis_DepMCMC", "Max_Sis_DepMCMC_Lh",
                           "True_IsMax_Sis_IndRJ", "Max_Sis_IndRJ_Lh", "True_IsMax_Sis_DepRJ", "Max_Sis_DepRJ_Lh",
                           "True_IsMax_Anc_IndMCMC", "Max_Anc_IndMCMC_Lh", "True_IsMax_Anc_DepMCMC", "Max_Anc_DepMCMC_Lh",
                           "True_IsMax_Anc_IndRJ", "Max_Anc_IndRJ_Lh", "True_IsMax_Anc_DepRJ", "Max_Anc_DepRJ_Lh",
                           "TBL", "Tree_Height", "Std_TBL")
}



# Start the for-loop that will run each prediction and fill it in
for (i in 1:trial_amount) {
  # First, call rTaxon by subsetting the rows where column 1 matches i
  matching_rows <- which(results[, "Trial"] == i)
  rTaxon <- results[matching_rows, "rTaxon"]

  # Call the data file for this iteration
  data_name <- paste0(data_path, ".", i, ".Full_data.txt")
  data <- read.table(data_name, skip = 1, sep = "\t")
  data[, 4] <- match(data[,4], c(0,1,10,11))


  # Pull the tree to read the terminal branch length
  treename <- paste0("Trees/Full_tree.", i, ".tre")
  tree_complete <- read.nexus(treename)

  # Loop through every taxon for which we are predicting
  for (j in 1:length(rTaxon)) {
    # This collects the numerical average value of the sister taxon's trait information
    AncConMat[matching_rows[j], "Avg_Sis_A"] <- sisterData(tree_complete, rTaxon[j], data, "A")
    AncConMat[matching_rows[j], "Avg_Sis_B"] <- sisterData(tree_complete, rTaxon[j], data, "B")

    # Reset the sister node number and find the other important nodes
    sisters <- getSisters(tree_complete, rTaxon[j], mode = "number")
    sister_tips <- getAllDescendantTips(tree_complete, sisters)
    grandparent <- getParent(tree_complete, getMRCA(tree_complete, c(rTaxon, sister_tips)))
    cousins <- getAllDescendantTips(tree_complete, grandparent)
    cousins <- cousins[cousins != rTaxon[j]]

    # Find the node of interest in a tree where you don't know the unknown's information
    edited_tree <- drop.tip(tree_complete, rTaxon[j])
    edited_data <- data[-rTaxon,]
    edited_cousins <- match(cousins, edited_tree$tip.label)
    grandparent_node_asr <- getMRCA(edited_tree, edited_cousins)



    ### Ancestral State Reconstructin with both traits
    if (as.numeric(sisters) > Ntip(tree_complete) && length(unique(data[, 4])) > 1) {
      # Adjust node number to row number of the output
      sister_asr <- getMRCA(edited_tree, sister_tips)
      sister_asr <- sister_asr - Ntip(edited_tree)

      ### First pull the Independent MCMC model, and perform ML Ancestral State Reconstruction on the nearest sister taxon
      if (RJmodel == "MCMC" || RJmodel == "BOTH") {
        IndMCMC_qmat <- estimated_qmat(type, i, IndependentRates = TRUE, ConstantRates = TRUE, RJ = FALSE, Predict = FALSE, Random = TRUE)
        ace_IndMCMC <- asr_mk_model(tree = edited_tree, tip_states = edited_data[, 4], Nstates = 4, transition_matrix = IndMCMC_qmat, reroot = FALSE, root_prior = "flat")
        # Pull the likelihood
        Sis_IndMCMC_Lhs <- round(ace_IndMCMC$ancestral_likelihoods[sister_asr, ], digits = 6)
        # Store them in the matrix
        AncConMat[matching_rows[j], c("Sis_IndMCMC_00", "Sis_IndMCMC_01", "Sis_IndMCMC_10", "Sis_IndMCMC_11")] <- Sis_IndMCMC_Lhs

      ### Next Dependent MCMC models
        DepMCMC_qmat <- estimated_qmat(type, i, IndependentRates = FALSE, ConstantRates = T, RJ = FALSE, Predict = FALSE, Random = TRUE)
        ace_DepMCMC <- asr_mk_model(tree = edited_tree, tip_states = edited_data[, 4], Nstates = 4, transition_matrix = DepMCMC_qmat, reroot = FALSE, root_prior = "flat")
        # Pull the likelihood
        Sis_DepMCMC_Lhs <- round(ace_DepMCMC$ancestral_likelihoods[sister_asr, ], digits = 6)
        # Store them in the matrix
        AncConMat[matching_rows[j], c("Sis_DepMCMC_00", "Sis_DepMCMC_01", "Sis_DepMCMC_10", "Sis_DepMCMC_11")] <- Sis_DepMCMC_Lhs
      }


      ### Next Ind RJ models
      if (RJmodel == "RJMCMC" || RJmodel == "BOTH") {
        IndRJ_qmat <- estimated_qmat(type, i, IndependentRates = TRUE, ConstantRates = T, RJ = T, Predict = FALSE, Random = TRUE)
        ace_IndRJ <- asr_mk_model(tree = edited_tree, tip_states = edited_data[, 4], Nstates = 4, transition_matrix = IndRJ_qmat, reroot = FALSE, root_prior = "flat")
        # Pull the likelihood
        Sis_IndRJ_Lhs <- round(ace_IndRJ$ancestral_likelihoods[sister_asr, ], digits = 6)
        # Store them in the matrix
        AncConMat[matching_rows[j], c("Sis_IndRJ_00", "Sis_IndRJ_01", "Sis_IndRJ_10", "Sis_IndRJ_11")] <- Sis_IndRJ_Lhs

      ### Finally, the Dep RJ models
        DepRJ_qmat <- estimated_qmat(type, i, IndependentRates = TRUE, ConstantRates = T, RJ = FALSE, Predict = FALSE, Random = TRUE)
        ace_DepRJ <- asr_mk_model(tree = edited_tree, tip_states = edited_data[, 4], Nstates = 4, transition_matrix = DepRJ_qmat, reroot = FALSE, root_prior = "flat")
        # Pull the likelihood
        Sis_DepRJ_Lhs <- round(ace_DepRJ$ancestral_likelihoods[sister_asr, ], digits = 6)
        # Store them in the matrix
        AncConMat[matching_rows[j], c("Sis_DepRJ_00", "Sis_DepRJ_01", "Sis_DepRJ_10", "Sis_DepRJ_11")] <- Sis_DepRJ_Lhs
      }


      ### Line accounting for non-variation in the tree or single tip prediction
    } else {

      # In case of non-variation in the tree, make this adjustment for sister-node reconstruction
      # This only impacts cases of no trait variation across the tree when the sister node is internal
      # Although wonky, it covers edge cases without impacting sister-leaf reconstruction
      if (as.numeric(sisters) > Ntip(tree_complete)) {sisters <- sisters - Ntip(tree_complete)}

      # First cover the MCMC models
      if (RJmodel == "MCMC" || RJmodel == "BOTH") {
        Sis_IndMCMC_Lhs <- NULL
        Sis_IndMCMC_Lhs[1] <- Sis_IndMCMC_Lhs[2] <- Sis_IndMCMC_Lhs[3] <- Sis_IndMCMC_Lhs[4] <- 0
        Sis_IndMCMC_Lhs[ data[sisters,4] ] <- 1
        AncConMat[matching_rows[j], c("Sis_IndMCMC_00", "Sis_IndMCMC_01", "Sis_IndMCMC_10", "Sis_IndMCMC_11")] <- Sis_IndMCMC_Lhs

        Sis_DepMCMC_Lhs <- NULL
        Sis_DepMCMC_Lhs[1] <- Sis_DepMCMC_Lhs[2] <- Sis_DepMCMC_Lhs[3] <- Sis_DepMCMC_Lhs[4] <- 0
        Sis_DepMCMC_Lhs[ data[sisters,4] ] <- 1
        AncConMat[matching_rows[j], c("Sis_DepMCMC_00", "Sis_DepMCMC_01", "Sis_DepMCMC_10", "Sis_DepMCMC_11")] <- Sis_DepMCMC_Lhs
      }

      # Next cover the RJMCMC models
      if (RJmodel == "RJMCMC" || RJmodel == "BOTH") {
        Sis_IndRJ_Lhs <- NULL
        Sis_IndRJ_Lhs[1] <- Sis_IndRJ_Lhs[2] <- Sis_IndRJ_Lhs[3] <- Sis_IndRJ_Lhs[4] <- 0
        Sis_IndRJ_Lhs[ data[sisters,4] ] <- 1
        AncConMat[matching_rows[j], c("Sis_IndRJ_00", "Sis_IndRJ_01", "Sis_IndRJ_10", "Sis_IndRJ_11")] <- Sis_IndRJ_Lhs

        Sis_DepRJ_Lhs <- NULL
        Sis_DepRJ_Lhs[1] <- Sis_DepRJ_Lhs[2] <- Sis_DepRJ_Lhs[3] <- Sis_DepRJ_Lhs[4] <- 0
        Sis_DepRJ_Lhs[ data[sisters,4] ] <- 1
        AncConMat[matching_rows[j], c("Sis_DepRJ_00", "Sis_DepRJ_01", "Sis_DepRJ_10", "Sis_DepRJ_11")] <- Sis_DepRJ_Lhs
      }
    }


    ### Now pull the information for the ancestral node (the grandpa)
    # First, adjust to match the internal node number
    if (!is.null(grandparent)) {
      grandparent_node_asr <- grandparent_node_asr - Ntip(edited_tree)
    }


    ### Perform ML Ancestral State Reconstruction on the parent node of the Phylogenetic Bracket
    # First the MCMC models
    if ((RJmodel == "MCMC" || RJmodel == "BOTH") && !is.null(grandparent)) {
      IndMCMC_qmat <- estimated_qmat(type, i, IndependentRates = TRUE, ConstantRates = TRUE, RJ = FALSE, Predict = FALSE, Random = TRUE)
      Anc_ace_IndMCMC <- asr_mk_model(tree = edited_tree, tip_states = edited_data[, 4], Nstates = 4, transition_matrix = IndMCMC_qmat, reroot = FALSE, root_prior = "flat")
      # Pull the likelihood
      Anc_IndMCMC_Lhs <- round(Anc_ace_IndMCMC$ancestral_likelihoods[grandparent_node_asr, ], digits = 6)
      # Store them in the matrix
      AncConMat[matching_rows[j], c("Anc_IndMCMC_00", "Anc_IndMCMC_01", "Anc_IndMCMC_10", "Anc_IndMCMC_11")] <- Anc_IndMCMC_Lhs

      ### Next Dependent MCMC models
      DepMCMC_qmat <- estimated_qmat(type, i, IndependentRates = FALSE, ConstantRates = T, RJ = FALSE, Predict = FALSE, Random = TRUE)
      Anc_ace_DepMCMC <- asr_mk_model(tree = edited_tree, tip_states = edited_data[, 4], Nstates = 4, transition_matrix = DepMCMC_qmat, reroot = FALSE, root_prior = "flat")
      # Pull the likelihood
      Anc_DepMCMC_Lhs <- round(Anc_ace_DepMCMC$ancestral_likelihoods[grandparent_node_asr, ], digits = 6)
      # Store them in the matrix
      AncConMat[matching_rows[j], c("Anc_DepMCMC_00", "Anc_DepMCMC_01", "Anc_DepMCMC_10", "Anc_DepMCMC_11")] <- Anc_DepMCMC_Lhs
    }

    ### Next RJ models, starting with the Ind model
    if ((RJmodel == "RJMCMC" || RJmodel == "BOTH") && !is.null(grandparent)) {
      IndRJ_qmat <- estimated_qmat(type, i, IndependentRates = TRUE, ConstantRates = T, RJ = T, Predict = FALSE, Random = TRUE)
      Anc_ace_IndRJ <- asr_mk_model(tree = edited_tree, tip_states = edited_data[, 4], Nstates = 4, transition_matrix = IndRJ_qmat, reroot = FALSE, root_prior = "flat")
      # Pull the likelihood
      Anc_IndRJ_Lhs <- round(Anc_ace_IndRJ$ancestral_likelihoods[grandparent_node_asr, ], digits = 6)
      # Store them in the matrix
      AncConMat[matching_rows[j], c("Anc_IndRJ_00", "Anc_IndRJ_01", "Anc_IndRJ_10", "Anc_IndRJ_11")] <- Anc_IndRJ_Lhs

      ### Finally, the Dep RJ models
      DepRJ_qmat <- estimated_qmat(type, i, IndependentRates = F, ConstantRates = T, RJ = T, Predict = FALSE, Random = TRUE)
      Anc_ace_DepRJ <- asr_mk_model(tree = edited_tree, tip_states = edited_data[, 4], Nstates = 4, transition_matrix = DepRJ_qmat, reroot = FALSE, root_prior = "flat")
      # Pull the likelihood
      Anc_DepRJ_Lhs <- round(Anc_ace_DepRJ$ancestral_likelihoods[grandparent_node_asr, ], digits = 6)
      # Store them in the matrix
      AncConMat[matching_rows[j], c("Anc_DepRJ_00", "Anc_DepRJ_01", "Anc_DepRJ_10", "Anc_DepRJ_11")] <- Anc_DepRJ_Lhs
    }


    ### Adjust the code if the parent node is the root
    # First the MCMC models
    if ((RJmodel == "MCMC" || RJmodel == "BOTH") && is.null(grandparent)) {
      # Reset the likelihoods to clarify this trial does not have a valid "grandparent"
      Anc_IndMCMC_Lhs <- c(NA, NA, NA, NA)
      # Store them in the matrix
      AncConMat[matching_rows[j], c("Anc_IndMCMC_00", "Anc_IndMCMC_01", "Anc_IndMCMC_10", "Anc_IndMCMC_11")] <- Anc_IndMCMC_Lhs

      ### Next Dependent MCMC models
      # Reset the likelihoods to clarify this trial does not have a valid "grandparent"
      Anc_DepMCMC_Lhs <- c(NA, NA, NA, NA)
      # Store them in the matrix
      AncConMat[matching_rows[j], c("Anc_DepMCMC_00", "Anc_DepMCMC_01", "Anc_DepMCMC_10", "Anc_DepMCMC_11")] <- Anc_DepMCMC_Lhs
    }

    ### Next RJ models, starting with the Ind model
    if ((RJmodel == "RJMCMC" || RJmodel == "BOTH") && is.null(grandparent)) {
      # Reset the likelihoods to clarify this trial does not have a valid "grandparent"
      Anc_IndRJ_Lhs <- c(NA, NA, NA, NA)
      # Store them in the matrix
      AncConMat[matching_rows[j], c("Anc_IndRJ_00", "Anc_IndRJ_01", "Anc_IndRJ_10", "Anc_IndRJ_11")] <- Anc_IndRJ_Lhs

      ### Finally, the Dep RJ models
      # Reset the likelihoods to clarify this trial does not have a valid "grandparent"
      Anc_DepRJ_Lhs <- c(NA, NA, NA, NA)
      # Store them in the matrix
      AncConMat[matching_rows[j], c("Anc_DepRJ_00", "Anc_DepRJ_01", "Anc_DepRJ_10", "Anc_DepRJ_11")] <- Anc_DepRJ_Lhs
    }


    ### Pull the predicted states and likelihoods for the tests performed
    if ((RJmodel == "MCMC" || RJmodel == "BOTH") && !is.null(grandparent)) {
      # First for the sister
      AncConMat[matching_rows[j], "Sis_IndMCMC_Pred"] <- which.max(Sis_IndMCMC_Lhs)
      AncConMat[matching_rows[j], "Sis_IndMCMC_Lh"] <- max(Sis_IndMCMC_Lhs)


      AncConMat[matching_rows[j], "Sis_DepMCMC_Pred"] <- which.max(Sis_DepMCMC_Lhs)
      AncConMat[matching_rows[j], "Sis_DepMCMC_Lh"] <- max(Sis_DepMCMC_Lhs)

      # Now for the grandparent
      AncConMat[matching_rows[j], "Anc_IndMCMC_Pred"] <- which.max(Anc_IndMCMC_Lhs)
      AncConMat[matching_rows[j], "Anc_IndMCMC_Lh"] <- max(Anc_IndMCMC_Lhs)

      AncConMat[matching_rows[j], "Anc_DepMCMC_Pred"] <- which.max(Anc_DepMCMC_Lhs)
      AncConMat[matching_rows[j], "Anc_DepMCMC_Lh"] <- max(Anc_DepMCMC_Lhs)

      # Separate out the likelihood of the "true" state of rTaxon
      results[matching_rows[j], "TrueState_Sister_IndMCMC_Lh"] <- Sis_IndMCMC_Lhs[results[matching_rows[j], "True_4States"]]
      results[matching_rows[j], "TrueState_Sister_DepMCMC_Lh"] <- Sis_DepMCMC_Lhs[results[matching_rows[j], "True_4States"]]

      results[matching_rows[j], "TrueState_Ancestor_IndMCMC_Lh"] <- Anc_IndMCMC_Lhs[results[matching_rows[j], "True_4States"]]
      results[matching_rows[j], "TrueState_Ancestor_DepMCMC_Lh"] <- Anc_DepMCMC_Lhs[results[matching_rows[j], "True_4States"]]
    }

    if ((RJmodel == "RJMCMC" || RJmodel == "BOTH") && !is.null(grandparent)) {
      # First for the sister
      AncConMat[matching_rows[j], "Sis_IndRJ_Pred"] <- which.max(Sis_IndRJ_Lhs)
      AncConMat[matching_rows[j], "Sis_IndRJ_Lh"] <- max(Sis_IndRJ_Lhs)

      AncConMat[matching_rows[j], "Sis_DepRJ_Pred"] <- which.max(Sis_DepRJ_Lhs)
      AncConMat[matching_rows[j], "Sis_DepRJ_Lh"] <- max(Sis_DepRJ_Lhs)

      # Now for the grandparent
      AncConMat[matching_rows[j], "Anc_IndRJ_Pred"] <- which.max(Anc_IndRJ_Lhs)
      AncConMat[matching_rows[j], "Anc_IndRJ_Lh"] <- max(Anc_IndRJ_Lhs)

      AncConMat[matching_rows[j], "Anc_DepRJ_Pred"] <- which.max(Anc_DepRJ_Lhs)
      AncConMat[matching_rows[j], "Anc_DepRJ_Lh"] <- max(Anc_DepRJ_Lhs)

      # Separate out the likelihood of the "true" state of rTaxon
      results[matching_rows[j], "TrueState_Sister_IndRJ_Lh"] <- Sis_IndMCMC_Lhs[results[matching_rows[j], "True_4States"]]
      results[matching_rows[j], "TrueState_Sister_DepRJ_Lh"] <- Sis_DepMCMC_Lhs[results[matching_rows[j], "True_4States"]]

      results[matching_rows[j], "TrueState_Ancestor_IndRJ_Lh"] <- Anc_IndMCMC_Lhs[results[matching_rows[j], "True_4States"]]
      results[matching_rows[j], "TrueState_Ancestor_DepRJ_Lh"] <- Anc_DepMCMC_Lhs[results[matching_rows[j], "True_4States"]]
    }


    ### A fix for when there is no grandparent node (parent is the root)
    # These are again separated by tests run
    if ((RJmodel == "MCMC" || RJmodel == "BOTH") && is.null(grandparent)) {
      # Pull the predicted states and likelihoods for the sister
      AncConMat[matching_rows[j], "Sis_IndMCMC_Pred"] <- which.max(Sis_IndMCMC_Lhs)
      AncConMat[matching_rows[j], "Sis_IndMCMC_Lh"] <- max(Sis_IndMCMC_Lhs)

      AncConMat[matching_rows[j], "Sis_DepMCMC_Pred"] <- which.max(Sis_DepMCMC_Lhs)
      AncConMat[matching_rows[j], "Sis_DepMCMC_Lh"] <- max(Sis_DepMCMC_Lhs)

      # Separate out the likelihood of the "true" state of rTaxon
      results[matching_rows[j], "TrueState_Sister_IndMCMC_Lh"] <- Sis_IndMCMC_Lhs[results[matching_rows[j], "True_4States"]]
      results[matching_rows[j], "TrueState_Sister_DepMCMC_Lh"] <- Sis_DepMCMC_Lhs[results[matching_rows[j], "True_4States"]]

      # However, none of the "ancestor" state reconstruction can be done without a known ancestor.
      AncConMat[matching_rows[j], "Anc_IndMCMC_Pred"] <- NA
      AncConMat[matching_rows[j], "Anc_IndMCMC_Lh"] <- NA

      AncConMat[matching_rows[j], "Anc_DepMCMC_Pred"] <- NA
      AncConMat[matching_rows[j], "Anc_DepMCMC_Lh"] <- NA

      results[matching_rows[j], "TrueState_Ancestor_IndMCMC_Lh"] <- NA
      results[matching_rows[j], "TrueState_Ancestor_DepMCMC_Lh"] <- NA
    }

    if ((RJmodel == "RJMCMC" || RJmodel == "BOTH") && is.null(grandparent)) {
      # Pull the predicted states and likelihoods for the sister
      AncConMat[matching_rows[j], "Sis_IndRJ_Pred"] <- which.max(Sis_IndRJ_Lhs)
      AncConMat[matching_rows[j], "Sis_IndRJ_Lh"] <- max(Sis_IndRJ_Lhs)

      AncConMat[matching_rows[j], "Sis_DepRJ_Pred"] <- which.max(Sis_DepRJ_Lhs)
      AncConMat[matching_rows[j], "Sis_DepRJ_Lh"] <- max(Sis_DepRJ_Lhs)

      # Separate out the likelihood of the "true" state of rTaxon
      results[matching_rows[j], "TrueState_Sister_IndRJ_Lh"] <- Sis_IndRJ_Lhs[results[matching_rows[j], "True_4States"]]
      results[matching_rows[j], "TrueState_Sister_DepRJ_Lh"] <- Sis_DepRJ_Lhs[results[matching_rows[j], "True_4States"]]

      # However, none of the "ancestor" state reconstruction can be done without a known ancestor.
      AncConMat[matching_rows[j], "Anc_IndRJ_Pred"] <- NA
      AncConMat[matching_rows[j], "Anc_IndRJ_Lh"] <- NA

      AncConMat[matching_rows[j], "Anc_DepRJ_Pred"] <- NA
      AncConMat[matching_rows[j], "Anc_DepRJ_Lh"] <- NA

      results[matching_rows[j], "TrueState_Ancestor_IndRJ_Lh"] <- NA
      results[matching_rows[j], "TrueState_Ancestor_DepRJ_Lh"] <- NA
    }


    ### These lines add to the results whether or not the ACR would agree with our predictions or the true, unknown state
    # First the MCMC models
    if (RJmodel == "MCMC" || RJmodel == "BOTH") {
      # First determine if the maximum likelihood from the ASR matches the unknown taxon's character state
      AncConMat[matching_rows[j], "True_IsMax_Sis_IndMCMC"] <- ifelse(Sis_IndMCMC_Lhs[results[matching_rows[j], "True_4States"]] == max(round(Sis_IndMCMC_Lhs, digits = 6)), 1, 0)
      AncConMat[matching_rows[j], "True_IsMax_Sis_DepMCMC"] <- ifelse(Sis_DepMCMC_Lhs[results[matching_rows[j], "True_4States"]] == max(round(Sis_DepMCMC_Lhs, digits = 6)), 1, 0)

      AncConMat[matching_rows[j], "True_IsMax_Anc_IndMCMC"] <- ifelse(!is.null(grandparent), ifelse(Anc_IndMCMC_Lhs[results[matching_rows[j], "True_4States"]] == max(round(Anc_IndMCMC_Lhs, digits = 6)), 1, 0), NA)
      AncConMat[matching_rows[j], "True_IsMax_Anc_DepMCMC"] <- ifelse(!is.null(grandparent), ifelse(Anc_DepMCMC_Lhs[results[matching_rows[j], "True_4States"]] == max(round(Anc_DepMCMC_Lhs, digits = 6)), 1, 0), NA)

      # Now save whichever is the maximum likelihood
      AncConMat[matching_rows[j], "Max_Sis_IndMCMC_Lh"] <- max(round(Sis_IndMCMC_Lhs, digits = 6))
      AncConMat[matching_rows[j], "Max_Sis_DepMCMC_Lh"] <- max(round(Sis_DepMCMC_Lhs, digits = 6))

      AncConMat[matching_rows[j], "Max_Anc_IndMCMC_Lh"] <- ifelse(!is.null(grandparent), max(round(Anc_IndMCMC_Lhs, digits = 6)), NA)
      AncConMat[matching_rows[j], "Max_Anc_DepMCMC_Lh"] <- ifelse(!is.null(grandparent), max(round(Anc_DepMCMC_Lhs, digits = 6)), NA)
    }

    # Next RJ models, starting with the Ind model
    if (RJmodel == "RJMCMC" || RJmodel == "BOTH") {
      # First determine if the maximum likelihood from the ASR matches the unknown taxon's character state
      AncConMat[matching_rows[j], "True_IsMax_Sis_IndRJ"] <- ifelse(Sis_IndRJ_Lhs[results[matching_rows[j], "True_4States"]] == max(round(Sis_IndRJ_Lhs, digits = 6)), 1, 0)
      AncConMat[matching_rows[j], "True_IsMax_Sis_DepRJ"] <- ifelse(Sis_DepRJ_Lhs[results[matching_rows[j], "True_4States"]] == max(round(Sis_DepRJ_Lhs, digits = 6)), 1, 0)

      AncConMat[matching_rows[j], "True_IsMax_Anc_IndRJ"] <- ifelse(!is.null(grandparent), ifelse(Anc_IndRJ_Lhs[results[matching_rows[j], "True_4States"]] == max(round(Anc_IndRJ_Lhs, digits = 6)), 1, 0), NA)
      AncConMat[matching_rows[j], "True_IsMax_Anc_DepRJ"] <- ifelse(!is.null(grandparent), ifelse(Anc_DepRJ_Lhs[results[matching_rows[j], "True_4States"]] == max(round(Anc_DepRJ_Lhs, digits = 6)), 1, 0), NA)

      # Now save whichever is the maximum likelihood
      AncConMat[matching_rows[j], "Max_Sis_IndRJ_Lh"] <- max(round(Sis_IndRJ_Lhs, digits = 6))
      AncConMat[matching_rows[j], "Max_Sis_DepRJ_Lh"] <- max(round(Sis_DepRJ_Lhs, digits = 6))

      AncConMat[matching_rows[j], "Max_Anc_IndRJ_Lh"] <- ifelse(!is.null(grandparent), max(round(Anc_IndRJ_Lhs, digits = 6)), NA)
      AncConMat[matching_rows[j], "Max_Anc_DepRJ_Lh"] <- ifelse(!is.null(grandparent), max(round(Anc_DepRJ_Lhs, digits = 6)), NA)
    }


    ### Pull Terminal Branch Length and Tree Height
    AncConMat[matching_rows[j], "TBL"] <- TBL <- unknown_branch_length(rTaxon[j], tree_complete)
    AncConMat[matching_rows[j], "Tree_Height"] <- tree_height <- max(node.depth.edgelength(tree_complete))
    AncConMat[matching_rows[j], "Std_TBL"] <- TBL/tree_height

    ### Finally, add the trial information to this table
    AncConMat[matching_rows[j], "Trial_#"] <- i
    AncConMat[matching_rows[j], "Unknown_Taxon"] <- rTaxon[j]
    AncConMat[matching_rows[j], "True_A"] <- data[rTaxon[j], 2]
    AncConMat[matching_rows[j], "True_B"] <- data[rTaxon[j], 3]
    AncConMat[matching_rows[j], "True_4States"] <- data[rTaxon[j], 4]
  }
}
### Save this table
write.table(AncConMat, file = matrix_name, quote = F, row.names = F, col.names = T, sep = "\t")


### Add the important information into the main results folder, starting with Averaged Sister Trait Values
results[, "Avg_Sister_A"] <- AncConMat[, "Avg_Sis_A"]
results[, "Avg_Sister_B"] <- AncConMat[, "Avg_Sis_B"]

### Get the important information from the ancestral reconstructions, and add it to the main results file
if (RJmodel == "MCMC" || RJmodel == "BOTH") {
  results[, "Sister_IndMCMC_MaxLh"] <- AncConMat[, "Sis_IndMCMC_Lh"]
  results[, "Sister_IndMCMC_Prediction"] <- AncConMat[, "Sis_IndMCMC_Pred"]
  results[, "Sister_DepMCMC_MaxLh"] <- AncConMat[, "Sis_DepMCMC_Lh"]
  results[, "Sister_DepMCMC_Prediction"] <- AncConMat[, "Sis_DepMCMC_Pred"]

  results[, "Ancestor_IndMCMC_MaxLh"] <- AncConMat[, "Anc_IndMCMC_Lh"]
  results[, "Ancestor_IndMCMC_Prediction"] <- AncConMat[, "Anc_IndMCMC_Pred"]
  results[, "Ancestor_DepMCMC_MaxLh"] <- AncConMat[, "Anc_DepMCMC_Lh"]
  results[, "Ancestor_DepMCMC_Prediction"] <- AncConMat[, "Anc_DepMCMC_Pred"]
}
if (RJmodel == "RJMCMC" || RJmodel == "BOTH") {
  results[, "Sister_IndRJ_MaxLh"] <- AncConMat[, "Sis_IndRJ_Lh"]
  results[, "Sister_IndRJ_Prediction"] <- AncConMat[, "Sis_IndRJ_Pred"]
  results[, "Sister_DepRJ_MaxLh"] <- AncConMat[, "Sis_DepRJ_Lh"]
  results[, "Sister_DepRJ_Prediction"] <- AncConMat[, "Sis_DepRJ_Pred"]

  results[, "Ancestor_IndRJ_MaxLh"] <- AncConMat[, "Anc_IndRJ_Lh"]
  results[, "Ancestor_IndRJ_Prediction"] <- AncConMat[, "Anc_IndRJ_Pred"]
  results[, "Ancestor_DepRJ_MaxLh"] <- AncConMat[, "Anc_DepRJ_Lh"]
  results[, "Ancestor_DepRJ_Prediction"] <- AncConMat[, "Anc_DepRJ_Pred"]
}

# Now save the changes to the results matrices
write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t")

# Finally, send a progress signal
print(paste0("Finished with the ancestral reconstructions for Random Rates"))