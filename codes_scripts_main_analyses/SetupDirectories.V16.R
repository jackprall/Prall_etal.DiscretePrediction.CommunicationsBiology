# Into/This Script is for...



### First get the working directory and establish what we need from the shell script
original_dir <- getwd()

# Get simulation bash file and any useful functions for this script
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


### Create the directories
# First, the directory that will contain the tree
if (!dir.exists("Trees")) {dir.create("Trees")}

# Second, the directory that will contain the Results
if (!dir.exists("Results")) {
  dir.create("Results/ConstantRates/Single/AncestralStateReconstruction", recursive = TRUE)
  dir.create("Results/Random/AncestralStateReconstruction", recursive = TRUE)
  # Remember to create directories to separate by settings
  if (multiple_prediction == TRUE) {dir.create("Results/ConstantRates/Multiple/AncestralStateReconstruction", recursive = TRUE)}
  if (clade_prediction == TRUE) {dir.create("Results/ConstantRates/Clade/AncestralStateReconstruction", recursive = TRUE)}
  if (variable_rates == TRUE) {
    dir.create("Results/VariableRates/Single/AncestralStateReconstruction", recursive = TRUE)
    if (multiple_prediction == TRUE) {dir.create("Results/VariableRates/Multiple/AncestralStateReconstruction", recursive = TRUE)}
    if (clade_prediction == TRUE) {dir.create("Results/VariableRates/Clade/AncestralStateReconstruction", recursive = TRUE)}
  }
}

# Third, a directory to contain the Randomly generated data
if (!dir.exists("Random")) {
  dir.create("Random")
  # Next, create the directories within that will hold the data and output files
  dir.create("Random/Data")
  dir.create("Random/Single")
  # Check if we need additional directories
  if (multiple_prediction == TRUE) {dir.create("Random/Multiple")}
  if (clade_prediction == TRUE) {dir.create("Random/Clade")}
}

# Check if we are adding variable rates to this run
# If so, add those directories
if (variable_rates == TRUE && !dir.exists("VariableRates")) {
  dir.create("VariableRates")
  for (type in types) {
    # First create the main folders
    path <- paste0("VariableRates/", type)
    if (!dir.exists(path)) {dir.create(path)}
    # Then, the necessary directories
    dir.create(paste0(path, "/Data"))
    dir.create(paste0(path, "/Single"))
    # Check if we need additional directories
    if (multiple_prediction == TRUE) {dir.create(paste0(path, "/Multiple"))}
    if (clade_prediction == TRUE) {dir.create(paste0(path, "/Clade"))}

  }
}

# Finally, create the constant rates directories
if (!dir.exists("ConstantRates")) {
  dir.create("ConstantRates")
  for (type in types) {
    # First create the main folders
    path <- paste0("ConstantRates/", type)
    if (!dir.exists(path)) {dir.create(path)}
    # Then, the necessary directories
    dir.create(paste0(path, "/Data"))
    dir.create(paste0(path, "/Single"))
    # Check if we need additional directories
    if (multiple_prediction == TRUE) {dir.create(paste0(path, "/Multiple"))}
    if (clade_prediction == TRUE) {dir.create(paste0(path, "/Clade"))}

  }
}

print("All necessary directories have been created.")
