#!/bin/bash
## Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
## All other lines are read by the shell.

################################################################################
## Written by Jack Prall, Liam Feigin, Emerald Bender and Chris Organ
## October 2024

################################################################################
## Running the script below will test multiple phylogenetic and non-phylogenetic prediction models. 
## Calls the files creates by 'TestModels'
## Runs BayesTraits phylogenetic analyses
## Calls and runs 'TestModels' and 'CompileResults' in R
## Written for the Tempest Research Cluster at Montana State University

## Outputs will include any files made by BayesTraits and three Results files
## Results are stored within a 'Results' folder found within the folder for each matrix.

################################################################################

##SBATCH --account=priority-chrisorgan   # The class account
##SBATCH --partition=priority               # queue partition to run the job in
##SBATCH --ntasks=25                         # number of mpi processes (1 for multi-threaded tasks)
##SBATCH --threads-per-core=2               # Turn off hyperthreading
##SBATCH --cpus-per-task=1
##SBATCH --mem=10G
##SBATCH --time=1-00:00:00                 # maximum job run time in days-hours:minutes:seconds
##SBATCH --job-name=JOB.BayesTraits
##SBATCH --output=BT-%j.out                 # standard output from job
##SBATCH --error=BT-%j.err                  # standard error from job
##SBATCH --mail-user=bender.er3@gmail.com
##SBATCH --mail-type=ALL
##SBATCH --no-requeue

## Simulation settings
## Simulation version
sim_version="V10"

## Set the number of iterations
num_iterations=100

## Set the number of taxa that are used to build the tree and in the prediction trials
##	  "Tree name" is the filename of the nexus file for the real reptile tree
##    "Pop size" is the number of sampled tips to test with, which must be equal to or less than tree_size
tree_name=rept_eggdanger_timetree_apr25_V2.nex
pop_size=50

## Specify the input data file
## trait 2 endangeredness, trait 1 oviparity
data_name=rept_eggdanger_timetree_apr25_V2.txt

## Use MCMC, RJMCMC, or test both
##RJmodel=MCMC
##RJmodel=RJMCMC
RJmodel=BOTH

## Define the types and trials arrays - unnecessary in real data set so replace with single name
types=("rept_eggdanger_100drop")

## Set the Rates- Unnecessary in real data tests
## low_rate=0.1
## medium_rate=0.25
## high_rate=0.5

## Set a vector of scale factors which will be used to create rate-heterogeneity when simulating data
##tree_scales=(1 1 1 1)
##tree_scales=(0.5 1 1 1 1 2 3)

## Set the degree of dependency (multiplier on rates)
##   and adjustment so that the sum of the rates remains constant
##depend_scale1=0.5
##depend_adj1=1.16667

##depend_scale2=0.1
##depend_adj2=1.3

## Set the settings for the prediction target
unknown_n=100
sample_type=Random
##sample_type=Clade

## First, take note of the current directory, so we can return to it later.
original_dir=$(pwd)

## Loop over each type and create their directories
for type in "${types[@]}"; do
  ## Create a directory and name it
  directoryname="${type}"
  mkdir -p "$directoryname"
done

## Load R and BayesTraits into Tempest
##module load R/4.3.2-gfbf-2023a
##module load BayesTraits/4.1.3-Linux

## Call the two scripts that define the necessary functions
Rscript DiscreteFunctions.Rept.${sim_version}.R
Rscript TestModels.Rept.${sim_version}.R

## Call the script that will run the whole simulation
Rscript SimInstructions.Rept.${sim_version}.R

## Function to run BayesTraits for both independent and dependent tests
run_rates() {
  local i=$1
  local trial_name=$2
  local RJmodel=$3

  ## Set the general settings
  trial_number="${trial_name}.${i}"
  treename="${trial_number}.edited_tree.tre"
  dataname="${trial_number}.edited_data.txt"
  
   echo "starting rates ${i}"


  # Run any MCMC tests you've set up
  if [[ "$RJmodel" == "MCMC" || "$RJmodel" == "BOTH" ]]; then
    ## Multistate MCMC rates
    ##settings="${trial_name}.Multistate.MCMC.${i}.Rates.In.txt"
    ##BayesTraitsV4 "${treename}" "${dataname}" <"${settings}"

    ## Independent MCMC rates
    settings="${trial_name}.Ind.MCMC.$i.Rates.In.txt"
    BayesTraitsV4 "${treename}" "${dataname}" <"${settings}"
	
	echo "done ind mcmc rates ${i}"

    ## Dependent MCMC rates
    settings="${trial_name}.Dep.MCMC.${i}.Rates.In.txt"
    BayesTraitsV4 "${treename}" "${dataname}" <"${settings}"
	echo "done dep mcmc rates ${i}"

  fi

  # Run any RJMCMC tests you've set up
  if [[ "$RJmodel" == "RJMCMC" || "$RJmodel" == "BOTH" ]]; then
    ## Multistate RJ rates
    ##settings="${trial_name}.Multistate.RJMCMC.${i}.Rates.In.txt"
    ##BayesTraitsV4 "${treename}" "${dataname}" <"${settings}"    
    
    ## Independent RJ rates
    settings="${trial_name}.Ind.RJMCMC.${i}.Rates.In.txt"
    BayesTraitsV4 "${treename}" "${dataname}" <"${settings}"
	echo "done indRJ mcmc rates ${i}"


    ## Dependent RJ rates
    settings="${trial_name}.Dep.RJMCMC.${i}.Rates.In.txt"
    BayesTraitsV4 "${treename}" "${dataname}" <"${settings}"
	echo "done depRJ mcmc rates ${i}"

  fi
}

## Function to run BayesTraits for both independent and dependent tests
run_prediction() {
  local i=$1
  local trial_name=$2
  local RJmodel=$3

  ## Set the general settings
  trial_number="${trial_name}.${i}"
  treename="${trial_number}.full_tree.tre"
  dataname="${trial_number}.predict_data.txt"
  
  echo "starting prediction ${i}"

  # Run any MCMC tests you've set up
  if [[ "$RJmodel" == "MCMC" || "$RJmodel" == "BOTH" ]]; then
    ## Multistate Multistate test
    ##settings="${trial_name}.Multistate.MCMC.${i}.Predict.In.txt"
    ##BayesTraitsV4 "${treename}" "${dataname}" <"${settings}"    
    
    ## Independent MCMC test
    settings="${trial_name}.Ind.MCMC.$i.Predict.In.txt"
    BayesTraitsV4 "${treename}" "${dataname}" <"${settings}"
	
	echo "done ind mcmc prediction ${i}"

    ## Dependent MCMC test
    settings="${trial_name}.Dep.MCMC.${i}.Predict.In.txt"
    BayesTraitsV4 "${treename}" "${dataname}" <"${settings}"
	
	echo "done dep mcmc prediction ${i}"

  fi

  # Run any RJMCMC tests you've set up
  if [[ "$RJmodel" == "RJMCMC" || "$RJmodel" == "BOTH" ]]; then
    ## Multistate Multistate test
    ##settings="${trial_name}.Multistate.RJMCMC.${i}.Predict.In.txt"
    ##BayesTraitsV4 "${treename}" "${dataname}" <"${settings}"

    ## Independent RJ rates
    settings="${trial_name}.Ind.RJMCMC.${i}.Predict.In.txt"
    BayesTraitsV4 "${treename}" "${dataname}" <"${settings}"
	
	echo "done indRJ mcmc prediction ${i}"

    ## Dependent RJ rates
    settings="${trial_name}.Dep.RJMCMC.${i}.Predict.In.txt"
    BayesTraitsV4 "${treename}" "${dataname}" <"${settings}"
	echo "done depRJ mcmc prediction ${i}"

	
  fi
}

export -f run_rates
export -f run_prediction


## Then we do all the phylogenetic stuff
## Loop over each type
##for type in "${types[@]}"; do
type=$types
  directoryname="${type}"
  cd "./$directoryname" || { echo "Failed to change to directory: $directoryname"; exit 1; }

  ## Set the trial name
  trial_name="${type}"

  ## Placeholder for additional operations with $trial_name if needed
  echo "Processing trial: $trial_name"

  ## Run the rate calculation jobs in parallel
 ##for i in $(seq 1 $num_iterations); do
 ##   run_rates $i $trial_name $RJmodel &
 ## done
  
  ## add --line-buffer to parallel to see output of BT
  seq $num_iterations | parallel -j15 run_rates {#} $trial_name $RJmodel
	##TEMPORARILY SKIP RUN_RATES 04/21/2025, add back everything with 3 ###
  ## Wait for all background jobs to complete
  wait

  ## Run the prediction jobs in parallel
  ##for i in $(seq 1 $num_iterations); do
  ##  run_prediction $i $trial_name $RJmodel &
  ##done
  
   seq $num_iterations | parallel -j15 run_prediction {#} $trial_name $RJmodel

  ## Wait for all background jobs to complete
  wait

  ## Change back to the original directory
  cd "$original_dir"

done

## This script will compile and summarize all our results
##Rscript CompileResults.${sim_version}.R
