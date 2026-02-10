# Discrete Trait Prediction

This repository contains R scripts and other supplementary files, such as Excel result files, for "Predicting Discrete Traits in Evolving Systems".

## Contents

- 1_codes_scripts_main_analyses

    - Discrete_Simulation.*.sh
        - This shell script calls the others to perform the whole simulation.
        - This script also runs the phylogenetic analysis through the software BayesTraits.
        - This holds the primary settings for this simulation.

    - DiscreteFunctions.*.R
        - This R script details all sub-functions used in the next following scripts.

    - SetupDirectories.*.R
        - This script sets up the directories necessary for the simulation to run

    - TreeGeneration.*.R
        - This script generates the phylogenetic trees used in the simulation

    - FullDataGeneration.*.R
        - This script geneates the data for every trial in the simulation based on settings in the shell script

    - SampleSingleTaxa.*.R
        - This script samples one taxon from each data and records its data to be predicted later

    - SampleMultipleTaxa.*.R
        - This script only runs if multiple_prediction or clade_prediction are "true" in the shell script
        - It functions the same as SampleSingleTaxa.*.R except it samples a number of taxa equal to the unknown_size setting in the shell script

    - ResultsMatrixGeneration.*.R
        - This script generates the final results table and fills in the taxon and tree information

    - FilesForBayesTraits.*.R
        - This script generates the input files necessary to run BayesTraits based on the simulated data and settings

    - BBandNBPrediction.*.R
        - This script runs the Beta Binomial and Naive Bayes predictions
        - Results are stored in the main results tables in the Results folder

    - RunBayesTraits.*.R
        - This script performs the BayesTraits runs for this simulation
        - The type of BayesTraits runs are determined in the settings in the shell script

    - CompileBayesTraits.*.R
        - This script compiles the results in the log files from BayesTraits
        - It stores this information in the main results tables in the Results folder

    - AncestralStateReconstruction.*.R
        - This state performs Ancestral State Reconstruction on the sister and grandparent nodes of the unknown taxon using models created by BayesTraits
        - The full results are stored in the "AncestralStateReconstruction" subfolders within the Results folder
        - Summary stats are pulled from this table into the main results table

    - SummarizeResults.*.R
        - This script summarizes the results into more digestable summary tables
        - These are stored in the Results folder

    - RemoveSuperfluousFiles.*.R
        - This script removes any additional files created in the middle of the simulation to decrease the storage costs of the simulation long term
        - It removes instruction, schedule, model, and particular data files which are only used as inputs for BayesTraits

- 2_codes_scripts_others_stochastic_mapping

	- 00_archive
		- 
	
	- 01_scripts
		- 
	
	- 02_inputs
		- 
	
	- 03_outputs
		- 

- results

## Details...

### System and Program Requirements

- R 4.4.3
    - This can be found at the following URL: https://www.r-project.org/
- BayesTraits 5.0.1
    - This can be found at the following URL: https://www.evolution.reading.ac.uk/BayesTraitsV5.0.1/BayesTraitsV5.0.1.html
- 250 Gigabytes of extra storage space on the device you run the simulation on (for the full simulation performed by the authors, but this number drops with fewer additional runs).

These scripts were run using the HPC research cluster at Montana State University known as Tempest. Using the settings found in 'Discrete_Simulation.sh', the code will download any required R packages, run the full simulation, and summarize the results in a specified file.
More information on the Tempest research cluster can be found at https://www.montana.edu/uit/rci/tempest/

### Instructions

To repeat the study, download all scripts into a folder on the computer which is to hold and run the simulation. Next, put all R scripts into a subdirectory called "Scripts" within this main folder. Make sure to add the BayesTraitsV5 executable into this main folder (See "File Organization" below for more detail).
Next, open the command line of your computer and navigate to the folder containing these scripts. Then run 'Discrete_Simulation.sh' using a line that should look like "bash Discrete_Simulation.sh"
Alternatively, you can run each R script individually, but they are written to run in the order that they are found in the shell script.
This took two weeks to run using the research cluster at Montana State University using the full settings described below, but may take longer on a less powerful machine.

### Settings

In the shell script "Discrete_Simulation.*.sh", the first 40 lines of code establish the 18 settings for the simulation. They are described below with suggests of how to use them.

- sim_version
    - This represents the version of the simulation that the simulation will look for and use
    - This should be set to "V16" or to match the other scripts

- num_iterations
    - This setting represents the number of trials will be performed for each method of simulating data
    - Our simulation had 1000 trials for each method of simulating data
    - We recommend testing with 5-10 trials, as this has the largest impact on run time, but full simulations should be done with more
    - This can be any integer greater than 1

- pop_size
    - This setting represents the number of tips on the tree being used
    - Our simulation used 50, 100, or 500 tips across three runs of the simulation
    - This also has a large impact on run time
    - This can be any integer greater than 1

- RJmodel
    - This setting determines whether we run the Reversible Jump algorithm in BayesTraits or not
    - You can use a standard MCMC algorithm, the RJ-MCMC, or test both.
    - Our simulation tested both MCMC and RJ-MCMC algorithms
    - Options for this setting are...
        - MCMC
        - RJMCMC
        - BOTH

- types
    - This setting represents the types of evolutionary scenarios used to simulate data
    - Our simulation used the full list of coded evolutionary scenarios which are...
        - "ER.L", "ER.M", "ER.H", "DR.LM", "DR.LH", "DR.MH", "DEP1.L", "DEP1.M", "DEP1.H", "DEP2.L", "DEP2.M", and "DEP2.H"
        - Models
            - "ER"      - This denotes an equal rates model
            - "DR"      - This denotes a direction rates model
            - "DEP"     - This denotes a dependent rates model
            - "Random"  - This denotes a model where data is simulated randomly which is automatically run and is not show in the settings
        - Rates
            - "L"   - This denotes that the base rate is LOW for an equal or dependent model.
            - "M"   - This denotes that the base rate is MEDIUM for an equal or dependent model.
            - "H"   - This denotes that the base rate is HIGH for an equal or dependent model.
            - "LM"  - This denotes that the beta rate is LOW and the alpha rate is MEDIUM for a directional model.
            - "LH"  - This denotes that the beta rate is LOW and the alpha rate is HIGH for a directional model.
            - "MH"  - This denotes that the beta rate is MEDIUM and the alpha rate is HIGH for a directional model.
    - These matrices are created relative to three rate parameters and two dependent scalars defined below

- low_rate
- medium_rate
- high_rate
    - These settings represent the rate parameters, L, M, and H
    - Our simulation set the rates as...
        - L = 0.1
        - M = 0.25
        - H = 0.5
    - These can be set to any number greater than 0
    - They do not have to follow L < M < H, but it is recommended to help the user interpret results

- variable_rates
    - This setting determines if you will simulate data using variable rates of evolution in addition to constant rates of evolution
    - This will roughly double the simulation time if set to true
    - This must be set to...
        - true
        - false

- tree_scales
    - This provides a list of values, one of which will be randomly applied to every branch length before simulating data for variable rates tests
    - Stretching and compressing branches is functionally equivalent to varying rates of evolution because the tree will be detransformed after data simulation, prior to prediction.
    - Our simulation used the following list of scalars...
        - 0.5, 1, 1, 1, 1, 2, and 3
    - This can be set to one or more values greater than 0
    - The simulation is not yet set up to sample from a continuous distribution

- depend_scale1
- depend_adj1
    - These values represent the effect on DEP1 or "low dependency" matrices
    - depend_scale1 multiplies transition rates q13 and q43
        - In our simulation, this was set to 0.5
    - depend_adj1 multiplies the rest of the transition rates, ideally to balance the total evolutionary change
        - In our simulation, this was set to 1.16667
    - These can be set to any numeric value
        - If depend_adj is set to 0, the simulation will not produce trait variation and will fail
        - We recommend that depend_adj is set to (1 + (1 - depend_scale)/3) to balance out the amount of evolutionary change between ER and DEP trials

- depend_scale2
- depend_adj2
    - These values represent the effect on DEP2 or "high dependency" matrices
    - These function the same as the settings above, but they impact different data simulations
    - In our simulation...
        - depend_scale2 = 0.1
        - depend_adj2 = 1.3
    - These can be set to any numeric value
        - If depend_adj is set to 0, the simulation will not produce trait variation and will fail
        - We recommend that depend_adj is set to (1 + (1 - depend_scale)/3) to balance out the amount of evolutionary change between ER and DEP trials
        - These do not have to be more "extreme" than the previous settings (scale2 < scale1 & adj2 > adj1), but it is recommended to help the user interpret results


- multiple_prediction
    - This setting determines if we are predicting for multiple, randomly selected tips
    - This setting increases run time significantly when true
    - This setting must be set to...
        - true
        - false

- clade_prediction
    - This setting determines if we are predicting for a random clade of tips
    - This setting increases run time significantly when true
    - This setting must be set to...
        - true
        - false

- unknown_size
    - This setting determines the amount of unknown tips to be predicted when either of the two previous settings are set to "true"
    - In our simulation, this was set to 10% of the total tips for each run
    - For clade prediction, the actual sampled clade may be slightly smaller than this number, if no clade is found of this size
    - This can be set to any integer between 2 and pop_size

- multistate_prediction
    - This setting determines if we add an additional BayesTraits prediction using the Multistate model
    - The two models currently used are Pagel's Independent and Dependent model, which are both 2-trait models
    - This model instead only looks at 1 trait, ignoring the other simulated trait
    - This setting must be set to...
        - true
        - false

### File Organization

To start the simulation, your main directory should follow this organizational setup before running.

- Discrete_Simulation.*.sh
- BayesTraitsV5.exe
- Scripts/
    - This folder will contain all of the R scripts, which are called by the shell script from this folder

The first R script, "SetupDirectories.*.R", will set up the necessary filing system for the rest of the simulation. It will follow this pattern in the main directory.

- Trees/
    - This folder will contain a number of simulated trees equal to the value of num_iterations
    - These will be used to simulated data for each combination of data simulation settings

- Results/
    - All results will be found sorted into one of the two to three following folders
    - Random/
        - This results folder will contain a text file for the prediction results for each prediction type that was run (Single, Multiple, Clade)
        - It will also contain a folder containing the extended results of the Ancestral State Reconstruction ("AncestralStateReconstruction/")
    - ConstantRates/
        - This results folder will contain a subfolder for the prediction results for each prediction type that was run (Single, Multiple, Clade)
        - As single-tip prediction is always run, we will use this to describe all three folders
        - Single/
            - This folder contains a text file for the prediction results for each matrix used to simulate data (e.g. "ER.L" or "DEP1.M")
            - It will also contain a folder containing the extended results of the Ancestral State Reconstruction for each matrix used to simulate data ("AncestralStateReconstruction/")
    - (VariableRates/)
        - If the variable_rates setting is "true", then this folder will be created to hold all of the files relevant for those trials
        - The structure will match the "ConstantRates" subfolder described above

- ConstantRates/
    - This folder will contain all input and output files for the trials which used data simulated with constant rates
    - They first get sorted into subfolders based on the transition matrix used to simulated data, we will use "ER.L" to describe the structure of all of the other subfolders
    - ER.L/
        - This folder will contain all the relevant files for trials using data simulated with Constant Rates and an ER.L transition matrix
        - They will get sorted into the following subfiles
        - Data/
            - This folder will contain the raw simulated data before sampling and an additional file describing the state distributions for that iteration
        - Single/
            - This folder contains the input and output files for single-tip prediction
            - The file created first will detail all sampled taxa to be predicted and the trait information for these simulation settings (e.g. Constant rates, ER.L, single-tip prediction)
            - The rest of the files will be input and output files created to run BayesTraits
        - (Multiple/)
            - When multiple_prediction is set to "true" in the settings, this folder will be created to match the "Single" subfolder above, but for multiple, randomly selected tips instead of single-tip prediction
        - (Clade/)
            - When clade_prediction is set to "true" in the settings, this folder will be created to match the "Single" subfolder above, but for a randomly selected clade instead of single-tip prediction

- Random/
    - This folder will contain all input and output files for the trials which used data that had been randomly simulated.
    - It will match the structure and contents of the "ER.L" file described above within the "ConstantRates" folder.

- (VariableRates/)
    - If the variable_rates setting is "true", then this folder will be created to hold all of the files relevant for those trials
    - It will match the structure and contents of the "ConstantRates" folder described above

### Expected Outputs

When you return to the main directory at the end of the simulation, you will find several new folders described above. There will also be a Results folder where you can find the full results table for each test. Below is a quick guide to interpret the column names.

- General Info
    - "Trial_#" - This is the trial number for that Matrix
    - "Taxon_#" - This is the tip sampled from that trial's tree
    - "Trait_A" - This is the unknown taxon's character state for the FIRST trait
    - "Trait_B" - This is the unknown taxon's character state for the SECOND trait
    - "Terminal_Branch_Length" - The length of the unknown taxon's terminal branch
    - "True_4States" - This is the unknown taxon's character state when combining A and B into 1 trait with 4 states
    - "nXX" - This is the number of taxa in that trial's data that has the character states of XX (e.g. 00, 01, 10, 11)
    - "Avg_Sis_A" - This is the mean value of the sister taxon or sister clade's FIRST trait
    - "Avg_Sis_B" - This is the mean value of the sister taxon or sister clade's SECOND trait

- Prediction models (only a part of the column name)
    - "Beta_bin" - These columns are results from the Beta Binomial prediction method
    - "Naive" - These columns are results from the Naive Bayes Classifier method
    - "Multi" - These columns are results from predictions using the Multistate model in BayesTraits (if multistate_prediction is set to "true")
    - "Ind" - These columns are results from predictions using the Independent model in BayesTraits
    - "Dep" - These columns are results from predictions using the Dependent model in BayesTraits
    - "MCMC" - These columns are results from predictions which don't use the Reversible Jump setting in BayesTraits (when applicable)
    - "RJ" - These columns are results from predictions which use the Reversible Jump setting in BayesTraits (when applicable)

- Prediction statistic (only a part of the column name)
    - "_Prob" - The model's predictive probability of a 1 for the unknown taxon's Trait B
    - "_acc" - The model's accuracy in predicting the correct state for the unknown taxon's Trait B
    - "_LL" - The log loss score of that prediction

- Ancestral State Reconstruction columns
    - "Sister_*_MaxLh" - The maximum, standardized likelihood from ancestral state reconstruction using the * model from BayesTraits for the sister node of the unknown taxon
    - "Sister_*_Prediction" - The state combination (1, 2, 3, or 4) which had the maximum likelihood for the sister node of the unknown taxon from the same reconstruction as above

    - "Ancestor_*_MaxLh" - The maximum, standardized likelihood from ancestral state reconstruction using the * model from BayesTraits for the grandparent node of the unknown taxon
    - "Ancestor_*_Prediction" - The state combination (1, 2, 3, or 4) which had the maximum likelihood for the grandparent node of the unknown taxon from the same reconstruction as above

### Most common error

If you notice that the simulation has failed, this was most likely due to a lack of variation in the simulated data of a particular iteration. BayesTraits needs to see a variation in the trait for which it is predicting (Trait B); otherwise, it will fail to run. Here are steps to fix that issue.

First, you can check the Results folder, and look under the nXX columns. These will show which trials are missing the required variation. This is most common in simulations with small trees, especially when paired with low rates of evolution. Other options include checking the most recently recreated log files from BayesTraits, which will most likely be the trial that failed, and finding where the Beta Binomial method predicted a perfect 0 or a 1, which also occurs with no trait variation.

Once the problematic trial(s) has been found, you can...
- Edit the 'predict_data' and 'edited_data' files to add the necessary variation to allow BayesTraits can run. If you do this, rerun the Beta Binomial and Naive Bayes predictions with this new data.
- Rerun all the data simulation with larger trees and/or faster rates, particularly if this is a common issue.
- Resimulate the data for that specific trial.

From here, you can retry running 'CompileResults' to summarize all results.
