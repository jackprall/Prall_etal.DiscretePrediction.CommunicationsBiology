# Discrete Trait Prediction

This repository contains R scripts and other supplementary files, such as Excel result files, for "Predicting Discrete Traits in Evolving Systems".

## Contents

- 1_codes_scripts_main_analyses

	- README.md
		- This file is a manual describing the analysis workflow involving the scripts below.

    - Discrete_Simulation.*.sh
        - This shell script calls the others to perform the whole simulation.
        - This script also runs the phylogenetic analysis through the software BayesTraits.
        - This holds the primary settings for this simulation.

    - DiscreteFunctions.*.R
        - This R script details all sub-functions used in the next following scripts.

    - SetupDirectories.*.R
        - This script sets up the directories necessary for the simulation to run.

    - TreeGeneration.*.R
        - This script generates the phylogenetic trees used in the simulation.

    - FullDataGeneration.*.R
        - This script generates the data for every trial in the simulation based on settings in the shell script.

    - SampleSingleTaxa.*.R
        - This script samples one taxon from each data and records its data to be predicted later.

    - SampleMultipleTaxa.*.R
        - This script only runs if multiple_prediction or clade_prediction are "true" in the shell script.
        - It functions the same as SampleSingleTaxa.*.R except it samples a number of taxa equal to the unknown_size setting in the shell script.

    - ResultsMatrixGeneration.*.R
        - This script generates the final results table and fills in the taxon and tree information.

    - FilesForBayesTraits.*.R
        - This script generates the input files necessary to run BayesTraits based on the simulated data and settings.

    - BBandNBPrediction.*.R
        - This script runs the Beta Binomial and Naive Bayes predictions.
        - Results are stored in the main results tables in the Results folder.

    - RunBayesTraits.*.R
        - This script performs the BayesTraits runs for this simulation.
        - The type of BayesTraits runs are determined in the settings in the shell script.

    - CompileBayesTraits.*.R
        - This script compiles the results in the log files from BayesTraits.
        - It stores this information in the main results tables in the Results folder.

    - AncestralStateReconstruction.*.R
        - This state performs Ancestral State Reconstruction on the sister and grandparent nodes of the unknown taxon using models created by BayesTraits.
        - The full results are stored in the "AncestralStateReconstruction" subfolders within the Results folder.
        - Summary stats are pulled from this table into the main results table.

    - SummarizeResults.*.R
        - This script summarizes the results into more digestable summary tables.
        - These are stored in the Results folder.

    - RemoveSuperfluousFiles.*.R
        - This script removes any additional files created in the middle of the simulation to decrease the storage costs of the simulation long term.
        - It removes instruction, schedule, model, and particular data files which are only used as inputs for BayesTraits.

- 2_codes_scripts_etc_supplementary_analyses

	- This folder contains files related to analyses in the Supplementary Information.

- 3_codes_scripts_etc_figures

	- This folder contains files related to generating all figures in the manuscript and Supplementary Information.

- 4_results_excel

	- This folder contains Excel result files: those that detail the results for each simulation trial and those that summarize comparisons of interest.
