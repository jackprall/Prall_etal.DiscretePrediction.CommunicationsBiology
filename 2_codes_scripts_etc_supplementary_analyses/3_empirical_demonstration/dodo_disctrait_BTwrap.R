# Prall et al. 2025
# Below script 1) predicts nesting behavior in dodos and 2) tests accuracy of behavior predictions in extant birds
# eubirds nesting behavior data from Storchova and Horak 2017
# Timetree.org topology, february 2025
# IMPORTANT: user must edit line 59 in 'bayestraits_btw_edited.R' to local directory with BayesTraitsV4

# for Emi: 
# setwd("~/Library/CloudStorage/OneDrive-MontanaStateUniversity/LabDrive/Papers/Prall_Discrete_Prediction/Code/TempestFiles/10/RealData")

library(tidyverse)
library(btw)
library(phytools)
library(geiger)
library(rotl)
library(traitdata)
library(janitor)
library(stringr)
source('parse_schedule_btw_edited.R')  # edited functions in package btw
source('bayestraits_btw_edited.R') # so that they allow directory creation, logfile naming, and running BayesTraitsV4
source('allometryFunctions_Apr25.R') # add updated taxonomy to data 

# preparing bird nesting data
original_home <- getwd()
original_dodo <- file.path(path=paste0(original_home, "/Dodo_Local"))
original_accuracy <- file.path(path=paste0(original_home, "/AccuracyTest_Nests"))

Sto17<- eubirds

Sto17.edit<- Sto17 %>% rename(
  binomial= scientificNameStd,
)

# use package rotl to update taxonomy
Sto17.edit.tax<- add_taxonomy(Sto17.edit)

# remove duplicates
# a few rows are identically duplicated now after synonymizing taxa 
Sto17.edit.tax<- Sto17.edit.tax[-which(duplicated(Sto17.edit.tax$binomial_rotl)),]

eubirds_nest<- Sto17.edit.tax %>% select(
  binomial_rotl,
  Nest.type,
  Nest.building
) %>% mutate(
  binomial_rotl = gsub(" ", "_", binomial_rotl),
)

# binning ground vs non-ground nests
eubirds_nest$Nest.type.bin[eubirds_nest$Nest.type %in% c("G", "GC", "G,OA", "G,H")]<- 0 # ground nests
eubirds_nest$Nest.type.bin[eubirds_nest$Nest.type %in% c("OA", "H", "OA,H", "CA")]<- 1 # off-ground nests

# binning sex-specific nest building
eubirds_nest$Nest.building.bin[eubirds_nest$Nest.building %in% c("M", "B")]<- 1 # males participate
eubirds_nest$Nest.building.bin[eubirds_nest$Nest.building %in% c("F")]<- 0 # females only
eubirds_nest$Nest.building.bin[eubirds_nest$Nest.building %in% c("N")]<- NA # neither builds
eubirds_nest<- na.omit(eubirds_nest)

# adding Dodo data row
eubirds_nest_dodo<-rbind(eubirds_nest, c("Raphus_cucullatus","G", "N", "0", "-"))

# add rownames for treedata function
rownames(eubirds_nest)<- eubirds_nest$binomial_rotl
rownames(eubirds_nest_dodo)<- eubirds_nest_dodo$binomial_rotl

bird.tree<- read.newick("birds_species.nwk")
# tip name matching takes a while for the tree
# but needs to happen BEFORE tree is reduced 
bird.names<- tnrs_match_names(bird.tree$tip.label)
bird.tree$tip.label<- gsub(" ", "_", bird.names$unique_name)
remove_tips<- c(which(duplicated(bird.tree$tip.label)), which(is.na(bird.tree$tip.label)))
bird.tree.trim<- drop.tip(bird.tree, remove_tips)


# Manually fixing some of the tip labels
# 1) adding full species names (monospecific genera) to taxa with only a genus name in timetree
# 2) Reducing subspecies in timetree to one species-level tip
# 3) Change conflicting spelling in Phylloscopus (rotl missed this one)
bird.tree.trim$tip.label[bird.tree.trim$tip.label == "Botaurus"]<- "Botaurus_stellaris"
bird.tree.trim$tip.label[bird.tree.trim$tip.label == "Ardeola"]<- "Ardeola_ralloides"
bird.tree.trim$tip.label[bird.tree.trim$tip.label == "Calandrella (genus in Deuterostomia)"]<- "Calandrella_brachydactyla"
bird.tree.trim$tip.label[bird.tree.trim$tip.label == "Coracias"]<- "Coracias_garrulus"
bird.tree.trim$tip.label[bird.tree.trim$tip.label == "Cyanistes_cyanus_tianschanicus"]<- "Cyanistes_cyanus"
bird.tree.trim$tip.label[bird.tree.trim$tip.label == "Elanus"]<- "Elanus_caeruleus"
bird.tree.trim$tip.label[bird.tree.trim$tip.label == "Eremophila_alpestris_atlas"]<- "Eremophila_alpestris"
bird.tree.trim$tip.label[bird.tree.trim$tip.label == "Gyps"]<- "Gyps_fulvus"
bird.tree.trim$tip.label[bird.tree.trim$tip.label == "Hydrobates_pelagicus"]<- "Hydrobates_castro"
bird.tree.trim$tip.label[bird.tree.trim$tip.label == "Ixobrychus"]<- "Ixobrychus_minutus"
bird.tree.trim$tip.label[bird.tree.trim$tip.label == "Lophophanes"]<- "Lophophanes_cristatus"
bird.tree.trim$tip.label[bird.tree.trim$tip.label == "Melanocorypha"]<- "Melanocorypha_calandra"
bird.tree.trim$tip.label[bird.tree.trim$tip.label == "Montifringilla"]<- "Montifringilla_nivalis"
bird.tree.trim$tip.label[bird.tree.trim$tip.label == "Oenanthe_hispanica_hispanica"]<- "Oenanthe_hispanica"
bird.tree.trim$tip.label[bird.tree.trim$tip.label == "Pernis"]<- "Pernis_apivorus"
bird.tree.trim$tip.label[bird.tree.trim$tip.label == "Phalacrocoracidae"]<- "Phalacrocorax_pygmaeus"
bird.tree.trim$tip.label[bird.tree.trim$tip.label == "Phoenicopterus"]<- "Phoenicopterus_roseus"
bird.tree.trim$tip.label[bird.tree.trim$tip.label == "Phylloscopus_sibilatrix"]<- "Phylloscopus_sibillatrix"
bird.tree.trim$tip.label[bird.tree.trim$tip.label == "Puffinus"]<- "Puffinus_mauretanicus"
bird.tree.trim$tip.label[bird.tree.trim$tip.label == "Recurvirostra"]<- "Recurvirostra_avosetta"
bird.tree.trim$tip.label[bird.tree.trim$tip.label == "Regulus_regulus_azoricus"]<- "Regulus_regulus"
bird.tree.trim$tip.label[bird.tree.trim$tip.label == "Remiz"]<- "Remiz_pendulinus"
bird.tree.trim$tip.label[bird.tree.trim$tip.label == "Scolopax"]<- "Scolopax_rusticola"
bird.tree.trim$tip.label[bird.tree.trim$tip.label == "Tachybaptus"]<- "Tachybaptus_ruficollis"
bird.tree.trim$tip.label[bird.tree.trim$tip.label == "Thalasseus_sandvicensis_acuflavidus"]<- "Thalasseus_sandvicensis"
bird.tree.trim$tip.label[bird.tree.trim$tip.label == "Xema"]<- "Xema_sabini"


saveRDS(list(bird.tree.trim, eubirds_nest, eubirds_nest_dodo), file = "birddatprepped.Rdata")
# rds<- readRDS("birddatprepped.Rdata") # if session is interrupted, skip waiting on rotl


# resulting object should have 419 taxa 
eubird.treedat<- treedata(data=na.omit(eubirds_nest), phy=bird.tree.trim)

# double check tree for duplicates
which(duplicated(eubird.treedat$phy$tip.label))

tree<- eubird.treedat$phy
dat_test<- data.frame(eubird.treedat$data[,c(1,4,5)])
dat_test[,2]<- as.numeric(dat_test[,2])
dat_test[,3]<- as.numeric(dat_test[,3])
sums<- dat_test[,2]+dat_test[,3]
table(sums) # observe proportion of 00 and 11 states in data 


# preparing dodo data 
eubird.treedat.dodo<- treedata(data=na.omit(eubirds_nest_dodo), phy=bird.tree.trim)

tree_dodo<- eubird.treedat.dodo$phy
dat_test_dodo<- data.frame(eubird.treedat.dodo$data[,c(1,4,5)])

dat_predict_dodo<- dat_test_dodo
dat_test_dodo$Nest.building.bin[dat_test_dodo$binomial_rotl %in% c("Raphus_cucullatus")]<- "-" # skip the '-' taxon
dat_predict_dodo$Nest.building.bin[dat_predict_dodo$binomial_rotl %in% c("Raphus_cucullatus")]<- "?" # predict the '?' taxon


#########
#########
######### running dodo predictions

# if run is interrupted, return to directory 'original' before restarting
original<- original_dodo
setwd(original)

# bayestraits settings
iter<- "iterations 11000000"
burnin<- "burnin 1000000"
commands_list<- list(
  #c("2", "2", iter, "sample 1000", burnin, "PriorAll exp 100"), # ind
  #c("3", "2", iter, "sample 1000", burnin, "PriorAll exp 100"), #dep
  c("2", "2", iter, "sample 1000", burnin, "RJHP exp 0 5"), #indRJ
  c("3", "2", iter, "sample 1000", burnin, "RJHP exp 0 5") #depRJ
)

run_names<- names(commands_list)<- paste0("eubird_nest_dodo", "_", c(#"ind", "dep", # omitted non-RJ runs for inferior performance
  "indRJ", 
  "depRJ"))
results_list<- list()

for (i in 1:length(commands_list)) {
  temp_command_vec<- commands_list[[i]]
  temp_list<- list()
  temp_command_vec_rate<-c(temp_command_vec, 
                           paste0("SaveModels ", run_names[i], ".bin"), 
                           paste0("logfile ", run_names[i], "rate"))
  temp_list[[1]]<- bayestraits(data = dat_test_dodo, tree = tree_dodo, commands= temp_command_vec_rate, 
                               remove_files = FALSE, 
                               dirname=(paste0(run_names[i], "rate")))
  
  loadPath<- normalizePath(list.files(path= paste0("./", run_names[i], "rate"), pattern= "bin", full.names = T, recursive = T))
  temp_command_vec_predict<- c(temp_command_vec[!grepl('RJHP*', temp_command_vec)], 
                               paste0("LoadModels ", run_names[i], ".bin"), 
                               paste0("logfile ", run_names[i], "precict")) 
  prior_tf<- grepl('RJHP*',temp_command_vec_predict)
  
  if(! TRUE %in% prior_tf) temp_command_vec_predict<- c(temp_command_vec_predict, "PriorAll exp 5")
  
  temp_list[[2]]<- bayestraits(data = dat_predict_dodo, tree = tree_dodo, commands= temp_command_vec_predict, 
                               remove_files = FALSE, 
                               dirname=(paste0(run_names[i], "predict")),
                               loadPath = loadPath)
  results_list[[i]]<- temp_list
}



#########
#########
######### looping over random tips to test accuracy of predictions
# results format is set of folders with species names
# each species folder contains subdirectories for 1) indRJ rates; 2) indRJ predict; 3) depRJ rates; and 4) depRJ predict

original<- original_accuracy
setwd(original)

iter<- "iterations 1000000"
burnin<- "burnin 100000"
commands_list<- list(
  #c("2", "2", iter, "sample 1000", burnin, "PriorAll exp 100"), # ind
  #c("3", "2", iter, "sample 1000", burnin, "PriorAll exp 100"), #dep
  c("2", "2", iter, "sample 1000", burnin, "RJHP exp 0 5"), #indRJ
  c("3", "2", iter, "sample 1000", burnin, "RJHP exp 0 5") #depRJ
)

run_names<- names(commands_list)<- paste0("eubird_nest", "_", c(#"ind", "dep", 
  "indRJ", 
  "depRJ"))

# randomly sample from topology 
nsamples<- 5
tip_list<- sample(419, size=nsamples)
samples_out<- list()

# R will not show that it is active while running loop; check logfile updates to make sure everything is running 

for (j in 1:nsamples) {

  temp_dat<- dat_test
  sample_tip<- tip_list[j]
  sample_tip_name<- temp_dat$binomial_rotl[sample_tip]
  results_list<- list()
  dir.create(file.path(original, sample_tip_name))
  setwd(file.path(original, sample_tip_name))
  
for (i in 1:length(commands_list)) {
  temp_command_vec<- commands_list[[i]]
  temp_list<- list()
  temp_dat_rate<- temp_dat
  temp_dat_rate[,3][sample_tip]<- "-"
  temp_dat_predict<- temp_dat
  temp_dat_predict[,3][sample_tip]<- "?"
  known_val<- temp_dat[,3][sample_tip]

  temp_command_vec_rate<-c(temp_command_vec, 
                           paste0("SaveModels ", sample_tip_name, run_names[i], ".bin"), 
                           paste0("logfile ", sample_tip_name, run_names[i], "rate"))
  temp_list[[1]]<- bayestraits(data = temp_dat_rate, tree = tree, commands= temp_command_vec_rate, 
                               remove_files = FALSE, 
                               dirname=(paste0(sample_tip_name, run_names[i], "rate")))
  
  loadPath<- normalizePath(list.files(path= paste0("./", sample_tip_name, run_names[i], "rate"), pattern= "bin", full.names = T, recursive = T))
  temp_command_vec_predict<- c(temp_command_vec[!grepl('RJHP*', temp_command_vec)], 
                               paste0("LoadModels ",sample_tip_name,  run_names[i], ".bin"), 
                               paste0("logfile ",sample_tip_name,  run_names[i], "predict")) 
  prior_tf<- grepl('RJHP*',temp_command_vec_predict)
  
  if(! TRUE %in% prior_tf) temp_command_vec_predict<- c(temp_command_vec_predict, "PriorAll exp 5")
  
  temp_list[[2]]<- bayestraits(data = temp_dat_predict, tree = tree, commands= temp_command_vec_predict, 
                               remove_files = FALSE, 
                               dirname=(paste0(sample_tip_name, run_names[i], "predict", "_", known_val)),
                               loadPath = loadPath)
  results_list[[i]]<- temp_list
}
  
  samples_out[[j]]<- results_list
  setwd(original)
  
}
