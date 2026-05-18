library(tidyverse)
library(phytools)
library(phangorn)
library(picante)

########## Load & process tree & data
# Load tree
TTree = read.tree("./taxa.nwk")

# Number of tips
print(Ntip(TTree))

# Load data
data = read.csv("./data.txt", sep="\t", header=FALSE, col.names = c("binomial","t1","t2"))

# Finds all species with data that aren't in the tree (setdiff finds all values in a that do not occur in b)
in_data_not_tree = data.frame(binomial = setdiff(data$binomial, TTree$tip.label))
length(in_data_not_tree$binomial)


########## Tree Cleaning: Convert soft polytomies to hard (trifurcation in newick format)
# likely artifacts, soft polytomies [0 branch length with massively inflates parameters]
# Test to see if the tree is binary (FALSE if polytomies are present)
cat("Tree is binary: ", is.binary(TTree))

# Identify internal edges
internal_edges <- which(TTree$edge[,2] > length(TTree$tip.label))  # Internal nodes have index > number of tips

# Find internal branches with zero length
zero_length_internal <- internal_edges[TTree$edge.length[internal_edges] == 0]
print(length(zero_length_internal))

# Deletes all branches smaller than tol and collapses the corresponding dichotomies into a multichotomy. 
TTree_HardPT = di2multi(TTree, tol = 1e-08)

# Retest to see if the tree is binary (FALSE if polytomies are present)
cat("Tree is binary: ", is.binary(TTree_HardPT))


########## Tree Cleaning: Remove zero (or close to zero) terminal branches [massively inflates variance]
# Identify terminal branches (tips)
terminal_edges <- which(TTree_HardPT$edge[,2] <= length(TTree_HardPT$tip.label))  # Terminal nodes have index ≤ number of tips

# Find & print terminal branches with very small lengths
small_terminal_branches <- terminal_edges[TTree_HardPT$edge.length[terminal_edges] < 1e-2]
print(length(small_terminal_branches))

# Get tip labels of small branches
small_tip_labels <- TTree_HardPT$tip.label[TTree_HardPT$edge[small_terminal_branches, 2]]

# Identify clades these taxa belong to
affected_nodes <- unique(TTree_HardPT$edge[small_terminal_branches, 1])

# Function to retain only one species per affected clade
prune_clades <- function(tree, affected_nodes) {
  tips_to_remove <- c()
  
  for (node in affected_nodes) {
    # Get all descendant tips of the node
    desc_tips <- extract.clade(tree, node)$tip.label
    
    # If multiple tips exist, keep one and mark others for removal
    if (length(desc_tips) > 1) {
      tips_to_remove <- c(tips_to_remove, desc_tips[-1])  # Remove all but one
    }
  }
  
  # Drop marked tips from the tree
  tree_pruned <- drop.tip(tree, tips_to_remove)
  
  return(tree_pruned)
}

# Apply pruning function
TTree_cleaned <- prune_clades(TTree_HardPT, affected_nodes)

# Report results
cat("Pruned", length(TTree_HardPT$tip.label) - length(TTree_cleaned$tip.label), "taxa from small-branch clades.\n")



########## Taxon matching between tree and data and export
# Sort & move binomial names to row names for PCM libraries
data_cleaned3=data
data_cleaned3 = data_cleaned3 %>% arrange(binomial)
data_cleaned3 = column_to_rownames(data_cleaned3, var = names(data_cleaned3)[1])

# Use Picante to automatically match up data and tree
data_tree_synced = match.phylo.data(TTree_cleaned, data_cleaned3)
cat("Number of taxa that match in data & tree:", nrow(data_tree_synced$data))

# Bring binomial back to a column & export data
df_export <- cbind(binomial = rownames(data_tree_synced$data), data_tree_synced$data)
write_tsv(df_export, "DataCleaned.txt", col_names=FALSE)



########## Move Hieraaetus_morphnoides to correct position
# https://academic.oup.com/biolinnean/article/144/2/blae028/7633854?login=false
# Taxon to move, then drop
taxon_to_move <- "Hieraaetus_morphnoides"
tree_pruned <- drop.tip(data_tree_synced$phy, taxon_to_move)

# Aquila clade tips & MRCA
aquila_tips <- grep("Aquila", tree_pruned$tip.label, value = TRUE)
aquila_mrca_pruned <- getMRCA(tree_pruned, aquila_tips)

# Now find the edge in the PRUNED tree leading to that MRCA
edge_row_pruned <- which(tree_pruned$edge[,2] == aquila_mrca_pruned)
edge_length_pruned <- tree_pruned$edge.length[edge_row_pruned]

# Now bind Hieraaetus_morphnoides mid-edge
TTree_new <- bind.tip(tree_pruned,
                     tip.label = taxon_to_move,
                     where = aquila_mrca_pruned,
                     position = edge_length_pruned / 2)


########## Add Hieraaetus_moorei
# https://www.sciencedirect.com/science/article/abs/pii/S1055790318306328?via%3Dihub
# Locate the node number of Hieraaetus_morphnoides
morphnoides_node <- which(TTree_new$tip.label == "Hieraaetus_morphnoides")

# Use bind.tip to graft moorei onto morphnoides with 2.2 branch length
TTree_final <- bind.tip(TTree_new,
                       tip.label = "Hieraaetus_moorei",
                       where = morphnoides_node,
                       position = 2.2)



########## Write the Tree in format readable by BayesTraits & remove taxa block
# Remove node labels and get the tree in nexus format
TTree_final$node.label <- NULL
nexus_lines=capture.output(write.nexus(TTree_final, file = ""))

# Remove lines between BEGIN TAXA and END;
start_taxa <- grep("BEGIN TAXA;", nexus_lines, ignore.case = TRUE)
end_taxa   <- grep("END;", nexus_lines, ignore.case = TRUE)
end_taxa   <- end_taxa[end_taxa > start_taxa][1]  # find correct END; after TAXA

# Remove the TAXA block
if (length(start_taxa) > 0 && length(end_taxa) > 0) {
  nexus_lines <- nexus_lines[-(start_taxa:end_taxa)]
}

# Save the modified file
writeLines(nexus_lines, "BirdTreeCleaned.tre")

