# Tutorial.

# Load packages. ----
library(export)
library(ggtree)
library(phytools)
library(tidyverse)

# Load the phylogenetic tree. ----
tree <- ape::read.nexus(file = "tutorial_tree.nex")

# Specify the test tip. ----
test_tip <- "t100"

# Load and prepare the trait dataset. ----
traits <- readr::read_tsv(file = "tutorial_traits.tsv", col_names = FALSE)
colnames(traits) <- c("tip", "X", "Y")
traits <- traits |>
    dplyr::select(tip, X, Y) |>
    dplyr::mutate(set = ifelse(tip == test_tip, "test", "training")) |>
    dplyr::mutate(set = as.factor(set))
traits <- as.data.frame(traits)
rownames(traits) <- traits$tip

# Create the plot. ----
set_levels <- levels(traits$set)
set_list <- vector(mode = "list", length = length(set_levels))
for (j in 1:length(set_levels)) {
  set_list[[j]] <- traits$tip[traits$set == set_levels[j]]
}
names(set_list) <- set_levels
tree_object <- ggtree::groupOTU(tree, set_list)
plot_tree_trait <- ggtree::ggtree(
  tr = tree_object,
  mapping = aes(color = .data$group),
  linewidth = 0.25
) +
  theme_tree2() +
  scale_color_manual(values = c("red", "dark gray")) +
  labs(
    caption = "Arbitrary time unit"
  )
plot_tree_trait <- ggtree::revts(treeview = plot_tree_trait)
plot_tree_trait <- ggtree::gheatmap(
  p = plot_tree_trait,
  data = traits[, c("X", "Y")],
  offset = 0,
  width = 0.075,
  low = "gray",
  high = "black",
  colnames = FALSE
) +
  scale_x_ggtree() +
  theme(legend.position = "none", plot.caption = element_text(hjust = 0))
print(plot_tree_trait)
export::graph2pdf(
  file = "tutorial_plot_tree_trait",
  width = 6.5 / 1.618,
  height = 6.5,
  font = "Arial",
  bg = "transparent",
  colormodel = "rgb"
)
