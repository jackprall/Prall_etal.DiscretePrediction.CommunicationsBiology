calc_log_loss <- function(true, pred) {
  if (pred == 1) {
    log_loss <- -1 * ((true * log(0.999999)) + ((1 - true) * log(1 - 0.999999)))
  } else if (pred == 0) {
    log_loss <- -1 * ((true * log(0.000001)) + ((1 - true) * log(1 - 0.000001)))
  } else {
    log_loss <- -1 * ((true * log(pred)) + ((1 - true) * log(1 - pred)))
  }
  return(log_loss)
}

post_pred_simmap <- function(tree, trait, test_taxon, post_Q) {
  ppd <- vector(length = length(post_Q))
  for (i in 1:length(ppd)) {
    temp <- phytools::make.simmap(
      tree = tree,
      x = trait,
      Q = post_Q[[i]],
      pi = "equal",
      nsim = 1,
      message = FALSE
    )
    final_state <- names(
      tail(
        x = temp$maps[[
          which(temp$edge[, 2] == which(temp$tip.label == test_taxon))
        ]],
        n = 1
      )
    )
    ppd[i] <- as.numeric(substr(x = final_state, start = 2, stop = 2))
  }
  return(ppd)
}

post_pred_simmap_v2 <- function(tree, trait, test_taxon, post_Q, nsim_perQ) {
  ppd <- vector(length = nsim_perQ)
  ppd <- rep(list(ppd), length(post_Q))
  for (i in 1:length(ppd)) {
    for (j in 1:nsim_perQ) {
      temp <- phytools::make.simmap(
        tree = tree,
        x = trait,
        Q = post_Q[[i]],
        pi = "equal",
        nsim = 1,
        message = FALSE
      )
      final_state <- names(
        tail(
          x = temp$maps[[
            which(temp$edge[, 2] == which(temp$tip.label == test_taxon))
          ]],
          n = 1
        )
      )
      ppd[[i]][j] <- as.numeric(substr(x = final_state, start = 2, stop = 2))
    }
  }
  ppd <- unlist(ppd)
  return(ppd)
}

getAllDescendantTips <- function(tree, nodes) {  # from Jack Prall
  tip_labels <- c()
  for (node in nodes) {
    if (node %in% tree$tip.label) {
      # Already a tip
      tip_labels <- c(tip_labels, node)
    } else {
      # Not a tip — get descendants
      desc_nodes <- getDescendants(tree, node)
      # Separate tips and internal nodes
      tips <- desc_nodes[desc_nodes <= Ntip(tree)]
      internals <- desc_nodes[desc_nodes > Ntip(tree)]
      # Add tip labels
      tip_labels <- c(tip_labels, tree$tip.label[tips])
      # Recurse into internal nodes (if any)
      if (length(internals) > 0) {
        tip_labels <- c(tip_labels, getAllDescendantTips(tree, internals))
      }
    }
  }
  return(unique(tip_labels))
}

find_grandparent_node <- function(
  tree_full,
  tree_train,
  test_tip
) {
  sisters_mrca <- phytools::getSisters(
    tree = tree_full,
    node = test_tip,
    mode = "label"
  )
  if (names(sisters_mrca) == "tips") {
    index <- which(tree_train$tip.label == sisters_mrca[[1]])
    tip_nodepath <- ape::nodepath(phy = tree_train)[[index]]
    grandparent_node <- tail(tip_nodepath, 2)[1]
  } else {
    sisters <- getAllDescendantTips(
      tree = tree_full,
      node = sisters_mrca[[1]]
    )
    sisters_mrca <- ape::getMRCA(phy = tree_train, tip = sisters)
    if (sisters_mrca == (ape::Ntip(phy = tree_train) + 1)) {
      grandparent_node <- NA
      print("The test taxon does not have a grandparent node.")
    } else {
      grandparent_node <- phytools::getParent(
        tree = tree_train,
        node = sisters_mrca
      )
    }
  }
  grandparent_node <- as.numeric(grandparent_node) - ape::Ntip(phy = tree_train)
  return(grandparent_node)
}

find_sister_node <- function(
  tree_full,
  tree_train,
  test_tip
) {
  sisters_mrca <- phytools::getSisters(
    tree = tree_full,
    node = test_tip,
    mode = "label"
  )
  if (names(sisters_mrca) == "tips") {
    sister_node <- which(tree_train$tip.label == sisters_mrca[[1]])
  } else {
    sisters <- getAllDescendantTips(
      tree = tree_full,
      node = sisters_mrca[[1]]
    )
    sister_node <- ape::getMRCA(phy = tree_train, tip = sisters)
  }
  return(sister_node)
}

get_scaled_condlh <- function(
  tree_train,
  tip_states,
  n_states,
  q_median,
  node
) {
  temp <- castor::asr_mk_model(
    tree = tree_train,
    tip_states = tip_states,
    Nstates = n_states,
    transition_matrix = q_median,
    include_ancestral_likelihoods = TRUE,
    reroot = FALSE,
    root_prior = "flat"
  )
  scaled_condlh <- temp$ancestral_likelihoods[node, ]
  return(scaled_condlh)
}
