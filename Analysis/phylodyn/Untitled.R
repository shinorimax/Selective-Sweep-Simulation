confusion_matrices <- list()  # Store confusion matrices for each scenario

# Iterate over each selection coefficient scenario
for (scenario_idx in seq_along(scenarios)) {
  cat(sprintf("Processing Scenario %d...\n", scenario_idx))
  
  scenario_list <- scenarios[[scenario_idx]]
  total_trees <- length(scenario_list) * trees_per_simulation  # Total trees in this scenario
  progress_step <- ceiling(total_trees / 100)  # Print progress every 1%
  current_count <- 0  # Counter for processed trees
  
  # Store trees for this scenario, but merge Tree Types 1&2 into Class 1 and 3&4 into Class 2
  tree_classes <- vector("list", 2)  # We now have 2 classes instead of 4
  
  # Collect trees per merged type (Class 1: {1,2}, Class 2: {3,4})
  for (sim_idx in seq_along(scenario_list)) {
    for (tree_idx in 1:trees_per_simulation) {
      class_label <- ifelse(tree_idx %in% c(1,2), 1, 2)  # Merge into two classes
      tree_classes[[class_label]] <- append(tree_classes[[class_label]], list(scenario_list[[sim_idx]][[tree_idx]]))
    }
  }
  
  trees_per_class <- length(tree_classes[[1]])  # Number of trees per merged class
  progress_step_centroid <- ceiling(trees_per_class^2 / 100)  # Print every 0.1%
  current_count_centroid <- 0  # Counter for progress tracking
  
  # Compute centroid for each merged class (mean distance matrix)
  centroids <- list()
  for (c in 1:2) {  # Now we have only 2 classes instead of 4
    centroid_matrix <- Reduce(`+`, lapply(tree_classes[[c]], function(tree) {
      sapply(tree_classes[[c]], function(other_tree) {
        dist_value <- dist_pairwise(tree, other_tree, dist.method = dist_method, weighted = weighted)
        current_count_centroid <<- current_count_centroid + 1
        if (current_count_centroid %% progress_step_centroid == 0) {
          cat(sprintf("Class %d Centroid Progress: %.1f%% completed\n", 
                      c, (current_count_centroid / trees_per_class^2) * 100))
          flush.console()  # Ensures immediate printing
        }
        return(dist_value)  # Return value for sapply()
      })
    })) / length(tree_classes[[c]])
    centroids[[c]] <- centroid_matrix
  }
}
