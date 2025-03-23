# Define the function
generate_distance_matrices <- function(scenarios, dist_method = "l1", weighted = FALSE, trees_per_simulation = 4) {
  # List to store distance matrices for each selection coefficient scenario
  distance_matrices <- list()
  
  # Iterate over each scenario in the input list
  for (scenario_idx in seq_along(scenarios)) {
    # Extract the F matrices list for the current selection scenario
    scenarios_list <- scenarios[[scenario_idx]]
    
    # Calculate the number of simulations per scenario
    scenario_length <- length(scenarios_list)
    scenario_num <- scenario_length / trees_per_simulation
    
    # Initialize a list to store distance matrices for each simulation in this scenario
    scenario_dmat_list <- list()
    
    # Iterate over each simulation in the scenario
    for (index in 1:scenario_length) {
      # Initialize the distance matrix
      dmat <- matrix(0, nrow = trees_per_simulation, ncol = trees_per_simulation)
      
      # Calculate pairwise distances
      for (i in 1:trees_per_simulation) {
        for (j in 1:i) {
          dmat[i, j] <- dist_pairwise(
            scenarios_list[[index]][[i]], 
            scenarios_list[[index]][[j]], 
            dist.method = dist_method, 
            weighted = weighted
          )
        }
      }
      
      # Mirror the lower triangle to the upper triangle
      tmp <- dmat[lower.tri(dmat)]
      dmat <- t(dmat)
      dmat[lower.tri(dmat)] <- tmp
      
      # Store the distance matrix for the current simulation
      scenario_dmat_list[[index]] <- dmat
    }
    
    # Store the distance matrices for the current scenario
    distance_matrices[[scenario_idx]] <- scenario_dmat_list
  }
  
  # Return the list of distance matrices
  return(distance_matrices)
}



# Define the function
plot_mds_scenarios <- function(distance_matrices, trees_per_simulation, tree_colors = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) {
  
  plot_list <- list()  # Store plots for each scenario
  
  # Iterate over each scenario in the distance matrices
  for (scenario_idx in seq_along(distance_matrices)) {
    scenario_dmat_list <- distance_matrices[[scenario_idx]]
    
    # Prepare a data frame to hold MDS points for all simulations in the scenario
    plot_data <- data.frame()
    variance_explained <- c()  # To store variance explained for each dimension
    
    for (sim_idx in seq_along(scenario_dmat_list)) {
      dmat <- scenario_dmat_list[[sim_idx]]
      
      # Check if the distance matrix is valid
      if (any(is.na(dmat)) || any(!is.finite(dmat)) || nrow(dmat) < 2 || ncol(dmat) < 2) {
        cat(sprintf("Skipping invalid distance matrix in Scenario %d, Simulation %d\n", scenario_idx, sim_idx))
        next  # Skip this simulation
      }
      
      # Perform MDS on the valid distance matrix
      mds <- cmdscale(dmat, eig=TRUE, k=2)
      
      # Ensure MDS produced valid results
      if (is.null(mds$points) || ncol(mds$points) < 2 || sum(mds$eig) == 0) {
        cat(sprintf("Skipping non-invertible matrix in Scenario %d, Simulation %d\n", scenario_idx, sim_idx))
        next  # Skip this simulation
      }
      
      # Calculate the variance explained by each dimension
      variance_explained <- mds$eig[1:2] / sum(mds$eig) * 100
      
      # Create a data frame for this simulation's MDS points
      sim_data <- data.frame(
        Dimension1 = mds$points[, 1],
        Dimension2 = mds$points[, 2],
        TreeType = factor(paste("Tree", 1:trees_per_simulation)),
        Simulation = sim_idx
      )
      
      # Combine with overall plot data
      plot_data <- rbind(plot_data, sim_data)
    }
    
    # If there are no valid data points, skip plotting
    if (nrow(plot_data) == 0) {
      cat(sprintf("Skipping Scenario %d due to no valid MDS data\n", scenario_idx))
      next
    }
    
    # Add centroids for each tree type
    centroids <- plot_data %>%
      group_by(TreeType) %>%
      summarize(
        Centroid1 = mean(Dimension1, na.rm = TRUE),
        Centroid2 = mean(Dimension2, na.rm = TRUE)
      )
    
    # Determine symmetric limits based on the max absolute value in plot_data
    # x_limits <- quantile(plot_data$Dimension1, probs = c(0.005, 0.995), na.rm = TRUE)
    # y_limits <- quantile(plot_data$Dimension2, probs = c(0.005, 0.995), na.rm = TRUE)
    x_limits <- quantile(plot_data$Dimension1, probs = c(0, 1), na.rm = TRUE)
    y_limits <- quantile(plot_data$Dimension2, probs = c(0, 1), na.rm = TRUE)
    x_range <- range(x_limits)
    y_range <- range(y_limits)
    
    # Create the ggplot for this scenario
    p <- ggplot(plot_data, aes(x=Dimension1, y=Dimension2, color=TreeType)) +
      geom_point(size=1, alpha=0.8) +
      geom_point(data=centroids, aes(x=Centroid1, y=Centroid2, color=TreeType), 
                 size=5, shape=17, stroke=1.2, inherit.aes=FALSE) +
      scale_color_manual(values=tree_colors) +
      coord_fixed(ratio=1, xlim=x_range, ylim=y_range) +
      labs(
        title=paste("MDS Plot - Scenario", scenario_idx),
        color="Tree Type",
        x=paste("Dimension 1 (", round(variance_explained[1], 2), "%)", sep=""),
        y=paste("Dimension 2 (", round(variance_explained[2], 2), "%)", sep="")
      ) +
      theme(legend.position="right")
    
    # Store the plot
    plot_list[[length(plot_list) + 1]] <- p
    
    print(p)
  }
  # Arrange plots in a 2x2 grid
  final_plot <- ggarrange(plotlist = plot_list, ncol = 2, nrow = 2, common.legend = TRUE, legend = "right")
  
  # Display the final grid of plots
  print(final_plot)
}


# Define the function
calculate_pairwise_stats <- function(distance_matrices, indices = c(1, 2, 3, 4), conf_level = 0.95) {
  # Check if indices are valid
  if (any(indices > length(distance_matrices))) {
    stop("Indices exceed the number of elements in the distance_matrices list.")
  }
  
  # Pairs to report
  pairs <- list(
    c(1, 2), c(1, 3), c(1, 4),
    c(2, 3), c(2, 4),
    c(3, 4)
  )
  
  # Iterate over specified indices
  for (index in indices) {
    # Get the list of matrices for the current index
    matrices_list <- distance_matrices[[index]]
    
    # Number of matrices
    num_matrices <- length(matrices_list)
    
    # Stack matrices for CI computation
    stacked_values <- array(0, dim = c(4, 4, num_matrices))
    
    # Sum each [i, j] element across all matrices
    for (i in seq_along(matrices_list)) {
      matrix <- matrices_list[[i]]
      stacked_values[, , i] <- matrix
    }
    
    # Function to calculate mean for bootstrapping
    calculate_mean <- function(data, indices) {
      apply(data[, , indices, drop = FALSE], c(1, 2), mean)
    }
    
    # Bootstrap confidence intervals
    boot_result <- boot(
      data = stacked_values,
      statistic = calculate_mean,
      R = 1000
    )
    
    # Print results for unique pairs
    cat("\nResults for [[", index, "]]:\n", sep = "")
    for (pair in pairs) {
      i <- pair[1]
      j <- pair[2]
      
      # Average and CI for the pair
      avg <- mean(stacked_values[i, j, ])
      ci <- quantile(boot_result$t[, i + (j - 1) * 4], probs = c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2))
      
      # Print the results
      cat(sprintf(
        "Pair [%d, %d]: Average = %.2f, CI = [%.2f, %.2f]\n",
        i, j, avg, ci[1], ci[2]
      ))
    }
  }
}


# Function to compute pairwise tree distances
dist_matrix <- function(trees, dist_method = "l2", weighted = TRUE) {
  n <- length(trees)
  dmat <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in i:n) {
      dmat[i, j] <- dist_pairwise(trees[[i]], trees[[j]], dist.method = dist_method, weighted = weighted)
    }
  }
  dmat <- dmat + t(dmat) - diag(diag(dmat)) # Mirror the lower triangle
  return(dmat)
}


generate_mds_confusion_matrices <- function(distance_matrices, trees_per_simulation = 4) {
  confusion_matrices <- list()  # Store confusion matrices for each scenario
  
  for (scenario_idx in seq_along(distance_matrices)) {
    cat(sprintf("\nProcessing Scenario %d...\n", scenario_idx))
    
    scenario_dmat_list <- distance_matrices[[scenario_idx]]
    num_sims <- length(scenario_dmat_list)
    total_trees <- num_sims * trees_per_simulation
    progress_step_sim <- ceiling(num_sims / 10)  # Print every 10% of MDS computations
    progress_step_class <- ceiling(total_trees / 10)  # Print every 10% of classifications
    
    current_count_sim <- 0  # Counter for MDS calculations
    current_count_class <- 0  # Counter for tree classifications
    
    # Data storage
    all_mds_data <- data.frame()
    centroids <- list()
    
    # Track MDS computation progress
    for (sim_idx in seq_along(scenario_dmat_list)) {
      dmat <- scenario_dmat_list[[sim_idx]]
      
      # Validate distance matrix
      if (any(is.na(dmat)) || any(!is.finite(dmat)) || nrow(dmat) < 2 || ncol(dmat) < 2) {
        # cat(sprintf("Skipping invalid distance matrix in Scenario %d, Simulation %d\n", scenario_idx, sim_idx))
        next
      }
      
      # Perform MDS
      mds <- cmdscale(dmat, eig=TRUE, k=2)
      
      # Ensure MDS results are valid
      if (is.null(mds$points) || ncol(mds$points) < 2 || sum(mds$eig) == 0) {
        # cat(sprintf("Skipping non-invertible matrix in Scenario %d, Simulation %d\n", scenario_idx, sim_idx))
        next
      }
      
      # Store MDS results
      sim_data <- data.frame(
        Dimension1 = mds$points[, 1],
        Dimension2 = mds$points[, 2],
        TreeType = factor(ifelse(rep(1:trees_per_simulation, each = nrow(mds$points) / trees_per_simulation) %in% c(1,2), 1, 2)), # Merge 1&2 -> Class 1, 3&4 -> Class 2
        Simulation = sim_idx
      )
      
      # Combine with all simulation data
      all_mds_data <- rbind(all_mds_data, sim_data)
      
      # Update MDS progress
      current_count_sim <- current_count_sim + 1
      if (current_count_sim %% progress_step_sim == 0) {
        # cat(sprintf("Scenario %d MDS Progress: %.1f%% completed\n", scenario_idx, (current_count_sim / num_sims) * 100))
        flush.console()
      }
    }
    
    # Skip if no valid data
    if (nrow(all_mds_data) == 0) {
      # cat(sprintf("Skipping Scenario %d due to no valid MDS data\n", scenario_idx))
      next
    }
    
    # Compute MDS-based centroids (2 groups: {1,2} and {3,4})
    # cat(sprintf("Computing centroids for Scenario %d...\n", scenario_idx))
    centroid_df <- all_mds_data %>%
      group_by(TreeType) %>%
      summarize(
        Centroid1 = mean(Dimension1, na.rm = TRUE),
        Centroid2 = mean(Dimension2, na.rm = TRUE)
      )
    
    centroids <- split(centroid_df[, c("Centroid1", "Centroid2")], centroid_df$TreeType)
    
    # Initialize confusion matrix (2x2 instead of 4x4)
    confusion_matrix <- matrix(0, nrow = 2, ncol = 2)
    
    # Track classification progress
    # cat(sprintf("Classifying trees for Scenario %d...\n", scenario_idx))
    for (i in 1:nrow(all_mds_data)) {
      tree_point <- all_mds_data[i, c("Dimension1", "Dimension2")]
      true_label <- as.integer(all_mds_data$TreeType[i])
      
      # Compute distances to MDS centroids
      distances_to_centroids <- sapply(1:2, function(c) {
        sum((tree_point - centroids[[c]])^2)  # Euclidean distance in MDS space
      })
      
      # Assign classification to the closest centroid
      predicted_label <- which.min(distances_to_centroids)
      
      # Update confusion matrix
      confusion_matrix[true_label, predicted_label] <- confusion_matrix[true_label, predicted_label] + 1
      
      # Update classification progress
      current_count_class <- current_count_class + 1
      if (current_count_class %% progress_step_class == 0) {
        # cat(sprintf("Scenario %d Classification Progress: %.1f%% completed\n", scenario_idx, (current_count_class / total_trees) * 100))
        flush.console()
      }
    }
    
    # Store confusion matrix for this scenario
    confusion_matrices[[scenario_idx]] <- confusion_matrix
  }
  
  return(confusion_matrices)
}
