library(phylodyn)

# Read all files matching a regex pattern; return a named list of 'multiPhylo'
# (If a file has a single tree, we wrap it into a multiPhylo of length 1)
read_files_by_pattern <- function(dir, pattern) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  files <- sort(files)
  out <- lapply(files, function(f) {
    tr <- read.tree(f)                         # phylo or multiPhylo
    if (inherits(tr, "phylo")) {
      tr <- structure(list(tr), class = "multiPhylo")
    }
    tr
  })
  names(out) <- basename(files)
  out
}


# --- Function to flatten nested tree lists ----------------------------------
flatten_trees <- function(tree_lists) {
  out <- list()
  idx <- 1L
  for (i in seq_along(tree_lists)) {
    trees_i <- tree_lists[[i]]
    for (j in seq_along(trees_i)) {
      out[[idx]] <- trees_i[[j]]
      attr(out[[idx]], "src_idx") <- c(i = i, j = j)  # keep track of source
      idx <- idx + 1L
    }
  }
  out
}


# --- Compute & save F-matrices while preserving nesting ---------------------
compute_and_save_Fs_nested <- function(obj, name_prefix, out_dir = "data", tol = 8L,
                                       show_progress = TRUE) {
  n_chr <- length(obj)
  total_trees <- sum(vapply(obj, length, integer(1L)))
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  if (show_progress) pb <- txtProgressBar(min = 0, max = total_trees, style = 3)
  
  Fs_nested <- vector("list", n_chr)
  if (!is.null(names(obj))) names(Fs_nested) <- names(obj)
  
  k <- 0L
  for (i in seq_along(obj)) {
    trees_i <- obj[[i]]
    Fi <- vector("list", length(trees_i))
    if (!is.null(names(trees_i))) names(Fi) <- names(trees_i)
    
    for (j in seq_along(trees_i)) {
      Fi[[j]] <- gen_Fmat(trees_i[[j]], tol = tol)
      attr(Fi[[j]], "src_idx") <- c(i = i, j = j)  # optional trace-back
      k <- k + 1L
      if (show_progress) setTxtProgressBar(pb, k)
    }
    Fs_nested[[i]] <- Fi
  }
  if (show_progress) close(pb)
  
  out_path <- file.path(out_dir, paste0("Fs_nested_", name_prefix, ".rds"))
  saveRDS(Fs_nested, out_path)
  cat("Saved nested F-matrices ->", out_path, "\n")
  
  invisible(Fs_nested)
}


compute_weighted_Fs <- function(Fs_nested, breaks_dir,
                                output_name = "F_weighted_neutrals.rds",
                                store_dir = "data") {
  #' Compute weighted-average F-matrices from a nested list of F-matrices
  #'
  #' @param Fs_nested   Nested list: Fs_nested[[i]][[j]] = F-matrix for jth tree of ith chromosome
  #' @param breaks_dir  Directory containing *_breaks.csv files
  #' @param output_name Output .rds filename (stored in `data_dir`)
  #' @param data_dir    Directory to store output file (must exist)
  #'
  #' @return List of weighted-average F-matrices (one per chromosome)
  
  # Assign names to Fs_nested if missing (use CSV stems)
  if (is.null(names(Fs_nested))) {
    csv_files <- list.files(breaks_dir, pattern = "_breaks\\.csv$")
    stems <- sub("_breaks\\.csv$", "", csv_files)
    names(Fs_nested) <- stems
  }
  
  # Helper: weighted average of F-matrices
  weighted_average_F <- function(F_list, weights) {
    stopifnot(length(F_list) == length(weights))
    dimF <- dim(F_list[[1]])
    F_sum <- matrix(0, nrow = dimF[1], ncol = dimF[2])
    for (k in seq_along(F_list)) {
      F_sum <- F_sum + weights[k] * F_list[[k]]
    }
    F_sum / sum(weights)
  }
  
  # Initialize result container
  F_weighted_list <- vector("list", length(Fs_nested))
  names(F_weighted_list) <- names(Fs_nested)
  
  # Iterate over each chromosome/scenario
  for (i in seq_along(Fs_nested)) {
    stem <- names(Fs_nested)[i]
    csv_path <- file.path(breaks_dir, paste0(sub("\\.tree$", "", stem), "_breaks.csv"))
    if (!file.exists(csv_path)) {
      warning("⚠️ Missing break file: ", csv_path)
      next
    }
    
    # Read genomic intervals
    df_breaks <- read_csv(csv_path, show_col_types = FALSE)
    df_breaks$width <- df_breaks$right - df_breaks$left
    
    # Extract F-matrices and align with tree order
    F_list <- Fs_nested[[i]]
    if (length(F_list) != nrow(df_breaks)) {
      warning(sprintf("⚠️ %s: mismatch (%d Fs vs %d CSV rows)",
                      stem, length(F_list), nrow(df_breaks)))
      next
    }
    
    df_breaks <- df_breaks[order(df_breaks$tree_index), ]
    weights <- df_breaks$width
    
    # Compute weighted mean
    F_weighted_list[[i]] <- weighted_average_F(F_list, weights)
    
    cat(sprintf("✅ %s: weighted F-matrix from %d trees computed.\n",
                stem, length(F_list)))
  }
  
  # Save the result as .rds
  output_path <- file.path(store_dir, output_name)
  saveRDS(F_weighted_list, output_path)
  cat(sprintf("💾 Saved weighted F-matrices to %s\n", output_path))
  
  invisible(F_weighted_list)
}

compute_empirical_kingman <- function(F_weighted_neutrals, output_path = "data/F_empirical_kingman.rds") {
  # --- 1. Take the first half of the trees -----------------------------------
  n_total <- length(F_weighted_neutrals)
  n_half  <- floor(n_total / 2)
  subset_Fs <- F_weighted_neutrals[1:n_half]
  
  # --- 2. Compute average F-matrix ------------------------------------------
  # Assuming each element is a numeric matrix of equal dimension
  F_avg <- Reduce("+", subset_Fs) / n_half
  
  # --- 3. Save to RDS -------------------------------------------------------
  saveRDS(F_avg, file = output_path)
  
  # --- 4. Return the averaged matrix ----------------------------------------
  return(F_avg)
}

compute_empirical_kingman_2 <- function(F_weighted_neutrals, output_path = "data/F_empirical_kingman.rds") {
  # --- 1. Take the first half of the trees -----------------------------------
  n_total <- length(F_weighted_neutrals)
  n_half  <- floor(n_total / 2)
  subset_Fs <- F_weighted_neutrals[(n_half + 1):n_total]
  
  # --- 2. Compute average F-matrix ------------------------------------------
  # Assuming each element is a numeric matrix of equal dimension
  F_avg <- Reduce("+", subset_Fs) / n_half
  
  # --- 3. Save to RDS -------------------------------------------------------
  saveRDS(F_avg, file = output_path)
  
  # --- 4. Return the averaged matrix ----------------------------------------
  return(F_avg)
}

break_coalescent_ties <- function(tr, eps = 1e-8) {
  n <- tr$Nnode + 1L
  depths <- ape::node.depth.edgelength(tr)
  ages   <- max(depths) - depths
  int    <- (n+1):(n + tr$Nnode)
  a      <- ages[int]
  
  # enforce strictly increasing internal ages
  o <- order(a)
  a_sorted <- a[o]
  for (i in 2:length(a_sorted)) {
    if (a_sorted[i] <= a_sorted[i-1]) a_sorted[i] <- a_sorted[i-1] + eps
  }
  a[o] <- a_sorted
  ages[int] <- a
  
  # rebuild edge lengths from parent/child ages
  el <- numeric(nrow(tr$edge))
  for (e in seq_len(nrow(tr$edge))) {
    p <- tr$edge[e,1]; c <- tr$edge[e,2]
    el[e] <- ages[p] - ages[c]
    if (el[e] < 0) stop("Negative branch length after jitter; decrease eps.")
  }
  tr$edge.length <- el
  tr
}
# usage:
# tr <- break_coalescent_ties(tr)
# gen.tr.data2(tr)

# --- Compute & save Weighted F-matrices (FWs) while preserving nesting -------
compute_and_save_FWs_nested <- function(obj, name_prefix, out_dir = "data", tol = 8L,
                                        show_progress = TRUE) {
  # obj: nested list of trees (e.g., from read_files_by_pattern)
  # name_prefix: prefix for output file name (e.g., "neutrals")
  # tol: numerical tolerance
  # out_dir: directory for saving results

  n_chr <- length(obj)
  total_trees <- sum(vapply(obj, length, integer(1L)))

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  if (show_progress) pb <- txtProgressBar(min = 0, max = total_trees, style = 3)

  FWs_nested <- vector("list", n_chr)
  if (!is.null(names(obj))) names(FWs_nested) <- names(obj)

  k <- 0L
  for (i in seq_along(obj)) {
    trees_i <- obj[[i]]
    FW_i <- vector("list", length(trees_i))
    if (!is.null(names(trees_i))) names(FW_i) <- names(trees_i)

    for (j in seq_along(trees_i)) {
      # --- Generate F and W matrices for this tree ---
      tr <- break_coalescent_ties(trees_i[[j]])
      tr_info <- phylodyn:::gen.tr.data2(tr, tol = tol)
      Fmat <- tr_info$Fmat0
      Wmat <- tr_info$Wmat0

      # --- Weighted F-matrix ---
      FW_i[[j]] <- Fmat * Wmat

      # optional trace-back attribute
      attr(FW_i[[j]], "src_idx") <- c(i = i, j = j)

      k <- k + 1L
      if (show_progress) setTxtProgressBar(pb, k)
    }
    FWs_nested[[i]] <- FW_i
  }
  if (show_progress) close(pb)

  out_path <- file.path(out_dir, paste0("FWs_nested_", name_prefix, ".rds"))
  saveRDS(FWs_nested, out_path)
  cat("Saved nested weighted F-matrices ->", out_path, "\n")

  invisible(FWs_nested)
}

compute_weighted_Fs_around_positive_mutation <- function(
    Fs_nested,
    breaks_dir,
    base,      # mutation position in bp
    range,     # genomic range in bp (±range)
    output_name = "F_weighted_filtered.rds",
    store_dir = "data"
) {
  if (!dir.exists(store_dir)) dir.create(store_dir, recursive = TRUE)
  
  # Weighted average
  weighted_average_F <- function(F_list, weights) {
    dimF <- dim(F_list[[1]])
    F_sum <- matrix(0, nrow = dimF[1], ncol = dimF[2])
    for (k in seq_along(F_list)) {
      F_sum <- F_sum + weights[k] * F_list[[k]]
    }
    F_sum / sum(weights)
  }
  
  F_weighted_list <- vector("list", length(Fs_nested))
  names(F_weighted_list) <- names(Fs_nested)
  
  for (i in seq_along(Fs_nested)) {
    
    stem <- names(Fs_nested)[i]
    stem2 <- sub("\\.tree$", "", stem)
    csv_path <- file.path(breaks_dir, paste0(stem2, "_breaks.csv"))
    
    if (!file.exists(csv_path)) {
      warning("⚠️ Missing break file: ", csv_path)
      next
    }
    
    df_breaks <- readr::read_csv(csv_path, show_col_types = FALSE)
    df_breaks$width <- df_breaks$right - df_breaks$left
    
    F_full <- Fs_nested[[i]]
    if (length(F_full) != nrow(df_breaks)) {
      stop(sprintf("❌ mismatch: %d F-matrices but %d CSV rows",
                   length(F_full), nrow(df_breaks)))
    }
    
    # ------------------------------
    # ✔ Correct window-based filtering
    # ------------------------------
    win_start <- base - range
    win_end   <- base + range
    
    df_breaks$keep <- (df_breaks$right >= win_start) &
      (df_breaks$left  <= win_end)
    
    idx_keep <- which(df_breaks$keep)
    
    if (length(idx_keep) == 0) {
      warning(sprintf("⚠️ No trees overlap window in %s", stem))
      next
    }
    
    # Filter F matrices and weights
    F_filtered <- F_full[idx_keep]
    weights <- df_breaks$width[idx_keep]
    
    # Compute weighted F matrix
    F_weighted_list[[i]] <- weighted_average_F(F_filtered, weights)
    
    cat(sprintf(
      "✅ %s: %d trees used (pos %d–%d).\n",
      stem, length(idx_keep), win_start, win_end
    ))
  }
  
  output_path <- file.path(store_dir, output_name)
  saveRDS(F_weighted_list, output_path)
  cat(sprintf("💾 Saved weighted F-matrices to %s\n", output_path))
  
  invisible(F_weighted_list)
}

?phylodyn:::gen.tr.data2
