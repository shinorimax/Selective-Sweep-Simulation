library(devtools)
# install_github("JuliaPalacios/phylodyn")
library(phylodyn)
library(future.apply)

# --- Helper: list files by regex pattern ---
list_tree_files <- function(data_dir, pattern) {
  files <- list.files(data_dir, pattern = pattern, full.names = TRUE)
  sort(files)
}

# --- Helper: read a single .tree file as multiPhylo ---
read_tree_file <- function(path) {
  tr <- ape::read.tree(path)
  if (inherits(tr, "phylo")) {
    tr <- structure(list(tr), class = "multiPhylo")
  }
  tr
}

# --- Helper: compute Fs for a multiPhylo ---
compute_Fs_for_treefile <- function(tr, tol = 8L, show_progress = TRUE, progressor = NULL) {
  Fs <- vector("list", length(tr))
  if (show_progress && is.null(progressor)) {
    pb <- utils::txtProgressBar(min = 0, max = length(tr), style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }
  for (j in seq_along(tr)) {
    Fs[[j]] <- gen_Fmat(tr[[j]], tol = tol)
    if (!is.null(progressor)) {
      progressor()
    } else if (show_progress) {
      utils::setTxtProgressBar(pb, j)
    }
  }
  Fs
}

# --- Helper: read breaks and compute weights ---
read_break_weights <- function(breaks_dir, stem) {
  csv_path <- file.path(breaks_dir, paste0(stem, "_breaks.csv"))
  if (!file.exists(csv_path)) {
    warning("Missing break file: ", csv_path)
    return(NULL)
  }
  df <- readr::read_csv(csv_path, show_col_types = FALSE)
  df <- df[order(df$tree_index), ]
  df$width <- df$right - df$left
  df$width
}

# --- Helper: weighted average F ---
weighted_average_F <- function(F_list, weights) {
  stopifnot(length(F_list) == length(weights))
  dimF <- dim(F_list[[1]])
  F_sum <- matrix(0, nrow = dimF[1], ncol = dimF[2])
  for (k in seq_along(F_list)) {
    F_sum <- F_sum + weights[k] * F_list[[k]]
  }
  F_sum / sum(weights)
}

# --- Worker: process a single file (safe to export to parallel workers) ---
process_one_file <- function(f, idx, total, fs_dir, data_dir, tol,
                             show_tree_progress, progressor = NULL, show_file_progress = TRUE) {
  stem <- sub("\\.tree$", "", basename(f))
  if (show_file_progress) {
    cat(sprintf("[%d/%d] %s\n", idx, total, basename(f)))
  }

  # Read trees and compute Fs
  tr <- read_tree_file(f)
  Fs <- compute_Fs_for_treefile(tr, tol = tol, show_progress = show_tree_progress, progressor = progressor)

  # Save per-file Fs
  saveRDS(Fs, file = file.path(fs_dir, paste0(stem, ".rds")))

  # Weighted-average F for this file
  weights <- read_break_weights(data_dir, stem)
  if (!is.null(weights) && length(weights) == length(Fs)) {
    wavg <- weighted_average_F(Fs, weights)
  } else {
    wavg <- NULL
  }

  # Free memory inside worker
  rm(tr, Fs)
  gc()

  list(name = basename(f), weighted = wavg)
}

# --- Main driver: process one pattern ---
# --- Main driver: process one pattern ---
process_pattern_to_Fs <- function(pattern, data_dir, out_dir, scenario_name, tol = 8L,
                                  use_parallel = TRUE, n_workers = NULL,
                                  show_tree_progress = TRUE,
                                  show_file_progress = TRUE) {
  
  files <- list_tree_files(data_dir, pattern)
  if (length(files) == 0L) {
    warning("No files found for pattern: ", pattern)
    return(list(files = 0L, weighted = NULL))
  }
  
  # Output folders
  fs_dir <- file.path(out_dir, paste0("Fs_nested_", scenario_name))
  if (!dir.exists(fs_dir)) dir.create(fs_dir, recursive = TRUE)
  
  # Store weighted averages (one per file)
  weighted_list <- vector("list", length(files))
  names(weighted_list) <- basename(files)
  
  if (use_parallel) {
    # Tree-level progress in parallel adds overhead; keep minimal progress updates.
    show_tree_progress <- FALSE
    
    if (requireNamespace("progressr", quietly = TRUE)) {
      progressr::handlers(progressr::handler_txtprogressbar(style = 3))
      total_steps <- length(files)
      p <- progressr::progressor(steps = total_steps)
    } else {
      p <- NULL
    }
    
    res <- progressr::with_progress({
      future.apply::future_lapply(seq_along(files), function(i) {
        out <- process_one_file(files[[i]], i, length(files), fs_dir, data_dir, tol,
                                show_tree_progress, progressor = NULL, show_file_progress = FALSE)
        if (!is.null(p) && show_file_progress) p()
        out
      })
    })
    
  } else {
    res <- lapply(seq_along(files), function(i) {
      process_one_file(files[[i]], i, length(files), fs_dir, data_dir, tol,
                       show_tree_progress, progressor = NULL, show_file_progress = show_file_progress)
    })
  }
  
  # Reconstruct weighted_list in original order
  for (i in seq_along(res)) {
    weighted_list[[i]] <- res[[i]]$weighted
  }
  
  # Save weighted averages list for the scenario
  weighted_path <- file.path(out_dir, paste0("F_weighted_average_", scenario_name, ".rds"))
  saveRDS(weighted_list, weighted_path)
  cat("Saved weighted averages ->", weighted_path, "\n")
  
  list(files = length(files), weighted_path = weighted_path, fs_dir = fs_dir)
}
