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
