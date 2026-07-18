#!/usr/bin/env Rscript
# ============================================================================
# Hyperparameter sensitivity grid (spec §4)
# ----------------------------------------------------------------------------
# Question: do the headline conclusions -- the confounding gap being large and
# positive, and the per-target residualization verdicts -- survive reasonable
# perturbations of the RF hyperparameters (num.trees, min.node.size, k)?
#
# For each grid cell we record:
#   * observed_gap from confoundGapTest()   -> is the gap stable & positive?
#   * per-target verdict from residualize() -> do any verdicts FLIP across cells?
#
# COMPUTE NOTE: the full grid at production settings is a multi-day job on the
# reference box with no remote compute available. Ships with a REDUCED grid that
# still spans the relevant ranges; widen via env vars to reproduce offline:
#   DURUM_SENS_NTREES (comma list), DURUM_SENS_NODESIZE, DURUM_SENS_K,
#   DURUM_SENS_NPERM
# ============================================================================

suppressMessages({
  if (requireNamespace("durumGEAr", quietly = TRUE)) library(durumGEAr) else
    devtools::load_all(".")
})
data(durumUnits)

getl <- function(v, d) { x <- Sys.getenv(v)
  if (nzchar(x)) as.integer(strsplit(x, ",")[[1]]) else d }
geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }

# ---- Reduced grid (override with env vars) -----------------------------------
grid_ntrees   <- getl("DURUM_SENS_NTREES",   c(80L, 160L, 300L))  # prod: 300-600
grid_nodesize <- getl("DURUM_SENS_NODESIZE", c(3L, 5L, 10L))
grid_k        <- getl("DURUM_SENS_K",        c(5L, 10L))
N_PERM        <- geti("DURUM_SENS_NPERM",    20L)                 # prod: 100
targets <- paste0("gPC", 1:5)

grid <- expand.grid(num.trees = grid_ntrees, min.node.size = grid_nodesize,
                    k = grid_k, KEEP.OUT.ATTRS = FALSE,
                    stringsAsFactors = FALSE)
cat(sprintf("Sensitivity grid: %d cells (n_perm=%d)\n", nrow(grid), N_PERM))

rows <- vector("list", nrow(grid))
t0 <- Sys.time()
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  rz <- residualize(durumUnits, seeds = 42:43, k = g$k,
                    num.trees = g$num.trees, min.node.size = g$min.node.size,
                    n_perm = N_PERM)
  gt <- confoundGapTest(durumUnits, seed = 42, n_perm = N_PERM,
                        num.trees = g$num.trees, min.node.size = g$min.node.size)
  vr <- setNames(as.character(rz$verdict), targets)
  rows[[i]] <- data.frame(
    num.trees = g$num.trees, min.node.size = g$min.node.size, k = g$k,
    observed_gap = gt$observed_gap, gap_p = gt$p_value,
    resid_meanR2 = rz$per_target_R2["gPC2"],   # the strongest within-country axis
    t(setNames(vr, paste0("verdict_", targets))),
    stringsAsFactors = FALSE)
  writeLines(sprintf("cell %d/%d done (%.0fs)", i, nrow(grid),
                     as.numeric(difftime(Sys.time(), t0, units = "secs"))),
             "inst/.sens_progress.txt")
}
res <- do.call(rbind, rows)

# Verdict-flip summary: for each target, how many DISTINCT verdicts appear?
flip <- sapply(targets, function(tg) length(unique(res[[paste0("verdict_", tg)]])))

out <- list(
  settings = list(grid_ntrees = grid_ntrees, grid_nodesize = grid_nodesize,
                  grid_k = grid_k, N_PERM = N_PERM, reduced = TRUE),
  grid = res,
  gap_range = range(res$observed_gap),
  gap_all_positive = all(res$observed_gap > 0),
  verdict_distinct_per_target = flip,
  n_cells = nrow(res))

saveRDS(out, "inst/cache_sensitivity_grid.rds")
cat("\n=== gap range across grid:", paste(round(out$gap_range, 4), collapse=" .. "),
    "| all positive:", out$gap_all_positive, "===\n")
cat("distinct verdicts per target (1 = stable):\n"); print(flip)
cat("SAVED inst/cache_sensitivity_grid.rds\n")
