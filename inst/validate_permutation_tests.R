#!/usr/bin/env Rscript
# ============================================================================
# Permutation-test calibration study (spec §3)
# ----------------------------------------------------------------------------
# Question: are the permutation-based p-values in residualize() and
# confoundGapTest() well-calibrated -- i.e. under a TRUE null (no real
# within-country climate signal), do they reject at ~= alpha?
#
# Design (Type I error simulation):
#   * residualize(): destroy any real within-country signal by shuffling each
#     target WITHIN country before the fit, then run the test. A calibrated test
#     rejects (p < 0.05) about 5% of the time across replicates.
#   * confoundGapTest(): the interpolation-vs-LOCO gap under a true null is
#     estimated by shuffling COUNTRY labels before computing the observed gap, so
#     there is no genuine country structure for the test to detect.
#
# COMPUTE NOTE: full production settings (n_perm = 100, num.trees >= 300, >= 200
# replicates) run for many hours on the reference box and there is no remote
# compute configured here. This script therefore ships with REDUCED defaults
# (see the block below) that finish in minutes and still demonstrate the
# calibration behaviour and the method. Override via the env vars to reproduce
# at production settings offline:
#   DURUM_CAL_NREP, DURUM_CAL_NPERM, DURUM_CAL_NTREES, DURUM_CAL_K
# ============================================================================

suppressMessages({
  if (requireNamespace("durumGEAr", quietly = TRUE)) library(durumGEAr) else
    devtools::load_all(".")
})
data(durumUnits)

geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }

# ---- Reduced-but-documented defaults (override with env vars) ----------------
# residualize Type I loop (fast: k-fold, single seed):
N_REP   <- geti("DURUM_CAL_NREP",   30L)   # production: >= 200
N_PERM  <- geti("DURUM_CAL_NPERM",  20L)   # production: 100
N_TREES <- geti("DURUM_CAL_NTREES", 100L)  # production: 300-600
K       <- geti("DURUM_CAL_K",       5L)
# confoundGapTest Type I loop tuned SEPARATELY -- each rep reruns LOCO over ~43
# countries n_perm times, so it is ~2 orders of magnitude heavier per rep. Its
# reduced defaults are correspondingly smaller; widen offline on bigger hardware.
N_REP_GAP  <- geti("DURUM_CAL_NREP_GAP",  8L)    # production: >= 200
N_PERM_GAP <- geti("DURUM_CAL_NPERM_GAP", 10L)   # production: 100
N_TREES_GAP<- geti("DURUM_CAL_NTREES_GAP",60L)   # production: 300-600
ALPHA   <- 0.05
targets <- paste0("gPC", 1:5)

cat(sprintf("Calibration: resid[N_REP=%d N_PERM=%d N_TREES=%d K=%d]  gap[N_REP=%d N_PERM=%d N_TREES=%d]\n",
            N_REP, N_PERM, N_TREES, K, N_REP_GAP, N_PERM_GAP, N_TREES_GAP))

# ---- residualize Type I: shuffle each target within country --------------------
shuffle_within_country <- function(df, seed) {
  set.seed(seed)
  for (tg in targets)
    for (ct in unique(df$Country)) {
      idx <- which(df$Country == ct)
      if (length(idx) > 1) df[[tg]][idx] <- sample(df[[tg]][idx])
    }
  df
}

rej_resid <- matrix(NA, N_REP, length(targets), dimnames = list(NULL, targets))
t0 <- Sys.time()
for (r in seq_len(N_REP)) {
  dfn <- shuffle_within_country(durumUnits, seed = 1000 + r)
  rz <- residualize(dfn, seeds = 42, k = K, num.trees = N_TREES, n_perm = N_PERM)
  rej_resid[r, ] <- as.integer(rz$p_value < ALPHA & rz$obs_R2_test > 0)
  if (r %% 5 == 0)
    writeLines(sprintf("resid rep %d/%d (%.0fs elapsed)", r, N_REP,
                       as.numeric(difftime(Sys.time(), t0, units = "secs"))),
               "inst/.cal_progress.txt")
}
type1_resid <- colMeans(rej_resid, na.rm = TRUE)

# ---- confoundGapTest Type I: shuffle country labels ---------------------------
rej_gap <- numeric(N_REP_GAP)
t0 <- Sys.time()
for (r in seq_len(N_REP_GAP)) {
  dfn <- durumUnits
  set.seed(2000 + r)
  dfn$Country <- sample(dfn$Country)          # break country structure
  gt <- confoundGapTest(dfn, seed = 42, n_perm = N_PERM_GAP,
                        num.trees = N_TREES_GAP, min.node.size = 5)
  rej_gap[r] <- as.integer(gt$p_value < ALPHA)
  writeLines(sprintf("gap rep %d/%d (%.0fs elapsed)", r, N_REP_GAP,
                     as.numeric(difftime(Sys.time(), t0, units = "secs"))),
             "inst/.cal_progress.txt")
  saveRDS(list(rej_gap = rej_gap[seq_len(r)]), "inst/.cal_gap_partial.rds")
}
type1_gap <- mean(rej_gap, na.rm = TRUE)

result <- list(
  settings = list(N_REP = N_REP, N_PERM = N_PERM, N_TREES = N_TREES, K = K,
                  N_REP_GAP = N_REP_GAP, N_PERM_GAP = N_PERM_GAP,
                  N_TREES_GAP = N_TREES_GAP, ALPHA = ALPHA, reduced = TRUE),
  type1_resid = type1_resid,      # per-target false-positive rate
  type1_resid_overall = mean(type1_resid),
  type1_gap = type1_gap,
  n_rep = N_REP, n_rep_gap = N_REP_GAP)

saveRDS(result, "inst/cache_permutation_calibration.rds")
cat("\n=== Type I error (should be ~", ALPHA, ") ===\n")
print(round(type1_resid, 3))
cat("residualize overall:", round(mean(type1_resid), 3),
    "| confoundGapTest:", round(type1_gap, 3), "\n")
cat("SAVED inst/cache_permutation_calibration.rds\n")
