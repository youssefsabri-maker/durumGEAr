## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>", fig.width = 6, fig.height = 4.5
)

## ----cache-load, include = FALSE----------------------------------------------
## This vignette displays REAL COMPUTED NUMBERS, not eval=FALSE stubs
## (adversarial-review issue #1). The cross-validation and permutation-test
## chunks below (spatialBlockCV, naiveCV, locoCV, confoundGapTest,
## residualize, driverAnalysis, scoreThenCluster) are individually cheap to
## call, but the FULL set of them together (in particular confoundGapTest's
## and residualize's permutation loops) takes on the order of 30-40 minutes
## even at reduced settings. Re-running all of that on every `R CMD build`
## is impractical, so the numbers are pre-computed once by
## `data-raw/build_vignette_cache.R` (using deliberately lighter settings
## than the package defaults -- fewer seeds/trees/permutations -- purely to
## keep the vignette build tractable) and cached to
## `inst/extdata/vignette_cache.rds`. Every code chunk below is still the
## exact call that produced the cached object; only the *execution* is
## replaced by a cache lookup, and the object printed is identical to what
## calling the function live would return.
cache_path <- system.file("extdata", "vignette_cache.rds", package = "durumGEAr")
if (!nzchar(cache_path)) {
  # package not yet installed (e.g. building the vignette in place during
  # development) -- fall back to the source-tree copy
  cache_path <- "../inst/extdata/vignette_cache.rds"
}
vc <- readRDS(cache_path)

## ----setup--------------------------------------------------------------------
library(durumGEAr)

## -----------------------------------------------------------------------------
data(durumRaw)
dim(durumRaw)

cu <- collapseUnits(durumRaw)
print(cu)

## -----------------------------------------------------------------------------
data(durumUnits)
dim(durumUnits)
stopifnot(nrow(cu$confident) == nrow(durumUnits))

## -----------------------------------------------------------------------------
ep <- computeEPC(durumUnits)
round(ep$var_explained, 3)

## -----------------------------------------------------------------------------
## sb <- spatialBlockCV(durumUnits, seeds = 42:44, num.trees = 300, min.node.size = 3)
## (cached; see "cache-load" chunk above -- lighter settings than the
##  num.trees = 600 package default, purely to keep the vignette build fast)
sb <- vc$sb
sb

## -----------------------------------------------------------------------------
## nv <- naiveCV(durumUnits, seeds = 42:44, num.trees = 300, min.node.size = 3)
nv <- vc$naive
nv

## -----------------------------------------------------------------------------
## lo <- locoCV(durumUnits, group = "Country", num.trees = 150, min.node.size = 3)
lo <- vc$loco
lo$metrics$mean_R2
lo$median_group_R2
nrow(lo$per_group)              # number of countries with >= 2 members
sum(lo$per_group$beats_baseline)   # how many countries beat a no-skill baseline

## lc <- locoCV(durumUnits, group = "Cluster", num.trees = 150, min.node.size = 3)
lc <- vc$loclu
lc$metrics$mean_R2
lc$per_group

## -----------------------------------------------------------------------------
## gt <- confoundGapTest(durumUnits, seed = 42, n_perm = 60,
##                        num.trees = 15, min.node.size = 3)
## (n_perm and num.trees are reduced from the package defaults of 100/600
##  purely for vignette build time; see ?confoundGapTest "Expected runtime")
gt <- vc$gaptest
gt

## -----------------------------------------------------------------------------
## rz <- residualize(durumUnits, seeds = 42:43, n_perm = 100,
##                    num.trees = 50, min.node.size = 3)
## (num.trees reduced from the default of 600 purely for vignette build
##  time; see ?residualize "Expected runtime")
rz <- vc$rz
rz

## -----------------------------------------------------------------------------
data.frame(target = names(rz$obs_R2_test),
           obs_R2_test = round(rz$obs_R2_test, 3),
           per_target_R2 = round(rz$per_target_R2, 3),
           p_value = round(rz$p_value, 3),
           verdict = rz$verdict,
           row.names = NULL)

## -----------------------------------------------------------------------------
## dr <- driverAnalysis(durumUnits, "gPC2", seeds = 42, n_perm = 40,
##                       num.trees = 150, min.node.size = 5)
## (n_perm and num.trees reduced from defaults of 200/300; see
##  ?driverAnalysis "Expected runtime")
dr <- vc$dr_gpc2
head(dr$univariate)
dr$p_value
dr$n_predictors_tested

## -----------------------------------------------------------------------------
## s2 <- scoreThenCluster(durumUnits, per_tree = FALSE, cv_ceiling = TRUE,
##                         num.trees = 300, min.node.size = 3)
s2 <- vc$s2
s2$accuracy      # honest out-of-fold cluster accuracy
s2$ceiling       # accuracy a regressor could reach, per ceiling_type below
s2$ceiling_type  # which ceiling estimator was used
s2$gap           # ceiling - accuracy = the cost of Stage-1 regression error
s2$confusion

## -----------------------------------------------------------------------------
s2$entropy_calibration
s2$state_accuracy
sum(s2$off_manifold)

## ----echo = FALSE, results = "asis"-------------------------------------------
targets <- names(rz$per_target_R2)
for (tg in targets) {
  cat(sprintf(
    "- **%s**: interpolation R2 = %.3f (spatial-block); within-country residual R2 = %.3f (p = %.3f) -> *%s*.\n",
    tg, vc$sb$metrics$per_target$R2[vc$sb$metrics$per_target$target == tg],
    rz$per_target_R2[tg], rz$p_value[tg], rz$verdict[tg]
  ))
}

## ----fig.alt = "Map of accession locations coloured by genetic cluster"-------
mapAccessions(durumUnits, color_by = "Cluster",
              main = "Durum wheat effective units")

## ----stability, echo = FALSE, results = "asis"--------------------------------
if (!is.null(vc$sensitivity)) {
  se <- vc$sensitivity
  cat(sprintf(
    "Across %d grid cells, the observed confounding gap ranged from %.3f to %.3f (median %.3f). ",
    se$n_cells, se$gap_range[1], se$gap_range[2], stats::median(se$grid$observed_gap)))
  flips <- se$verdict_distinct_per_target
  stable <- names(flips)[flips == 1]
  unstable <- names(flips)[flips > 1]
  cat(sprintf("Verdict stability: %s held a single verdict across every cell",
              if (length(stable)) paste(stable, collapse = ", ") else "no target"))
  if (length(unstable))
    cat(sprintf("; %s flipped across the grid and should be treated as borderline.",
                paste(unstable, collapse = ", "))) else
    cat(" — no target flipped.")
  cat("\n")
} else {
  cat("*The sensitivity-grid cache is not bundled in this build. ",
      "Run `Rscript inst/sensitivity_analysis.R` (env vars `DURUM_SENS_*` set the ",
      "grid) to generate `inst/cache_sensitivity_grid.rds`, then rebuild the ",
      "vignette cache.*\n", sep = "")
}

## ----calibration, echo = FALSE, results = "asis"------------------------------
if (!is.null(vc$calibration)) {
  ca <- vc$calibration
  cat(sprintf(
    "At the reduced study settings (residualize: N_REP=%d, N_PERM=%d, num.trees=%d; gap: N_REP=%d, N_PERM=%d, num.trees=%d), the empirical Type I error was %.3f for `residualize()` and %.3f for `confoundGapTest()`, against a nominal alpha of %.2f. ",
    ca$settings$N_REP, ca$settings$N_PERM, ca$settings$N_TREES,
    ca$settings$N_REP_GAP, ca$settings$N_PERM_GAP, ca$settings$N_TREES_GAP,
    ca$type1_resid_overall, ca$type1_gap, ca$settings$ALPHA))
  cat("A rate near alpha indicates the permutation test is calibrated; a rate ",
      "well above alpha would warn that the null is too easy to beat at these ",
      "settings.\n", sep = "")
} else {
  cat("*The permutation-calibration cache is not bundled in this build. ",
      "Run `Rscript inst/validate_permutation_tests.R` (env vars `DURUM_CAL_*` ",
      "set replicate/permutation counts) to generate ",
      "`inst/cache_permutation_calibration.rds`, then rebuild the vignette cache.*\n",
      sep = "")
}

## ----deploy-fit---------------------------------------------------------------
data(durum_residualize_results)
mod <- fitGeneticScoreModel(
  durumUnits,
  residualize_result = durum_residualize_results,
  num.trees = 100                    # light for the vignette; default is 300
)
mod

## ----deploy-predict-----------------------------------------------------------
newx <- durumUnits[1:5, ]
pr   <- predict(mod, newx)
round(pr, 3)
attr(pr, "gated_targets")           # axes the shipped validation does not certify

## ----deploy-robust------------------------------------------------------------
rs <- robustGeneticScore(mod, newx, keep = c("gPC2", "gPC5"))
data.frame(robust_score = round(rs, 3),
           axes = paste(attr(rs, "axes"), collapse = ", "))

## ----deploy-extrap------------------------------------------------------------
er <- checkExtrapolationRisk(mod, newx)
er[, c("mahalanobis", "md_flag", "n_out_of_range", "extrapolation")]

## ----deploy-cluster, eval = FALSE---------------------------------------------
# s2 <- scoreThenCluster(durumUnits, per_tree = FALSE, fit_final = TRUE,
#                        num.trees = 300, min.node.size = 3)
# pc <- predictCluster(mod, newx, s2$final_classifier)
# pc                                    # predicted cluster, posterior, entropy

