#' Permutation Test for the Interpolation-vs-Extrapolation Confounding Gap
#'
#' The headline diagnostic of this package is the gap between
#' \code{\link{spatialBlockCV}} (interpolation, grouped by site) and
#' \code{\link{locoCV}} (extrapolation, grouped by country): a large gap is
#' the signature of a model that has memorized country identity rather than
#' learned a transferable climate relationship. This function tests that gap
#' against a null distribution built by \strong{shuffling the country labels}
#' (holding the spatial-block interpolation model and folds fixed) and
#' re-running \code{\link{locoCV}} on the shuffled labels, so the null answers:
#' "how big a gap would we see from a model with no real country-identity
#' memorization, purely from finite-sample noise in re-assigning countries?"
#'
#' \strong{Design note:} only \code{locoCV} needs to be rerun per permutation
#' (interpolation skill does not depend on the country labels), which keeps the
#' permutation loop to one call to \code{.fitPredictMO()} (internal helper) per held-out
#' pseudo-country per permutation - still the dominant cost, so \code{n_perm}
#' defaults modestly and a lower \code{num.trees} is recommended for the null
#' loop relative to the headline run (this does NOT reintroduce the
#' hyperparameter-mismatch problem of issue #5, because here both the observed
#' gap and the null gap use \code{\link{locoCV}} with identical settings -
#' only \code{\link{spatialBlockCV}}, which is permutation-invariant, is
#' computed once).
#'
#' @param data A confident modelling frame (e.g. \code{durumUnits}).
#' @param predictors,targets Predictor and target columns.
#' @param site_group Column to group interpolation folds by. Default \code{"SiteCode"}.
#' @param country_group Column identifying the coarse group shuffled under the
#'   null and used for LOCO. Default \code{"Country"}.
#' @param seed Seed for the (unpermuted) spatial-block and LOCO fits. Default 42.
#' @param n_perm Number of country-label permutations for the null. Default 100.
#' @param num.trees,min.node.size Passed to the internal ExtraTrees fitter (used
#'   for both the observed and null runs, so they are directly comparable).
#'
#' @return An object of class \code{"geaGapTest"} with \code{interp_R2},
#'   \code{loco_R2}, \code{observed_gap}, \code{observed_gap_CI} (bootstrap 95\%
#'   CI, resolution capped by \code{n_perm}), \code{null_gap} (vector of length
#'   \code{n_perm}), \code{null_mean}, \code{null_sd}, \code{p_value}
#'   (one-sided: probability the null gap is as large as observed), and
#'   \code{p_value_CI} (bootstrap 95\% CI of that p-value).
#'
#' @section Expected runtime (this is the slowest function in the package):
#' Each permutation reruns the full \code{\link{locoCV}} loop (one ExtraTrees
#' fit per held-out country), so cost scales as
#' \code{n_perm * (number of countries)} tree fits. Benchmarked on the
#' 1,060-row \code{durumUnits} data (single core): \code{locoCV()} alone takes
#' about 130s at \code{num.trees=100} and about 230s at \code{num.trees=200}.
#' At the \strong{defaults} (\code{num.trees=600}, \code{n_perm=100}), a full
#' call is expected to take on the order of \strong{several hours}. For
#' interactive exploration, drop to \code{num.trees=100, n_perm=20} (a few
#' minutes); the package vignette caches a run at these lighter settings
#' rather than the full defaults. Reserve the full defaults for a final,
#' one-off, cached/offline confirmatory run.
#'
#' @examples
#' \donttest{
#' data(durumUnits)
#' gt <- confoundGapTest(durumUnits, seed = 42, n_perm = 20, num.trees = 100)
#' gt$observed_gap
#' gt$p_value
#' }
#' @seealso \code{\link{spatialBlockCV}}, \code{\link{locoCV}}
#' @export
confoundGapTest <- function(data,
                            predictors = c(paste0("BIO", 1:19), "Altitude",
                                           "Latitude", "Longitude"),
                            targets = paste0("gPC", 1:5),
                            site_group = "SiteCode", country_group = "Country",
                            seed = 42, n_perm = 100,
                            num.trees = 600, min.node.size = 3) {
  predictors <- intersect(predictors, names(data))
  Y <- as.matrix(data[, targets, drop = FALSE])
  n <- nrow(data)

  # Observed interpolation skill (permutation-invariant to country labels)
  sb <- spatialBlockCV(data, predictors = predictors, targets = targets,
                       group = site_group, k = 5, seeds = seed,
                       num.trees = num.trees, min.node.size = min.node.size)
  interp_R2 <- sb$metrics$mean_R2

  # Observed extrapolation skill under the TRUE country labels
  lo_obs <- locoCV(data, predictors = predictors, targets = targets,
                   group = country_group, seed = seed,
                   num.trees = num.trees, min.node.size = min.node.size)
  loco_R2 <- lo_obs$metrics$mean_R2
  observed_gap <- interp_R2 - loco_R2

  # Null: shuffle country labels, rerun locoCV under the SAME hyperparameters
  true_country <- as.character(data[[country_group]])
  null_gap <- numeric(n_perm)
  null_loco <- numeric(n_perm)
  data_perm <- data
  for (b in seq_len(n_perm)) {
    set.seed(2000 + b)
    data_perm[[country_group]] <- sample(true_country)
    lo_null <- locoCV(data_perm, predictors = predictors, targets = targets,
                      group = country_group, seed = seed,
                      num.trees = num.trees, min.node.size = min.node.size)
    null_loco[b] <- lo_null$metrics$mean_R2
    null_gap[b] <- interp_R2 - null_loco[b]
  }
  # one-sided: how often does shuffling produce a gap at least as large as observed?
  p_value <- (sum(null_gap >= observed_gap) + 1) / (n_perm + 1)

  # ---- Confidence intervals (issue #14) ----
  # Bootstrap-resample the null gap distribution (>=200 resamples). The gap CI is
  # the bootstrap CI of the null mean shifted to be centred on the observed gap
  # (i.e. observed_gap +/- the null's sampling spread); the p-value CI is the
  # bootstrap distribution of the exceedance p. Both are capped in resolution by
  # n_perm (the number of null_gap draws).
  n_boot <- 200L
  set.seed(7100)
  boot_nullmean <- vapply(seq_len(n_boot), function(b)
    mean(sample(null_gap, length(null_gap), replace = TRUE)), numeric(1))
  margin <- stats::quantile(boot_nullmean - mean(null_gap), c(0.025, 0.975),
                            names = FALSE)
  observed_gap_CI <- c(observed_gap + margin[1], observed_gap + margin[2])
  names(observed_gap_CI) <- c("2.5%", "97.5%")
  boot_p <- vapply(seq_len(n_boot), function(b) {
    rs <- sample(null_gap, length(null_gap), replace = TRUE)
    (sum(rs >= observed_gap) + 1) / (n_perm + 1)
  }, numeric(1))
  p_value_CI <- stats::quantile(boot_p, c(0.025, 0.975), names = FALSE)
  names(p_value_CI) <- c("2.5%", "97.5%")

  structure(list(
    interp_R2 = interp_R2, loco_R2 = loco_R2, observed_gap = observed_gap,
    observed_gap_CI = observed_gap_CI, p_value_CI = p_value_CI,
    null_gap = null_gap, null_loco = null_loco,
    null_mean = mean(null_gap), null_sd = stats::sd(null_gap),
    p_value = p_value, n_perm = n_perm,
    type = "country-label permutation test of the interpolation-vs-LOCO gap"),
    class = "geaGapTest")
}

#' @export
#' @method print geaGapTest
print.geaGapTest <- function(x, ...) {
  cat("<geaGapTest> ", x$type, "\n", sep = "")
  cat(sprintf("  Interpolation R2 (spatial-block) : %+.4f\n", x$interp_R2))
  cat(sprintf("  Extrapolation R2 (LOCO, true labels): %+.4f\n", x$loco_R2))
  cat(sprintf("  Observed gap                      : %+.4f\n", x$observed_gap))
  if (!is.null(x$observed_gap_CI))
    cat(sprintf("    95%% CI (bootstrap)              : [%+.4f, %+.4f]\n",
                x$observed_gap_CI[1], x$observed_gap_CI[2]))
  cat(sprintf("  Null gap (country-shuffled, n=%d)  : mean=%+.4f  SD=%.4f\n",
              x$n_perm, x$null_mean, x$null_sd))
  cat(sprintf("  p-value (observed gap vs null)     : %.4f\n", x$p_value))
  if (!is.null(x$p_value_CI))
    cat(sprintf("    95%% CI (bootstrap)              : [%.4f, %.4f]\n",
                x$p_value_CI[1], x$p_value_CI[2]))
  # Issue #9: warn if p-value is at the resolution floor
  floor_p <- 1 / (x$n_perm + 1)
  if (abs(x$p_value - floor_p) < 1e-10) {
    cat(sprintf("  WARNING: p-value equals resolution floor (1/%d=%.4f).\n", 
                x$n_perm + 1, floor_p))
    cat("           Rerun with higher n_perm before treating as strongly significant.\n")
  }
  invisible(x)
}
