# Internal helper: verdict logic (issue #8), isolated for unit-testing without
# paying for any ranger fit or permutation loop. The verdict must require BOTH a
# low permutation p-value AND a positive observed residual R2 -- a significant
# p-value alone is not enough, because the permutation null can be beaten by a
# target whose true residual R2 is negative (worse than the country-mean
# baseline), which is definitionally NOT "real within-country signal".
.residualVerdict <- function(p_value, obs_R2_test, alpha) {
  ifelse(p_value < alpha & obs_R2_test > 0,
         "real within-country signal", "country-identity artifact")
}

# Internal helper: fragility score for a single per-target residualization
# verdict (issue #15). Formalizes the three shakiness conditions that the print
# methods warned about in v1.1.1 (issues #8/#9) into one graded field. In
# increasing severity: (1) near-zero observed R2 (abs < 0.01) -> "marginal";
# (2) resolution-floor p-value (p == 1/(n_perm+1)) -> "marginal"; (3) sign
# disagreement between the single-seed obs_R2 (which drives the verdict) and the
# multi-seed per_target_R2 headline -> "fragile" (the sign itself is
# seed-dependent). Later conditions override earlier ones. per_target_R2 must be
# the multi-seed value for this same target.
.fragility_score <- function(obs_R2, p_value, per_target_R2, n_perm) {
  score <- 0L  # robust
  if (abs(obs_R2) < 0.01) score <- 1L                              # marginal: tiny effect
  if (abs(p_value - 1 / (n_perm + 1)) < 1e-10) score <- 1L         # marginal: floor p
  if (!isTRUE(all.equal(sign(obs_R2), sign(per_target_R2),
                        tolerance = 0.1))) score <- 2L             # fragile: sign flip
  c("robust", "marginal", "fragile")[score + 1L]
}

#' Fold-Safe Country-Mean Residualization (Confound Diagnosis)
#'
#' Measures how much within-country climate signal each target carries, by
#' removing between-country variance and asking whether climate still predicts
#' what remains. For each fold, each country's mean target is computed from
#' \emph{training rows only}, subtracted from the target, and an ExtraTrees model
#' is fit to the residual. A residual-scale R2 at or below zero means the target
#' is a pure country-identity artifact; a positive residual R2 means genuine
#' within-country climate signal.
#'
#' \strong{Intuition:} to test whether climate explains more than "which country
#' am I in", first subtract each country's average, then predict the leftover.
#' If the model can still explain the leftover, that is real, transferable signal;
#' if not, the earlier skill was just a country lookup.
#'
#' The fold-safe design is essential: computing country means on the full data
#' (including test rows) leaks the answer and inflates the residual R2. Optionally,
#' \code{shrink = TRUE} applies James-Stein / empirical-Bayes shrinkage of each
#' country mean toward the global mean, which is algebraically the random-intercept
#' BLUP of a mixed model and stabilizes small-country means.
#'
#' \strong{Significance test (not just a raw cutoff):} a hardcoded threshold on
#' the residual R2 (e.g. "R2 > 0.02 means real signal") has no distributional
#' justification - a positive residual R2 can arise from noise alone,
#' especially with the sparse per-country sample sizes in this dataset (median
#' 4 rows/country). This function therefore also runs a permutation test: the
#' target is shuffled \strong{within each country} (preserving country-level
#' structure and sample sizes exactly, but destroying any true climate-target
#' relationship), the identical fold-safe residualization + ExtraTrees fit is
#' rerun on the shuffled target, and this is repeated \code{n_perm} times to
#' build a null distribution of residual R2 under "no real within-country
#' signal". The verdict is based on the resulting p-value (default threshold
#' 0.05), with the raw R2 reported alongside for context rather than used alone.
#'
#' @param data A confident modelling frame (e.g. \code{durumUnits}).
#' @param predictors,targets Predictor and target columns.
#' @param country Column identifying the confounding group. Default \code{"Country"}.
#' @param group Column to build spatial-block folds on. Default \code{"SiteCode"}.
#' @param k,seeds Fold count and seeds (averaged). Defaults 5 and \code{42:46}.
#' @param shrink Logical; apply empirical-Bayes shrinkage of country means. Default FALSE.
#' @param num.trees,min.node.size Passed to the internal ExtraTrees fitter.
#' @param n_perm Number of within-country permutations for the significance test,
#'   per target. Default 100. The permutation test uses the fold assignment and
#'   hyperparameters of \code{seeds[1]} only (not averaged across all
#'   \code{seeds}), so the observed statistic it is compared against
#'   (\code{obs_R2_test}) is a single-seed value - reported alongside the
#'   multi-seed-averaged \code{per_target_R2} used for the headline number.
#' @param alpha Significance threshold for the verdict. Default 0.05.
#'
#' @return An object of class \code{"geaResidual"} with \code{metrics} (residual-scale
#'   \code{geaMetrics}), \code{per_target_R2} (multi-seed-averaged),
#'   \code{per_seed_R2} (the full seeds x targets matrix behind that mean, exposed
#'   so callers can compute their own across-seed summaries),
#'   \code{per_target_R2_CI} (2 x targets; 2.5\%/97.5\% empirical quantiles of the
#'   per-seed R2 - coarse at the default 5 seeds), \code{obs_R2_test}
#'   (single-seed statistic used in the permutation test), \code{p_value} (per
#'   target), \code{p_value_CI} (2 x targets; bootstrap CI from resampling the
#'   permutation null, resolution capped by \code{n_perm}), \code{null_mean},
#'   \code{null_sd} (per target), \code{verdict}
#'   (per-target "real within-country signal" / "country-identity artifact",
#'   based on BOTH \code{p_value < alpha} AND \code{obs_R2_test > 0}),
#'   \code{verdict_fragility} (per-target "robust"/"marginal"/"fragile" grade of
#'   how shaky that verdict is) and \code{shrink}.
#'
#' @examples
#' \donttest{
#' data(durumUnits)
#' rz <- residualize(durumUnits, seeds = 42, n_perm = 50)
#' rz$verdict     # which targets keep signal after removing country identity
#' rz$p_value
#' }
#' @section Expected runtime:
#' The headline \code{per_target_R2} averages \code{length(seeds)} full
#' fold-safe fits across all \code{length(targets)} targets; the permutation
#' test then reruns that same per-target fit \code{n_perm} times (single
#' seed). Benchmarked on \code{durumUnits} at the defaults
#' (\code{num.trees=600}, 5 targets, 1 seed for the permutation loop): about
#' 75s per permutation, so \code{n_perm=100} (the default) takes on the order
#' of \strong{2 hours}. For interactive exploration, drop to
#' \code{num.trees=100-150, n_perm=20-30} (several minutes); the package
#' vignette caches a run at lighter settings rather than the full default.
#'
#' @seealso \code{\link{spatialBlockCV}}, \code{\link{driverAnalysis}}
#' @export
residualize <- function(data,
                        predictors = c(paste0("BIO", c(2,5,6,8,9,12,14,15,16,17,19)),
                                       "Altitude"),
                        targets = paste0("gPC", 1:5),
                        country = "Country", group = "SiteCode",
                        k = 5, seeds = 42:46, shrink = FALSE,
                        num.trees = 600, min.node.size = 3,
                        n_perm = 100, alpha = 0.05) {
  predictors <- intersect(predictors, names(data))
  Y <- as.matrix(data[, targets, drop = FALSE])
  ctry <- as.character(data[[country]])
  grp <- data[[group]]

  # Internal: one fold-safe country-mean residual R2 for an arbitrary target
  # vector y_vec (the true targets, or a within-country-shuffled version of
  # them for the permutation null), under a fixed fold seed s.
  .residR2 <- function(y_vec, s) {
    folds <- .groupFolds(grp, k = k, seed = s)
    oof <- numeric(nrow(data)); true_res <- numeric(nrow(data))
    for (f in sort(unique(folds))) {
      te <- which(folds == f); tr <- which(folds != f)
      ytr <- y_vec[tr]
      cm <- tapply(ytr, ctry[tr], mean); gm <- mean(ytr)
      if (shrink) {
        nj <- tapply(ytr, ctry[tr], length)
        within_var <- tapply(ytr, ctry[tr], function(v) if (length(v) > 1) stats::var(v) else 0)
        s2w <- mean(within_var, na.rm = TRUE)
        s2b <- max(0, stats::var(cm) - s2w / mean(nj))
        B <- s2b / (s2b + s2w / nj); B[is.na(B)] <- 0
        cm <- gm + B * (cm - gm)
      }
      cm_tr <- ifelse(ctry[tr] %in% names(cm), cm[ctry[tr]], gm)
      cm_te <- ifelse(ctry[te] %in% names(cm), cm[ctry[te]], gm)
      ytr_res <- ytr - cm_tr; yte_res <- y_vec[te] - cm_te
      df <- data.frame(data[tr, predictors, drop = FALSE], .y = ytr_res)
      rf <- ranger::ranger(dependent.variable.name = ".y", data = df,
                           num.trees = num.trees, min.node.size = min.node.size,
                           splitrule = "extratrees", num.random.splits = 1,
                           mtry = length(predictors), seed = s, num.threads = 1)
      oof[te] <- stats::predict(rf, data = data[te, predictors, drop = FALSE])$predictions
      true_res[te] <- yte_res
    }
    ss_tot <- sum((true_res - mean(true_res))^2)
    1 - sum((true_res - oof)^2) / ss_tot
  }

  per_target_acc <- matrix(0, length(seeds), length(targets),
                           dimnames = list(NULL, targets))
  for (si in seq_along(seeds)) {
    s <- seeds[si]
    for (j in seq_along(targets)) per_target_acc[si, j] <- .residR2(Y[, j], s)
  }
  mean_pt <- colMeans(per_target_acc)
  metrics <- structure(list(
    per_target = data.frame(target = targets, R2 = mean_pt, RMSE = NA_real_,
                            stringsAsFactors = FALSE),
    mean_R2 = mean(mean_pt), mean_RMSE = NA_real_), class = "geaMetrics")

  # ---- Permutation significance test (issue #6): shuffle target WITHIN each
  # country (preserves country structure/sizes exactly), rerun the identical
  # fold-safe residualization + fit, build a null distribution per target, and
  # base the verdict on the resulting p-value rather than a raw cutoff alone.
  s1 <- seeds[1]
  obs_R2_test <- vapply(seq_along(targets), function(j) .residR2(Y[, j], s1), numeric(1))
  names(obs_R2_test) <- targets
  null_mat <- matrix(NA_real_, n_perm, length(targets), dimnames = list(NULL, targets))
  for (b in seq_len(n_perm)) {
    set.seed(3000 + b)
    for (j in seq_along(targets)) {
      yp <- Y[, j]
      for (g in unique(ctry)) { ix <- which(ctry == g); yp[ix] <- sample(yp[ix]) }
      null_mat[b, j] <- .residR2(yp, s1)
    }
  }
  p_value <- vapply(seq_along(targets), function(j)
    (sum(null_mat[, j] >= obs_R2_test[j]) + 1) / (n_perm + 1), numeric(1))
  names(p_value) <- targets
  null_mean <- colMeans(null_mat); null_sd <- apply(null_mat, 2, stats::sd)

  # ---- Confidence intervals (issue #14) ----
  # (a) Residual R2 CI from the across-seed distribution. With the default 5
  #     seeds this is a coarse empirical quantile; widen `seeds` for a finer CI.
  if (length(seeds) >= 2) {
    per_target_R2_CI <- apply(per_target_acc, 2, stats::quantile,
                              probs = c(0.025, 0.975), names = FALSE)
  } else {
    per_target_R2_CI <- matrix(NA_real_, 2L, length(targets))
  }
  dimnames(per_target_R2_CI) <- list(c("2.5%", "97.5%"), targets)

  # (b) p-value CI by bootstrap resampling of the permutation null (>=200
  #     resamples). NOTE: resolution is capped by n_perm (the number of rows in
  #     null_mat); with n_perm=100 the interval cannot be finer than ~0.01.
  n_boot <- 200L
  p_value_CI <- matrix(NA_real_, 2L, length(targets),
                       dimnames = list(c("2.5%", "97.5%"), targets))
  set.seed(7000)
  for (j in seq_along(targets)) {
    nj <- null_mat[, j]
    boot_p <- vapply(seq_len(n_boot), function(b) {
      rs <- sample(nj, length(nj), replace = TRUE)
      (sum(rs >= obs_R2_test[j]) + 1) / (n_perm + 1)
    }, numeric(1))
    p_value_CI[, j] <- stats::quantile(boot_p, c(0.025, 0.975), names = FALSE)
  }

  # ---- Fragility scoring (issue #15): grade each verdict's shakiness, reusing
  # the floor-p / sign-disagreement conditions the v1.1.1 print method already
  # warned about rather than adding a parallel warning path.
  verdict_fragility <- vapply(seq_along(targets), function(j)
    .fragility_score(obs_R2_test[j], p_value[j], mean_pt[j], n_perm),
    character(1))
  names(verdict_fragility) <- targets

  verdict <- .residualVerdict(p_value, obs_R2_test, alpha)
  names(verdict) <- targets
  structure(list(metrics = metrics, per_target_R2 = mean_pt,
                 per_seed_R2 = per_target_acc,
                 per_target_R2_CI = per_target_R2_CI,
                 obs_R2_test = obs_R2_test, p_value = p_value,
                 p_value_CI = p_value_CI,
                 null_mean = null_mean, null_sd = null_sd, null_mat = null_mat,
                 n_perm = n_perm,
                 alpha = alpha, verdict = verdict,
                 verdict_fragility = verdict_fragility,
                 shrink = shrink,
                 type = sprintf("fold-safe country-mean residual%s",
                                if (shrink) " (EB-shrunk)" else "")),
            class = "geaResidual")
}

#' @export
#' @method print geaResidual
print.geaResidual <- function(x, ...) {
  cat("<geaResidual> ", x$type, "\n", sep = "")
  cat(sprintf("  Residual-scale per-target-mean R2 (multi-seed headline): %+.4f\n", x$metrics$mean_R2))
  cat(sprintf("  Verdict: BOTH p < %.2f AND R2 > 0 required (n_perm=%d, within-country shuffle):\n",
              x$alpha, x$n_perm))
  cat("  (R2 below is obs_R2_test, the single-seed statistic the verdict is actually\n")
  cat("   computed from -- NOT per_target_R2 -- because it can differ in sign from the\n")
  cat("   multi-seed-averaged headline number above; see per_target_R2 for that context.)\n")
  floor_p <- 1 / (x$n_perm + 1)
  any_floor <- FALSE
  any_disagree <- FALSE
  has_ci <- !is.null(x$per_target_R2_CI) && !is.null(x$p_value_CI)
  has_frag <- !is.null(x$verdict_fragility)
  for (t in names(x$verdict)) {
    cat(sprintf("    %-6s obs_R2_test=%+.4f  (per_target_R2=%+.4f)  p=%.4f  -> %s\n",
                t, x$obs_R2_test[t], x$per_target_R2[t], x$p_value[t], x$verdict[t]))
    if (has_ci && all(is.finite(x$per_target_R2_CI[, t]))) {
      cat(sprintf("           95%% CI: R2 [%+.4f, %+.4f]  p [%.4f, %.4f]\n",
                  x$per_target_R2_CI[1, t], x$per_target_R2_CI[2, t],
                  x$p_value_CI[1, t], x$p_value_CI[2, t]))
    }
    if (has_frag) {
      frag <- x$verdict_fragility[t]
      cat(sprintf("           FRAGILITY: %s\n", frag))
      if (frag == "fragile")
        cat("           WARNING: verdict is shaky; do not use this target for downstream analysis without more seeds/n_perm.\n")
    }
    if (abs(x$p_value[t] - floor_p) < 1e-10) any_floor <- TRUE
    if (sign(x$obs_R2_test[t]) != sign(x$per_target_R2[t])) any_disagree <- TRUE
  }
  # Issue #9: warn if any p-value is at the resolution floor
  if (any_floor) {
    cat(sprintf("  WARNING: one or more p-values equal the resolution floor (1/%d=%.4f).\n",
                x$n_perm + 1, floor_p))
    cat("           Rerun with higher n_perm before treating those targets as strongly significant.\n")
  }
  if (any_disagree) {
    cat("  WARNING: for at least one target, obs_R2_test (single-seed, drives the verdict)\n")
    cat("           and per_target_R2 (multi-seed-averaged headline) disagree in sign.\n")
    cat("           Treat the verdict as fragile for that target; increase seeds/n_perm to confirm.\n")
  }
  invisible(x)
}
