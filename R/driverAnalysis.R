#' Univariate + Permutation-Importance Driver Analysis
#'
#' Ranks which climate predictors drive a target's within-country signal, on the
#' fold-safe residual scale (so the confound removed by \code{\link{residualize}}
#' cannot re-enter). Three complementary lines of evidence are computed:
#' univariate residual R2 (each predictor alone), multivariate permutation
#' importance (from the full joint model), and a within-country permutation
#' significance test (shuffle the target within each country to build a null).
#'
#' \strong{Intuition:} a single-predictor model isolates each variable's marginal
#' effect; permutation importance measures how much the joint model degrades when
#' one predictor is scrambled; the within-country shuffle tests whether the whole
#' association could have arisen by chance while preserving country structure.
#'
#' \strong{Multiple-comparison correction (max-statistic null):} the top-ranked
#' predictor is chosen \emph{after} looking at all univariate R2 values, so
#' testing it against a null built from re-fitting \emph{only that one
#' predictor} under permutation would be a selection-biased test - by chance
#' alone, the best of \code{length(predictors)} candidates will look better
#' than a single fixed candidate's null. Instead, for every permutation this
#' function recomputes the univariate residual R2 for \strong{every} candidate
#' predictor under the same shuffled target and records the \strong{maximum}
#' across predictors; the null distribution is therefore the distribution of
#' "best R2 you would see among this many predictors by chance alone", and the
#' observed top predictor's R2 is compared against that corrected null. This
#' follows the max-statistic / max-T approach to multiple-testing correction
#' (Westfall & Young, 1993).
#'
#' @param data A confident modelling frame.
#' @param target A single target column, e.g. \code{"gPC2"}.
#' @param predictors Candidate predictor columns.
#' @param country,group Confounding group and fold group (see \code{\link{residualize}}).
#' @param k,seeds Fold count and seeds.
#' @param n_perm Number of within-country permutations for the null. Default 200.
#' @param shrink Passed to the residualization (default TRUE, EB-shrunk).
#' @param num.trees,min.node.size,num.random.splits,mtry Ranger hyperparameters used
#'   \strong{identically} for both the observed univariate R2 and every null
#'   permutation, so the null is not biased noisier (or cleaner) than the
#'   observed statistic. Defaults \code{300}, \code{5}, \code{1}, \code{1}.
#'
#' @return An object of class \code{"geaDrivers"} with \code{univariate}
#'   (data frame of per-predictor residual R2, sorted, with \code{correlation}
#'   and \code{sign} columns giving each predictor's marginal direction and, when
#'   \code{length(seeds) >= 2}, \code{R2_CI_low}/\code{R2_CI_high} across-seed 95\%
#'   CI columns), \code{importance} (permutation importance from the joint model),
#'   \code{observed_R2}, \code{observed_R2_CI} (top driver, across-seed),
#'   \code{per_seed_uni} (seeds x predictors univariate R2 matrix, or NULL for a
#'   single seed), \code{correlation}, \code{sign} (named per predictor),
#'   \code{null_max} (vector of per-permutation maxima across all predictors,
#'   the corrected null), \code{null_mean}, \code{null_sd}, \code{p_value}
#'   (against the max-statistic null), \code{verdict_fragility}
#'   ("robust"/"marginal" grade of the driver call) and \code{target}.
#'
#' @references
#' Westfall, P.H. & Young, S.S. (1993). \emph{Resampling-Based Multiple
#'   Testing: Examples and Methods for p-Value Adjustment}. Wiley.
#'
#' @section Expected runtime:
#' Each permutation refits a univariate ExtraTrees model for every candidate
#' predictor (all under fold-safe OOF), so cost scales as
#' \code{n_perm * length(predictors)} small fits. Benchmarked on
#' \code{durumUnits} at the defaults (\code{num.trees=300}, 12 predictors):
#' about 10s per permutation, so \code{n_perm=200} (the default) takes on the
#' order of \strong{30-35 minutes} per target. For interactive exploration,
#' drop to \code{n_perm=20-50} (a few minutes); the package vignette caches a
#' run at lighter settings rather than the full default.
#'
#' @examples
#' \donttest{
#' data(durumUnits)
#' dr <- driverAnalysis(durumUnits, "gPC2", seeds = 42, n_perm = 50)
#' head(dr$univariate)      # top climate drivers of gPC2
#' dr$p_value
#' }
#' @seealso \code{\link{residualize}}
#' @export
driverAnalysis <- function(data, target,
                           predictors = c(paste0("BIO", c(2,5,6,8,9,12,14,15,16,17,19)),
                                          "Altitude"),
                           country = "Country", group = "SiteCode",
                           k = 5, seeds = 42, n_perm = 200, shrink = TRUE,
                           num.trees = 300, min.node.size = 5,
                           num.random.splits = 1, mtry = 1) {
  predictors <- intersect(predictors, names(data))
  stopifnot(length(target) == 1)
  ctry <- as.character(data[[country]]); grp <- data[[group]]
  y <- data[[target]]

  # build fold-safe EB-shrunk residual target once per seed, average
  make_resid <- function(s) {
    folds <- .groupFolds(grp, k = k, seed = s)
    r <- numeric(nrow(data))
    for (f in sort(unique(folds))) {
      te <- which(folds == f); tr <- which(folds != f)
      ytr <- y[tr]; cm <- tapply(ytr, ctry[tr], mean); gm <- mean(ytr)
      if (shrink) {
        nj <- tapply(ytr, ctry[tr], length)
        wv <- tapply(ytr, ctry[tr], function(v) if (length(v) > 1) stats::var(v) else 0)
        s2w <- mean(wv, na.rm = TRUE); s2b <- max(0, stats::var(cm) - s2w / mean(nj))
        B <- s2b / (s2b + s2w / nj); B[is.na(B)] <- 0; cm <- gm + B * (cm - gm)
      }
      r[te] <- y[te] - ifelse(ctry[te] %in% names(cm), cm[ctry[te]], gm)
    }
    r
  }
  resid <- rowMeans(sapply(seeds, make_resid))

  # 1) univariate residual R2 per predictor (spatial-block OOF, first seed)
  #    NOTE: num.trees/min.node.size/splitrule/num.random.splits/mtry are the
  #    same hyperparameters used below for the permutation null (issue #5) -
  #    do not change one without changing the other.
  s <- seeds[1]; folds <- .groupFolds(grp, k = k, seed = s)
  .univariateR2 <- function(y_vec, p) {
    oof <- numeric(nrow(data))
    for (f in sort(unique(folds))) {
      te <- which(folds == f); tr <- which(folds != f)
      df <- data.frame(x = data[[p]][tr], .y = y_vec[tr])
      rf <- ranger::ranger(dependent.variable.name = ".y", data = df,
                           num.trees = num.trees, min.node.size = min.node.size,
                           splitrule = "extratrees", num.random.splits = num.random.splits,
                           mtry = mtry, seed = s, num.threads = 1)
      oof[te] <- stats::predict(rf, data = data.frame(x = data[[p]][te]))$predictions
    }
    ss_tot <- sum((y_vec - mean(y_vec))^2)
    1 - sum((y_vec - oof)^2) / ss_tot
  }
  uni <- sapply(predictors, function(p) .univariateR2(resid, p))

  # Direction of each predictor's marginal effect (issue #19 / spec §8): Pearson
  # correlation between the predictor and the fold-safe residual target, on the
  # residual scale the driver ranking lives on. `sign` is just its sign, exposed
  # as a convenience for the ecological narrative ("warmer origin -> higher gPCk").
  corr <- vapply(predictors, function(p)
    suppressWarnings(stats::cor(data[[p]], resid)), numeric(1))
  sgn <- ifelse(is.na(corr), NA_character_, ifelse(corr >= 0, "positive", "negative"))

  # Per-seed univariate R2 (issue #14): recompute the univariate R2 of the top
  # predictor under each seed's own residual + folds, so a CI can be formed
  # across seeds. Only meaningful with >=2 seeds; single-seed runs get NA.
  uni_df <- data.frame(predictor = predictors, R2 = as.numeric(uni),
                       correlation = as.numeric(corr), sign = sgn,
                       stringsAsFactors = FALSE)
  uni_df <- uni_df[order(-uni_df$R2), ]; rownames(uni_df) <- NULL

  if (length(seeds) >= 2) {
    per_seed_uni <- matrix(NA_real_, length(seeds), length(predictors),
                           dimnames = list(NULL, predictors))
    for (si in seq_along(seeds)) {
      ss <- seeds[si]; rr <- make_resid(ss)
      folds_s <- .groupFolds(grp, k = k, seed = ss)
      for (p in predictors) {
        oof <- numeric(nrow(data))
        for (f in sort(unique(folds_s))) {
          te <- which(folds_s == f); tr <- which(folds_s != f)
          df <- data.frame(x = data[[p]][tr], .y = rr[tr])
          rf <- ranger::ranger(dependent.variable.name = ".y", data = df,
                               num.trees = num.trees, min.node.size = min.node.size,
                               splitrule = "extratrees", num.random.splits = num.random.splits,
                               mtry = mtry, seed = ss, num.threads = 1)
          oof[te] <- stats::predict(rf, data = data.frame(x = data[[p]][te]))$predictions
        }
        ss_tot <- sum((rr - mean(rr))^2)
        per_seed_uni[si, p] <- 1 - sum((rr - oof)^2) / ss_tot
      }
    }
    R2_CI <- apply(per_seed_uni, 2, stats::quantile, probs = c(0.025, 0.975),
                   names = FALSE)
    dimnames(R2_CI) <- list(c("2.5%", "97.5%"), predictors)
    uni_df$R2_CI_low  <- R2_CI[1, uni_df$predictor]
    uni_df$R2_CI_high <- R2_CI[2, uni_df$predictor]
  } else {
    per_seed_uni <- NULL
    uni_df$R2_CI_low <- NA_real_; uni_df$R2_CI_high <- NA_real_
  }

  # 2) joint-model permutation importance
  dfull <- data.frame(data[, predictors, drop = FALSE], .y = resid)
  rf_full <- ranger::ranger(dependent.variable.name = ".y", data = dfull,
                            num.trees = 600, min.node.size = 3,
                            splitrule = "extratrees", num.random.splits = 1,
                            mtry = length(predictors), importance = "permutation",
                            seed = s, num.threads = 1)
  imp <- sort(ranger::importance(rf_full), decreasing = TRUE)

  # 3) within-country permutation significance of the top univariate predictor,
  #    using a MAX-STATISTIC null across ALL candidate predictors (not just the
  #    single predictor that happened to rank first on the real data) to
  #    correct for the multiple-comparisons/selection effect of choosing the
  #    top predictor after seeing all univariate R2s (issue #4). The same
  #    hyperparameters as the observed univariate fit are used throughout
  #    (issue #5), so the null is neither noisier nor cleaner than the signal
  #    it is compared against.
  top <- uni_df$predictor[1]
  obs_r2 <- uni_df$R2[1]
  null_max <- numeric(n_perm)
  for (b in seq_len(n_perm)) {
    set.seed(1000 + b)
    yp <- resid
    for (g in unique(ctry)) { ix <- which(ctry == g); yp[ix] <- sample(yp[ix]) }
    r2_all_preds <- sapply(predictors, function(p) .univariateR2(yp, p))
    null_max[b] <- max(r2_all_preds)
  }
  p_value <- (sum(null_max >= obs_r2) + 1) / (n_perm + 1)

  # CI for the top driver's observed R2 (issue #14): the across-seed quantiles
  # already computed above for that predictor (NA if single-seed).
  observed_R2_CI <- if (length(seeds) >= 2)
    c(uni_df$R2_CI_low[1], uni_df$R2_CI_high[1]) else c(NA_real_, NA_real_)
  names(observed_R2_CI) <- c("2.5%", "97.5%")

  # Fragility flag (issue #15): a floor p-value or a top-driver R2 below the
  # small-effect threshold (0.02, per spec §2) makes the driver call marginal.
  frag <- "robust"
  if (abs(p_value - 1 / (n_perm + 1)) < 1e-10) frag <- "marginal"
  if (obs_r2 < 0.02) frag <- "marginal"

  structure(list(target = target, univariate = uni_df, importance = imp,
                 top_driver = top, observed_R2 = obs_r2,
                 observed_R2_CI = observed_R2_CI, per_seed_uni = per_seed_uni,
                 correlation = stats::setNames(corr, predictors),
                 sign = stats::setNames(sgn, predictors),
                 null_max = null_max,
                 null_mean = mean(null_max), null_sd = stats::sd(null_max),
                 p_value = p_value, verdict_fragility = frag,
                 n_perm = n_perm, n_predictors_tested = length(predictors)),
            class = "geaDrivers")
}

#' @export
#' @method print geaDrivers
print.geaDrivers <- function(x, ...) {
  cat("<geaDrivers> target = ", x$target, "\n", sep = "")
  cat("  Top univariate residual R2 (sign = direction of marginal effect):\n")
  top <- utils::head(x$univariate, 5)
  has_ci <- !is.null(top$R2_CI_low) && any(is.finite(top$R2_CI_low))
  for (i in seq_len(nrow(top))) {
    sgn_i <- if (!is.null(top$sign)) top$sign[i] else NA
    if (has_ci && is.finite(top$R2_CI_low[i])) {
      cat(sprintf("    %-9s %+.4f  [%+.4f, %+.4f]  (%s)\n",
                  top$predictor[i], top$R2[i], top$R2_CI_low[i], top$R2_CI_high[i],
                  sgn_i))
    } else {
      cat(sprintf("    %-9s %+.4f  (%s)\n", top$predictor[i], top$R2[i], sgn_i))
    }
  }
  cat(sprintf("  Top permutation-importance driver: %s\n", names(x$importance)[1]))
  cat(sprintf("  Significance (max-stat null over %d predictors, n_perm=%d):\n",
              x$n_predictors_tested, x$n_perm))
  cat(sprintf("    observed R2=%.4f vs null max %.4f (SD %.4f) -> p=%.4f\n",
              x$observed_R2, x$null_mean, x$null_sd, x$p_value))
  if (!is.null(x$verdict_fragility))
    cat(sprintf("    FRAGILITY: %s\n", x$verdict_fragility))
  invisible(x)
}
