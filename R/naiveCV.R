#' Naive (Ungrouped) Random K-Fold Cross-Validation (Leaky Baseline)
#'
#' Evaluates the same multi-output ExtraTrees model as \code{\link{spatialBlockCV}}
#' but with \strong{ordinary random k-fold assignment} - individual rows, not whole
#' spatial groups, are randomly split into folds. Because many rows share an
#' identical collection site (and therefore identical bioclimatic predictors and
#' near-duplicate genetic targets), this naive scheme routinely places
#' near-duplicate rows in both the training and test fold, letting the model be
#' graded on data it has effectively already seen.
#'
#' \strong{Purpose:} this function exists only as a deliberately leaky reference
#' point. Comparing its R2 against \code{\link{spatialBlockCV}} (grouped by site)
#' and \code{\link{locoCV}} (grouped by country) makes the inflation from
#' pseudoreplication a visible, quantified number rather than an assertion: the
#' gap \code{naiveCV - spatialBlockCV} is the pseudoreplication inflation, and
#' the gap \code{spatialBlockCV - locoCV} is the geographic-confounding
#' inflation (see the package vignette).
#'
#' @inheritParams spatialBlockCV
#' @param k Number of folds. Default 5.
#' @param seeds Integer vector of random seeds; results are averaged over seeds.
#'   Default \code{42:46}.
#'
#' @return An object of class \code{"geaCV"} (same shape as
#'   \code{\link{spatialBlockCV}}), with \code{type} recording that the folds
#'   were \strong{not} grouped.
#'
#' @examples
#' \donttest{
#' data(durumUnits)
#' nv <- naiveCV(durumUnits, seeds = 42)
#' nv$metrics$mean_R2      # inflated by pseudoreplication - not the honest estimate
#' }
#' @seealso \code{\link{spatialBlockCV}}, \code{\link{locoCV}}
#' @export
naiveCV <- function(data,
                    predictors = c(paste0("BIO", 1:19), "Altitude",
                                   "Latitude", "Longitude"),
                    targets = paste0("gPC", 1:5),
                    k = 5, seeds = 42:46,
                    num.trees = 600, min.node.size = 3) {
  predictors <- intersect(predictors, names(data))
  Y <- as.matrix(data[, targets, drop = FALSE])
  n <- nrow(data)
  per_seed_R2 <- numeric(length(seeds))
  first_oof <- NULL
  per_target_acc <- matrix(0, length(seeds), length(targets),
                           dimnames = list(NULL, targets))
  for (si in seq_along(seeds)) {
    s <- seeds[si]
    set.seed(s)
    # ordinary random assignment of INDIVIDUAL ROWS to folds - no grouping
    folds <- sample(rep(seq_len(k), length.out = n))
    oof <- matrix(NA_real_, n, length(targets), dimnames = list(NULL, targets))
    for (f in sort(unique(folds))) {
      te <- which(folds == f); tr <- which(folds != f)
      oof[te, ] <- .fitPredictMO(data[tr, ], data[te, ], predictors, targets,
                                 num.trees = num.trees, min.node.size = min.node.size,
                                 seed = s)
    }
    m <- getMetrics(Y, oof)
    per_seed_R2[si] <- m$mean_R2
    per_target_acc[si, ] <- m$per_target$R2
    if (si == 1) first_oof <- oof
  }
  mean_pt <- colMeans(per_target_acc)
  metrics <- structure(list(
    per_target = data.frame(target = targets, R2 = mean_pt,
                            RMSE = NA_real_, stringsAsFactors = FALSE),
    mean_R2 = mean(per_seed_R2), mean_RMSE = NA_real_), class = "geaMetrics")
  structure(list(oof = first_oof, metrics = metrics,
                 per_seed_mean_R2 = per_seed_R2, seeds = seeds,
                 type = "naive random k-fold (NOT grouped by site) - leaky baseline"),
            class = "geaCV")
}
