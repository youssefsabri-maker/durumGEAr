#' Spatial-Block Cross-Validation (Interpolation Skill)
#'
#' Evaluates multi-output regression skill under \emph{grouped} k-fold
#' cross-validation in which every whole spatial group (by default the collection
#' site) is kept entirely on one side of each train/test split. This is the
#' honest interpolation estimate: it prevents near-duplicate rows from the same
#' site leaking across the split, which is the standard mechanism by which
#' spatial autocorrelation inflates naive random-CV scores (Roberts et al., 2017).
#'
#' \strong{Intuition:} random CV can place a site's twin rows in both train and
#' test, so the model is quietly graded on data it has already seen. Grouping the
#' folds by site forces prediction of genuinely unseen locations.
#'
#' @param data A confident modelling frame (e.g. \code{durumUnits}).
#' @param predictors Character vector of predictor columns. Default BIO1-BIO19 +
#'   Altitude + Latitude + Longitude.
#' @param targets Character vector of target columns. Default gPC1-gPC5.
#' @param group Column name to group folds by. Default \code{"SiteCode"}.
#' @param k Number of folds. Default 5.
#' @param seeds Integer vector of random seeds; results are averaged over seeds.
#'   Default \code{42:46}.
#' @param num.trees,min.node.size Passed to the internal ExtraTrees fitter.
#'
#' @return An object of class \code{"geaCV"}: a list with \code{oof} (out-of-fold
#'   predictions for the first seed), \code{metrics} (a \code{geaMetrics} object
#'   averaged across seeds), \code{per_seed_mean_R2}, and \code{type}.
#'
#' @examples
#' \donttest{
#' data(durumUnits)
#' cv <- spatialBlockCV(durumUnits, seeds = 42)   # single seed for speed
#' cv$metrics
#' }
#'
#' @references
#' Roberts, D.R. et al. (2017). Cross-validation strategies for data with
#'   temporal, spatial, hierarchical, or phylogenetic structure.
#'   \emph{Ecography}, 40(8), 913-929.
#' @seealso \code{\link{locoCV}}, \code{\link{getMetrics}}
#' @export
spatialBlockCV <- function(data,
                           predictors = c(paste0("BIO", 1:19), "Altitude",
                                          "Latitude", "Longitude"),
                           targets = paste0("gPC", 1:5),
                           group = "SiteCode", k = 5, seeds = 42:46,
                           num.trees = 600, min.node.size = 3) {
  predictors <- intersect(predictors, names(data))
  Y <- as.matrix(data[, targets, drop = FALSE])
  grp <- data[[group]]
  per_seed_R2 <- numeric(length(seeds))
  first_oof <- NULL
  per_target_acc <- matrix(0, length(seeds), length(targets),
                           dimnames = list(NULL, targets))
  for (si in seq_along(seeds)) {
    s <- seeds[si]
    folds <- .groupFolds(grp, k = k, seed = s)
    oof <- matrix(NA_real_, nrow(data), length(targets), dimnames = list(NULL, targets))
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
                 type = "spatial-block (grouped by site) - interpolation"),
            class = "geaCV")
}

#' @export
#' @method print geaCV
print.geaCV <- function(x, ...) {
  cat("<geaCV> ", x$type, "\n", sep = "")
  cat(sprintf("  Per-target-mean R2: %+.4f", x$metrics$mean_R2))
  if (length(x$per_seed_mean_R2) > 1)
    cat(sprintf("  (SD %.4f over %d seeds)", stats::sd(x$per_seed_mean_R2),
                length(x$per_seed_mean_R2)))
  cat("\n")
  pt <- x$metrics$per_target
  cat("  ", paste(sprintf("%s=%+.3f", pt$target, pt$R2), collapse = "  "), "\n", sep = "")
  invisible(x)
}
