# ============================================================================
# Deployment layer (issue #20): fit a per-target predictor on all confident
# data and score new observations, inheriting the fragility/verdict gates that
# residualize() produced during validation. See ?fitGeneticScoreModel.
# ============================================================================

# Internal: derive a per-target reliability label from a residualize() result.
# The gate uses BOTH fields residualize() exposes: a target is "gated" (its
# predictions should be treated as unreliable) when its verdict is a pure
# country-identity artifact OR its verdict is shaky (fragile / marginal).
# "artifact" is NOT a verdict_fragility grade (those are robust/marginal/
# fragile) - it lives in $verdict as "country-identity artifact".
.reliabilityFromResidualize <- function(residualize_result, targets) {
  out <- data.frame(target = targets,
                    verdict = NA_character_, fragility = NA_character_,
                    residual_R2 = NA_real_, gated = FALSE,
                    stringsAsFactors = FALSE)
  if (is.null(residualize_result)) {
    out$reliability <- "ungated (no residualize_result supplied)"
    return(out)
  }
  v  <- residualize_result$verdict
  vf <- residualize_result$verdict_fragility
  r2 <- residualize_result$per_target_R2
  for (i in seq_along(targets)) {
    tg <- targets[i]
    out$verdict[i]     <- if (tg %in% names(v))  v[[tg]]  else NA_character_
    out$fragility[i]   <- if (tg %in% names(vf)) vf[[tg]] else NA_character_
    out$residual_R2[i] <- if (tg %in% names(r2)) r2[[tg]] else NA_real_
    is_artifact <- isTRUE(out$verdict[i] == "country-identity artifact")
    is_shaky    <- isTRUE(out$fragility[i] %in% c("fragile", "marginal"))
    out$gated[i] <- is_artifact || is_shaky
  }
  out$reliability <- ifelse(out$gated, "gated (do not trust)", "trusted")
  out
}

#' Fit a deployable genetic-score predictor
#'
#' Trains one ExtraTrees (\code{ranger}) regressor per genetic PC on
#' \strong{all} supplied rows - no cross-validation split, because validation
#' already happened in \code{\link{spatialBlockCV}} / \code{\link{residualize}}.
#' The returned model carries per-target reliability labels inherited from a
#' \code{\link{residualize}} result so that \code{\link{predict.geneticScoreModel}}
#' can flag which axes are trustworthy out of sample.
#'
#' @details
#' \strong{This is a deployment convenience, not a new validation.} The honest
#' out-of-sample skill of these predictions is whatever \code{spatialBlockCV()}
#' and \code{locoCV()} already reported: near-zero transferable skill for most
#' genetic axes across countries, with genuine within-domain signal essentially
#' only on \code{gPC2} (residual R2 ~ 0.18) and weakly \code{gPC5} (~ 0.06).
#' Passing a \code{residualize_result} makes the model label the remaining axes
#' as gated so downstream calls refuse to present them as reliable.
#'
#' @param data A data frame of training units (e.g. \code{durumUnits}) with the
#'   predictor and target columns.
#' @param targets Character vector of target columns. Default \code{gPC1:gPC5}.
#' @param predictors Character vector of predictor columns. Defaults to the
#'   bioclim + altitude set actually used in validation; any not present in
#'   \code{data} are dropped.
#' @param num.trees,min.node.size ranger hyperparameters, matching the values
#'   used during validation.
#' @param residualize_result Optional object from \code{\link{residualize}}.
#'   When supplied, each target inherits a reliability label: a target is
#'   \strong{gated} when its residualization verdict is
#'   \code{"country-identity artifact"} OR its fragility grade is
#'   \code{"fragile"} / \code{"marginal"}.
#' @param seed Integer seed passed to ranger for reproducibility.
#'
#' @return An object of class \code{"geneticScoreModel"}: a list with
#'   \code{models} (one ranger fit per target), \code{targets},
#'   \code{predictors}, \code{reliability} (a data frame, one row per target),
#'   \code{train_predictors} (the training predictor matrix, retained for
#'   \code{\link{checkExtrapolationRisk}}), and \code{hyper}.
#'
#' @examples
#' \donttest{
#' data(durumUnits)
#' rz  <- residualize(durumUnits, seeds = 42, n_perm = 20, num.trees = 100)
#' mod <- fitGeneticScoreModel(durumUnits, residualize_result = rz)
#' mod
#' }
#' @seealso \code{\link{predict.geneticScoreModel}},
#'   \code{\link{checkExtrapolationRisk}}, \code{\link{robustGeneticScore}}
#' @aliases geneticScoreModel
#' @export
fitGeneticScoreModel <- function(data,
                                 targets = paste0("gPC", 1:5),
                                 predictors = c(paste0("BIO", 1:19), "Altitude"),
                                 num.trees = 300, min.node.size = 3,
                                 residualize_result = NULL, seed = 42) {
  if (!requireNamespace("ranger", quietly = TRUE))
    stop("Package 'ranger' is required.")
  predictors <- intersect(predictors, names(data))
  if (length(predictors) == 0L) stop("No predictors found in `data`.")
  miss <- setdiff(targets, names(data))
  if (length(miss)) stop("Targets not found in `data`: ",
                         paste(miss, collapse = ", "))

  models <- stats::setNames(vector("list", length(targets)), targets)
  for (t in targets) {
    df <- data.frame(data[, predictors, drop = FALSE], .y = data[[t]])
    models[[t]] <- ranger::ranger(
      dependent.variable.name = ".y", data = df,
      num.trees = num.trees, min.node.size = min.node.size,
      splitrule = "extratrees", num.random.splits = 1,
      mtry = length(predictors), respect.unordered.factors = "order",
      seed = seed, num.threads = 1)
  }

  reliability <- .reliabilityFromResidualize(residualize_result, targets)

  structure(list(
    models = models, targets = targets, predictors = predictors,
    reliability = reliability,
    train_predictors = as.matrix(data[, predictors, drop = FALSE]),
    hyper = list(num.trees = num.trees, min.node.size = min.node.size,
                 seed = seed),
    type = "per-target ExtraTrees genetic-score predictor"),
    class = "geneticScoreModel")
}

#' @export
#' @method print geneticScoreModel
print.geneticScoreModel <- function(x, ...) {
  cat("<geneticScoreModel> ", x$type, "\n", sep = "")
  cat(sprintf("  targets    : %s\n", paste(x$targets, collapse = ", ")))
  cat(sprintf("  predictors : %d (%s...)\n", length(x$predictors),
              paste(utils::head(x$predictors, 3), collapse = ", ")))
  cat(sprintf("  trees      : %d, min.node.size = %d\n",
              x$hyper$num.trees, x$hyper$min.node.size))
  rl <- x$reliability
  cat("  reliability (from residualize):\n")
  for (i in seq_len(nrow(rl)))
    cat(sprintf("    %-5s %-20s frag=%-8s R2=%s -> %s\n",
                rl$target[i],
                ifelse(is.na(rl$verdict[i]), "-", rl$verdict[i]),
                ifelse(is.na(rl$fragility[i]), "-", rl$fragility[i]),
                ifelse(is.na(rl$residual_R2[i]), "  NA",
                       sprintf("%+.3f", rl$residual_R2[i])),
                rl$reliability[i]))
  trusted <- rl$target[!rl$gated]
  cat(sprintf("  trusted axes: %s\n",
              if (length(trusted)) paste(trusted, collapse = ", ") else "(none)"))
  invisible(x)
}

#' Predict genetic scores on new observations
#'
#' @param object A \code{\link{geneticScoreModel}}.
#' @param newdata A data frame containing the model's predictor columns.
#' @param gated_to_na Logical. If \code{TRUE}, columns for gated (untrusted)
#'   targets are returned as \code{NA} rather than silently presenting numbers
#'   the validation does not support. Default \code{FALSE} (return all columns,
#'   but attach reliability metadata).
#' @param ... Unused.
#'
#' @return A data frame with one column per target. Attribute
#'   \code{"reliability"} carries the per-target label table and attribute
#'   \code{"gated_targets"} the character vector of gated targets. The reported
#'   numbers are point predictions only - they are \strong{not} accompanied by
#'   a formal predictive confidence interval.
#' @seealso \code{\link{fitGeneticScoreModel}}, \code{\link{robustGeneticScore}}
#' @export
#' @method predict geneticScoreModel
predict.geneticScoreModel <- function(object, newdata, gated_to_na = FALSE, ...) {
  miss <- setdiff(object$predictors, names(newdata))
  if (length(miss)) stop("`newdata` is missing predictor columns: ",
                         paste(miss, collapse = ", "))
  out <- matrix(NA_real_, nrow(newdata), length(object$targets),
                dimnames = list(NULL, object$targets))
  for (t in object$targets)
    out[, t] <- stats::predict(
      object$models[[t]],
      data = newdata[, object$predictors, drop = FALSE])$predictions
  out <- as.data.frame(out)
  gated <- object$reliability$target[object$reliability$gated]
  if (isTRUE(gated_to_na) && length(gated))
    out[, gated] <- NA_real_
  attr(out, "reliability") <- object$reliability
  attr(out, "gated_targets") <- gated
  out
}

#' Flag out-of-domain observations before trusting a prediction
#'
#' The validation shows transferable climate-to-genetics skill is near zero
#' outside the training domain, so predictions on climate combinations far from
#' the training distribution are extrapolation and should be refused. This
#' function flags such rows two complementary ways - neither is a convex hull,
#' which is degenerate in the ~20-dimensional predictor space:
#' \enumerate{
#'   \item \strong{Per-variable range} (range-based applicability domain): a row
#'     is out-of-range if any predictor falls outside the training
#'     \code{[q_lo, q_hi]} band. With the default \code{envelope_prob = 1}
#'     this is the training min/max, so no training row is ever out-of-range;
#'     tighten it (e.g. 0.99) to also flag rows in the extreme tails.
#'   \item \strong{Multivariate distance}: the Mahalanobis distance of the row
#'     to the training predictor mean/covariance, flagged against the
#'     \strong{empirical} \code{mahalanobis_prob} quantile of the training
#'     distances. Using the empirical quantile (rather than a chi-square
#'     approximation that assumes multivariate normality) keeps the in-domain
#'     false-flag rate near \code{1 - mahalanobis_prob} on the real,
#'     collinear climate data. This catches unusual \emph{combinations} of
#'     otherwise in-range values.
#' }
#' A row is \code{extrapolation} if it trips either test.
#'
#' @param object A \code{\link{geneticScoreModel}} (carries the training
#'   predictor matrix), or a numeric matrix / data frame of training predictors.
#' @param newdata Data frame containing the predictor columns to check.
#' @param predictors Character vector of predictor columns. Defaults to those
#'   stored on \code{object} (required if \code{object} is a bare matrix).
#' @param envelope_prob Central probability mass of the per-variable band.
#'   Default 1 (training min/max). Set below 1 (e.g. 0.99) to also flag tails.
#' @param mahalanobis_prob Empirical training-distance quantile used as the
#'   multivariate threshold. Default 0.975.
#'
#' @return A data frame with one row per \code{newdata} row: \code{mahalanobis}
#'   (distance), \code{md_flag} (logical), \code{n_out_of_range} (count of
#'   predictors outside the band), \code{range_flag} (logical), and
#'   \code{extrapolation} (either flag). The distance threshold and band
#'   probability are returned as attributes \code{"md_threshold"} /
#'   \code{"envelope_prob"}.
#' @examples
#' \donttest{
#' data(durumUnits)
#' mod <- fitGeneticScoreModel(durumUnits)
#' risk <- checkExtrapolationRisk(mod, durumUnits[1:10, ])
#' risk$extrapolation
#' }
#' @seealso \code{\link{fitGeneticScoreModel}}, \code{\link{predict.geneticScoreModel}}
#' @export
checkExtrapolationRisk <- function(object, newdata, predictors = NULL,
                                   envelope_prob = 1,
                                   mahalanobis_prob = 0.975) {
  if (inherits(object, "geneticScoreModel")) {
    train <- object$train_predictors
    if (is.null(predictors)) predictors <- object$predictors
  } else {
    train <- as.matrix(object)
    if (is.null(predictors)) predictors <- colnames(train)
    if (is.null(predictors)) stop("Provide `predictors` when `object` is a bare matrix.")
  }
  predictors <- intersect(predictors, colnames(train))
  miss <- setdiff(predictors, names(newdata))
  if (length(miss)) stop("`newdata` is missing predictor columns: ",
                         paste(miss, collapse = ", "))
  train <- train[, predictors, drop = FALSE]
  X <- as.matrix(newdata[, predictors, drop = FALSE])

  # ---- Per-variable quantile envelope --------------------------------------
  tail_p <- (1 - envelope_prob) / 2
  q_lo <- apply(train, 2, stats::quantile, probs = tail_p,     na.rm = TRUE)
  q_hi <- apply(train, 2, stats::quantile, probs = 1 - tail_p, na.rm = TRUE)
  below <- sweep(X, 2, q_lo, "<")
  above <- sweep(X, 2, q_hi, ">")
  n_out <- rowSums(below | above)
  range_flag <- n_out > 0

  # ---- Multivariate Mahalanobis distance -----------------------------------
  mu <- colMeans(train)
  S  <- stats::cov(train)
  # Ridge-regularize in case a predictor is (near-)collinear at this scale.
  S  <- S + diag(1e-8 * mean(diag(S)), ncol(S))
  md       <- stats::mahalanobis(X,     center = mu, cov = S)
  md_train <- stats::mahalanobis(train, center = mu, cov = S)
  # Empirical threshold: the mahalanobis_prob quantile of TRAINING distances,
  # so the in-domain false-flag rate is ~ 1 - mahalanobis_prob regardless of
  # how non-normal / collinear the real climate data is (a chi-square quantile
  # would badly over-flag here).
  thr <- stats::quantile(md_train, probs = mahalanobis_prob, na.rm = TRUE)
  md_flag <- md > thr

  out <- data.frame(mahalanobis = md, md_flag = md_flag,
                    n_out_of_range = n_out, range_flag = range_flag,
                    extrapolation = md_flag | range_flag)
  attr(out, "md_threshold") <- thr
  attr(out, "envelope_prob") <- envelope_prob
  out
}

#' Predict cluster membership from a fitted genetic-score model
#'
#' Chains \code{\link{predict.geneticScoreModel}} (climate -> genetic scores)
#' into the deployable QDA classifier (genetic scores -> cluster) produced by
#' \code{scoreThenCluster(..., fit_final = TRUE)}. The two stages must have been
#' trained on the same target set.
#'
#' \strong{Honesty note.} The realistic accuracy of this end-to-end path is the
#' leakage-free \code{accuracy} that \code{\link{scoreThenCluster}} reported -
#' bounded well below the \code{ceiling}, because Stage-1 genetic scores are
#' predicted from climate with limited transferable skill. Treat the returned
#' cluster as a best guess, weighted by its posterior probability, not a
#' certainty.
#'
#' @param object A \code{\link{geneticScoreModel}}.
#' @param newdata Data frame with the model's predictor columns.
#' @param classifier A \code{"geaQDA"} object, i.e. the \code{final_classifier}
#'   element of a \code{scoreThenCluster(..., fit_final = TRUE)} result.
#' @param ... Unused.
#'
#' @return A data frame with \code{predicted} (cluster label), \code{posterior}
#'   (probability of the assigned cluster) and \code{entropy} (of the full
#'   posterior, a per-row uncertainty measure). The full posterior matrix is
#'   attached as attribute \code{"posterior"}.
#' @examples
#' \donttest{
#' data(durumUnits)
#' mod <- fitGeneticScoreModel(durumUnits)
#' s2  <- scoreThenCluster(durumUnits, per_tree = FALSE, fit_final = TRUE)
#' pc  <- predictCluster(mod, durumUnits[1:10, ], s2$final_classifier)
#' pc$predicted
#' }
#' @seealso \code{\link{scoreThenCluster}}, \code{\link{fitGeneticScoreModel}}
#' @export
predictCluster <- function(object, newdata, classifier, ...) {
  if (!inherits(object, "geneticScoreModel"))
    stop("`object` must be a geneticScoreModel.")
  if (!inherits(classifier, "geaQDA"))
    stop("`classifier` must be a 'geaQDA' object from ",
         "scoreThenCluster(..., fit_final = TRUE)$final_classifier.")
  if (!setequal(object$targets, classifier$targets))
    stop("Model targets and classifier targets differ; they must be trained ",
         "on the same genetic PCs.")
  scores <- stats::predict(object, newdata)          # data frame, per-target
  S <- as.matrix(scores[, classifier$targets, drop = FALSE])
  post <- .qdaPosterior(classifier, S)               # n x n_classes
  colnames(post) <- classifier$classes
  idx <- max.col(post, ties.method = "first")
  H <- -rowSums(ifelse(post > 0, post * log(post), 0))
  out <- data.frame(predicted = classifier$classes[idx],
                    posterior = post[cbind(seq_len(nrow(post)), idx)],
                    entropy = H, stringsAsFactors = FALSE)
  attr(out, "posterior") <- post
  out
}

#' Collapse trustworthy genetic axes into a single robust score
#'
#' Convenience wrapper that drops the untrustworthy (gated) genetic axes and
#' collapses what remains into one number per observation. By default the axes
#' to drop are taken from the model's own reliability metadata (inherited from
#' \code{\link{residualize}}), \strong{not} hardcoded - so if a re-run at
#' production \code{n_perm} promotes an axis from gated to trusted, this follows
#' automatically.
#'
#' @param object A \code{\link{geneticScoreModel}}.
#' @param newdata Data frame with the model's predictor columns.
#' @param method How to collapse the retained axes: \code{"mean"} (default),
#'   \code{"median"}, or \code{"pca1"} (first principal component score of the
#'   retained-axis predictions, sign-aligned to their mean).
#' @param keep Optional character vector of axes to force-keep, overriding the
#'   reliability-based default. Use \code{keep = object$targets} to collapse all
#'   axes regardless of gating.
#' @param ... Unused.
#'
#' @return A numeric vector, one robust score per row of \code{newdata}.
#'   Attributes \code{"axes"} (which targets were combined) and \code{"method"}.
#'   Errors if no axis survives gating and none is forced via \code{keep}.
#' @examples
#' \donttest{
#' data(durumUnits)
#' rz  <- residualize(durumUnits, seeds = 42, n_perm = 20, num.trees = 100)
#' mod <- fitGeneticScoreModel(durumUnits, residualize_result = rz)
#' robustGeneticScore(mod, durumUnits[1:10, ])
#' }
#' @seealso \code{\link{fitGeneticScoreModel}}, \code{\link{predict.geneticScoreModel}}
#' @export
robustGeneticScore <- function(object, newdata, method = c("mean", "median", "pca1"),
                               keep = NULL, ...) {
  method <- match.arg(method)
  if (!inherits(object, "geneticScoreModel"))
    stop("`object` must be a geneticScoreModel.")
  if (is.null(keep)) {
    rl <- object$reliability
    keep <- rl$target[!rl$gated]
  }
  keep <- intersect(keep, object$targets)
  if (length(keep) == 0L)
    stop("No trustworthy axis survives gating. Re-run residualize() at ",
         "production n_perm, or force axes with `keep=`.")
  scores <- stats::predict(object, newdata)
  M <- as.matrix(scores[, keep, drop = FALSE])
  val <- switch(method,
    mean   = rowMeans(M),
    median = apply(M, 1, stats::median),
    pca1   = {
      if (ncol(M) == 1L) M[, 1] else {
        pc <- stats::prcomp(M, center = TRUE, scale. = TRUE)
        s  <- pc$x[, 1]
        # sign-align to the mean so the score's direction is interpretable
        if (stats::cor(s, rowMeans(M)) < 0) s <- -s
        s
      }
    })
  attr(val, "axes") <- keep
  attr(val, "method") <- method
  val
}
