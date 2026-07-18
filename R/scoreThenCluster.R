# Internal: self-contained Quadratic Discriminant Analysis with feature
# standardisation and a small ridge shrinkage of each class covariance toward
# the identity (in the standardised space). Kept internal so the package needs
# no MASS dependency and the regularisation is fully controlled. Not exported.
.qdaFit <- function(X, y, reg = 1e-3) {
  X <- as.matrix(X)
  ctr <- colMeans(X)
  scl <- apply(X, 2, stats::sd); scl[scl == 0 | !is.finite(scl)] <- 1
  Xs <- scale(X, center = ctr, scale = scl)
  classes <- sort(unique(as.character(y)))
  y <- as.character(y)
  params <- lapply(classes, function(cl) {
    Xc <- Xs[y == cl, , drop = FALSE]
    mu <- colMeans(Xc)
    S <- stats::cov(Xc)
    S <- (1 - reg) * S + reg * diag(nrow(S))       # shrink toward identity
    list(mu = mu, prec = solve(S),
         logdet = as.numeric(determinant(S, logarithm = TRUE)$modulus),
         prior = nrow(Xc) / nrow(Xs))
  })
  names(params) <- classes
  list(center = ctr, scale = scl, classes = classes, params = params)
}

# Internal: QDA class posteriors (softmax over the per-class log-likelihoods).
.qdaPosterior <- function(fit, X) {
  Xs <- scale(as.matrix(X), center = fit$center, scale = fit$scale)
  ll <- vapply(fit$params, function(p) {
    d <- sweep(Xs, 2, p$mu)
    quad <- rowSums((d %*% p$prec) * d)
    -0.5 * quad - 0.5 * p$logdet + log(p$prior)
  }, numeric(nrow(Xs)))
  if (is.null(dim(ll))) ll <- matrix(ll, nrow = 1)
  ll <- ll - apply(ll, 1, max)
  e <- exp(ll)
  post <- e / rowSums(e)
  colnames(post) <- fit$classes
  post
}

#' Leakage-Free Score-then-Cluster Genetic Assignment (Stage 2)
#'
#' Predicts the discrete genetic cluster of each accession by a two-stage
#' \emph{Score-then-Cluster} procedure: Stage 1 regresses the continuous genetic
#' PCA scores (gPC1-gPC5) from environment with ExtraTrees, and Stage 2 fits a
#' quadratic discriminant classifier that maps genetic scores to cluster
#' membership. Crucially, the Stage-2 classifier is trained \strong{only on the
#' true genetic scores of the training fold} and is then applied to the
#' \emph{predicted} scores of the held-out fold, so no held-out genetic
#' information ever enters training. The whole two-stage pipeline is wrapped in
#' the same spatial-block folds as \code{\link{spatialBlockCV}}.
#'
#' \strong{Intuition:} clustering directly on true genetic scores tells you the
#' best a perfect regressor could do (the \emph{ceiling}). The honest
#' out-of-fold accuracy is always below that ceiling, and the entire gap is
#' attributable to Stage-1 regression error, not to the classifier. Reporting
#' both makes explicit how much of the assignment skill is genuine environmental
#' signal versus classifier optimism.
#'
#' \strong{Caveat on the ceiling (important):} by default \code{ceiling} is
#' computed \strong{in-sample} - the QDA classifier is fit on the true genetic
#' scores of ALL units and then evaluated on those same units, with no held-out
#' split. With 5 classes separated in only 5 dimensions and light ridge
#' shrinkage (\code{reg = 1e-3}), an in-sample QDA fit can be somewhat optimistic
#' - it is an \emph{upper bound}, not a validated estimate, and the true
#' achievable ceiling under sampling variability could be a little lower. As a
#' result, \code{gap = ceiling - accuracy} is only an \strong{approximation} to
#' "Stage-1 regression error alone"; some of the gap could reflect optimism in
#' the in-sample ceiling itself, not purely Stage-1 error. Set
#' \code{cv_ceiling = TRUE} to instead compute the ceiling out-of-fold (k-fold
#' QDA on the true scores, using the same spatial-block folds as the OOF
#' accuracy), which is slower (an extra k QDA fits) but removes this optimism;
#' the OOF ceiling is always \code{<=} the in-sample ceiling and is the more
#' defensible number to report as "the best a perfect regressor could do".
#'
#' Two calibration diagnostics accompany the accuracy: (i) posterior
#' \emph{entropy}, whose quartiles should show monotonically declining accuracy
#' if the posteriors are well calibrated; and (ii) a \emph{typicality} test that
#' flags accessions whose predicted scores fall outside every cluster's
#' 95\% Mahalanobis ellipsoid (\code{off_manifold}), i.e. predictions that land
#' in no known genetic neighbourhood and should not be trusted.
#'
#' @param data A confident modelling frame (e.g. \code{durumUnits}).
#' @param predictors Character vector of predictor columns. Default BIO1-BIO19 +
#'   Altitude + Latitude + Longitude.
#' @param targets Character vector of genetic-score columns. Default gPC1-gPC5.
#' @param cluster Name of the discrete cluster column to predict. Default
#'   \code{"Cluster"}.
#' @param group Column grouping the spatial-block folds. Default \code{"SiteCode"}.
#' @param k Number of spatial-block folds. Default 5.
#' @param seed Random seed for both the fold assignment and the ExtraTrees
#'   fitter. Default 42.
#' @param reg Ridge shrinkage for the QDA class covariances. Default \code{1e-3}.
#' @param num.trees,min.node.size Passed to the internal ExtraTrees fitter.
#' @param per_tree Logical. If \code{TRUE} (default), Stage-2 posteriors are
#'   averaged over the per-tree Stage-1 predictions (marginalising regression
#'   uncertainty through the classifier, matching the reference pipeline); if
#'   \code{FALSE}, the mean Stage-1 prediction is classified once (faster).
#' @param cv_ceiling Logical. If \code{FALSE} (default), \code{ceiling} is the
#'   fast in-sample QDA-on-true-scores accuracy (see Caveat above). If
#'   \code{TRUE}, the ceiling is instead computed out-of-fold using the same
#'   spatial-block folds as the honest accuracy (an additional QDA fit per
#'   fold - slower, but not in-sample-optimistic).
#' @param fit_final Logical. If \code{TRUE}, additionally fit the Stage-2 QDA
#'   classifier on \strong{all} true scores and return it as
#'   \code{final_classifier} (class \code{"geaQDA"}), the object
#'   \code{\link{predictCluster}} needs to map predicted genetic scores to
#'   cluster membership. This is a deployment convenience only - the reported
#'   \code{accuracy} remains the leakage-free out-of-fold estimate. Default
#'   \code{FALSE}.
#'
#' @return An object of class \code{"geaStage2"}: a list with \code{accuracy}
#'   (leakage-free out-of-fold), \code{ceiling} (QDA on the true scores; in-sample
#'   by default, or out-of-fold if \code{cv_ceiling = TRUE} - see Caveat above),
#'   \code{ceiling_type} (records which), \code{confusion} (a table, rows = true,
#'   cols = predicted), \code{posterior} (out-of-fold class posteriors),
#'   \code{predicted}, \code{entropy}, \code{entropy_calibration} (accuracy by
#'   entropy quartile), \code{state} (per-unit confidence state),
#'   \code{off_manifold}, \code{state_accuracy}, and \code{final_classifier}
#'   (\code{NULL} unless \code{fit_final = TRUE}; otherwise the deployable QDA
#'   classifier fit on all true scores).
#'
#' @examples
#' \donttest{
#' data(durumUnits)
#' s2 <- scoreThenCluster(durumUnits, per_tree = FALSE)  # FALSE = faster
#' s2$accuracy      # honest out-of-fold cluster accuracy
#' s2$ceiling       # accuracy achievable from perfect genetic scores
#' s2$confusion
#' }
#' @references
#' Hurlbert, S.H. (1984). Pseudoreplication and the design of ecological field
#'   experiments. \emph{Ecological Monographs}, 54(2), 187-211.
#'
#' Roberts, D.R. et al. (2017). Cross-validation strategies for data with
#'   temporal, spatial, hierarchical, or phylogenetic structure.
#'   \emph{Ecography}, 40(8), 913-929.
#' @seealso \code{\link{spatialBlockCV}}, \code{\link{locoCV}}
#' @export
scoreThenCluster <- function(data,
                             predictors = c(paste0("BIO", 1:19), "Altitude",
                                            "Latitude", "Longitude"),
                             targets = paste0("gPC", 1:5),
                             cluster = "Cluster", group = "SiteCode",
                             k = 5, seed = 42, reg = 1e-3,
                             num.trees = 600, min.node.size = 3,
                             per_tree = TRUE, cv_ceiling = FALSE,
                             fit_final = FALSE) {
  if (!requireNamespace("ranger", quietly = TRUE))
    stop("Package 'ranger' is required.")
  predictors <- intersect(predictors, names(data))
  Y <- as.matrix(data[, targets, drop = FALSE])
  y_clu <- as.character(data[[cluster]])
  n <- nrow(data)

  # ---- Leakage-free spatial-block folds (shared by ceiling and accuracy) ---
  folds <- .groupFolds(data[[group]], k = k, seed = seed)
  classes <- sort(unique(y_clu))

  # ---- Stage-2 ceiling: QDA on the TRUE scores ------------------------------
  if (cv_ceiling) {
    # Out-of-fold ceiling: same folds as the honest accuracy below, so the
    # ceiling is not in-sample-optimistic (an extra k QDA fits on true scores).
    ceil_pred <- rep(NA_character_, n)
    for (f in sort(unique(folds))) {
      te <- which(folds == f); tr <- which(folds != f)
      if (length(tr) < 5 || length(te) < 1) next
      cf <- .qdaFit(Y[tr, , drop = FALSE], y_clu[tr], reg = reg)
      cp <- .qdaPosterior(cf, Y[te, , drop = FALSE])
      ceil_pred[te] <- cf$classes[max.col(cp, ties.method = "first")]
    }
    ok <- !is.na(ceil_pred)
    ceiling <- mean(ceil_pred[ok] == y_clu[ok])
    ceiling_type <- "out-of-fold (spatial-block, leakage-free)"
  } else {
    ceil_fit <- .qdaFit(Y, y_clu, reg = reg)
    ceil_post <- .qdaPosterior(ceil_fit, Y)
    ceil_pred <- ceil_fit$classes[max.col(ceil_post, ties.method = "first")]
    ceiling <- mean(ceil_pred == y_clu)
    ceiling_type <- "in-sample (optimistic upper bound - see Caveat in ?scoreThenCluster)"
  }

  # ---- Leakage-free spatial-block OOF (Stage-1 -> Stage-2 accuracy) --------
  post <- matrix(0, n, length(classes), dimnames = list(NULL, classes))
  for (f in sort(unique(folds))) {
    te <- which(folds == f); tr <- which(folds != f)
    if (length(tr) < 5 || length(te) < 1) next
    qda <- .qdaFit(Y[tr, , drop = FALSE], y_clu[tr], reg = reg)  # TRUE scores only
    if (per_tree) {
      # per-target per-tree Stage-1 predictions -> average QDA posteriors
      tree_pred <- vector("list", length(targets))
      for (jt in seq_along(targets)) {
        df <- data.frame(data[tr, predictors, drop = FALSE], .y = Y[tr, jt])
        rf <- ranger::ranger(dependent.variable.name = ".y", data = df,
                             num.trees = num.trees, min.node.size = min.node.size,
                             splitrule = "extratrees", num.random.splits = 1,
                             mtry = length(predictors), seed = seed, num.threads = 1)
        tree_pred[[jt]] <- stats::predict(rf, data = data[te, predictors, drop = FALSE],
                                          predict.all = TRUE)$predictions  # n_te x n_trees
      }
      nt <- ncol(tree_pred[[1]])
      acc_post <- matrix(0, length(te), length(classes))
      for (b in seq_len(nt)) {
        Sb <- sapply(tree_pred, function(m) m[, b])          # n_te x n_targets
        acc_post <- acc_post + .qdaPosterior(qda, Sb)
      }
      post[te, ] <- acc_post / nt
    } else {
      pred_gpc <- .fitPredictMO(data[tr, ], data[te, ], predictors, targets,
                                num.trees = num.trees, min.node.size = min.node.size,
                                seed = seed)
      post[te, ] <- .qdaPosterior(qda, pred_gpc)
    }
  }
  pred <- classes[max.col(post, ties.method = "first")]
  accuracy <- mean(pred == y_clu)
  confusion <- table(true = y_clu, predicted = factor(pred, levels = classes))

  # ---- Entropy calibration --------------------------------------------------
  H <- -rowSums(ifelse(post > 0, post * log(post), 0))
  correct <- as.integer(pred == y_clu)
  qs <- stats::quantile(H, probs = seq(0, 1, 0.25), na.rm = TRUE)
  qs[1] <- -Inf; qs[length(qs)] <- Inf
  hq <- cut(H, breaks = unique(qs), labels = FALSE, include.lowest = TRUE)
  entropy_calibration <- data.frame(
    quartile = paste0("Q", sort(unique(hq))),
    accuracy = tapply(correct, hq, mean)[order(unique(sort(hq)))],
    n = as.integer(tapply(correct, hq, length)[order(unique(sort(hq)))]),
    row.names = NULL)

  # ---- Typicality: off-manifold via per-cluster Mahalanobis on OOF scores ---
  oof_pred_gpc <- if (per_tree) {
    # reconstruct mean OOF prediction for the Mahalanobis check
    op <- matrix(NA_real_, n, length(targets), dimnames = list(NULL, targets))
    for (f in sort(unique(folds))) {
      te <- which(folds == f); tr <- which(folds != f)
      if (length(tr) < 5 || length(te) < 1) next
      op[te, ] <- .fitPredictMO(data[tr, ], data[te, ], predictors, targets,
                                num.trees = num.trees, min.node.size = min.node.size,
                                seed = seed)
    }
    op
  } else post   # placeholder; replaced below if !per_tree
  if (!per_tree) {
    op <- matrix(NA_real_, n, length(targets), dimnames = list(NULL, targets))
    for (f in sort(unique(folds))) {
      te <- which(folds == f); tr <- which(folds != f)
      if (length(tr) < 5 || length(te) < 1) next
      op[te, ] <- .fitPredictMO(data[tr, ], data[te, ], predictors, targets,
                                num.trees = num.trees, min.node.size = min.node.size,
                                seed = seed)
    }
    oof_pred_gpc <- op
  }
  thr <- stats::qchisq(0.95, df = length(targets))
  mahal_min <- rep(Inf, n)
  for (cl in classes) {
    m <- y_clu == cl
    mu <- colMeans(Y[m, , drop = FALSE])
    S <- stats::cov(Y[m, , drop = FALSE])
    md <- stats::mahalanobis(oof_pred_gpc, center = mu, cov = S, tol = 1e-20)
    mahal_min <- pmin(mahal_min, md)
  }
  off_manifold <- mahal_min > thr
  medH <- stats::median(H, na.rm = TRUE)
  state <- ifelse(off_manifold, "off_manifold",
           ifelse(H > medH, "fuzzy_boundary", "confident_in_cluster"))
  state_accuracy <- data.frame(
    state = names(tapply(correct, state, mean)),
    accuracy = as.numeric(tapply(correct, state, mean)),
    n = as.integer(tapply(correct, state, length)),
    row.names = NULL)

  # ---- Optional deployable classifier (issue #21) ---------------------------
  # The QDA fit on ALL true scores is the object predictCluster() needs to map
  # predicted genetic scores -> cluster membership. It is deployment-only: the
  # honest accuracy above is still the leakage-free out-of-fold number. When
  # cv_ceiling = FALSE we already fit exactly this object (ceil_fit); reuse it.
  final_classifier <- NULL
  if (isTRUE(fit_final)) {
    final_classifier <- if (!cv_ceiling && exists("ceil_fit"))
      ceil_fit else .qdaFit(Y, y_clu, reg = reg)
    final_classifier$targets <- targets
    final_classifier$reg <- reg
    class(final_classifier) <- "geaQDA"
  }

  structure(list(
    accuracy = accuracy, ceiling = ceiling, ceiling_type = ceiling_type,
    gap = ceiling - accuracy,
    confusion = confusion, posterior = post, predicted = pred,
    entropy = H, entropy_calibration = entropy_calibration,
    off_manifold = off_manifold, state = state,
    state_accuracy = state_accuracy, classes = classes,
    per_tree = per_tree, final_classifier = final_classifier,
    type = "Score-then-Cluster (leakage-free, spatial-block OOF)"),
    class = "geaStage2")
}

#' @export
#' @method print geaStage2
print.geaStage2 <- function(x, ...) {
  cat("<geaStage2> ", x$type, "\n", sep = "")
  cat(sprintf("  Out-of-fold accuracy : %.4f\n", x$accuracy))
  cat(sprintf("  Ceiling (true scores): %.4f  [%s]\n", x$ceiling, x$ceiling_type))
  cat(sprintf("  Gap (Stage-1 error)  : %.4f%s\n", x$gap,
              if (identical(x$ceiling_type, "in-sample (optimistic upper bound - see Caveat in ?scoreThenCluster)"))
                "  (approx.; ceiling is in-sample, gap may include ceiling optimism)" else ""))
  cat("  Confusion (rows=true, cols=predicted):\n")
  print(x$confusion)
  cat("  Accuracy by entropy quartile (monotone decline = calibrated):\n")
  ec <- x$entropy_calibration
  for (i in seq_len(nrow(ec)))
    cat(sprintf("    %-4s acc=%.3f  n=%d\n", ec$quartile[i], ec$accuracy[i], ec$n[i]))
  cat(sprintf("  off_manifold units: %d / %d\n", sum(x$off_manifold), length(x$off_manifold)))
  invisible(x)
}
