#' Fit an Environmental PCA (ePC) on a Confident Modelling Frame
#'
#' Standardizes the bioclimatic predictors and fits a principal-component
#' rotation, returning both the scores and a re-usable projector so that the
#' \emph{same} scaling and rotation can be applied to held-out data without
#' leakage.
#'
#' \strong{Intuition:} the 19 bioclimatic variables are heavily redundant
#' (temperatures move together, rainfall variables move together). PCA rotates
#' them into a few uncorrelated axes that carry the same information more
#' compactly and stabilize downstream models.
#'
#' Note that in the durum workflow the PCA rotation is a dimensionality-reduction
#' convenience only: models built on ePC scores and on the raw BIO variables
#' achieve near-identical skill (see the package vignette), so the rotation adds
#' no predictive power - it only aids interpretation and conditioning.
#'
#' @param data A data frame containing the climate predictors.
#' @param vars Character vector of predictor columns. Default BIO1-BIO19 + Altitude.
#' @param n_comp Number of components to retain. Default 5.
#'
#' @return An object of class \code{"ePCfit"}: a list with \code{scores}
#'   (n x n_comp matrix), \code{rotation}, \code{center}, \code{scale},
#'   \code{var_explained} (proportion per axis) and \code{vars}. Use
#'   \code{\link{predict.ePCfit}} to project new data.
#'
#' @examples
#' data(durumUnits)
#' ep <- computeEPC(durumUnits)
#' round(ep$var_explained, 3)     # variance explained per axis
#' head(ep$scores)
#'
#' @seealso \code{\link{predict.ePCfit}}, \code{\link{collapseUnits}}
#' @export
computeEPC <- function(data,
                       vars = c(paste0("BIO", 1:19), "Altitude"),
                       n_comp = 5) {
  vars <- intersect(vars, names(data))
  X <- as.matrix(data[, vars, drop = FALSE])
  storage.mode(X) <- "double"
  ctr <- colMeans(X)
  scl <- apply(X, 2, stats::sd)
  scl[scl == 0] <- 1
  Xs <- scale(X, center = ctr, scale = scl)
  pc <- stats::prcomp(Xs, center = FALSE, scale. = FALSE)
  keep <- seq_len(min(n_comp, ncol(pc$rotation)))
  ve <- (pc$sdev^2 / sum(pc$sdev^2))[keep]
  scores <- pc$x[, keep, drop = FALSE]
  colnames(scores) <- paste0("ePC", keep)
  structure(list(scores = scores,
                 rotation = pc$rotation[, keep, drop = FALSE],
                 center = ctr, scale = scl,
                 var_explained = ve, vars = vars),
            class = "ePCfit")
}

#' Project New Data onto a Fitted Environmental PCA
#'
#' @param object An \code{"ePCfit"} from \code{\link{computeEPC}}.
#' @param newdata A data frame containing the same predictor columns.
#' @param ... Unused.
#' @return A matrix of ePC scores for \code{newdata}.
#' @examples
#' data(durumUnits)
#' ep <- computeEPC(durumUnits[1:800, ])
#' newscores <- predict(ep, durumUnits[801:1060, ])
#' @seealso \code{\link{computeEPC}}
#' @export
#' @method predict ePCfit
predict.ePCfit <- function(object, newdata, ...) {
  X <- as.matrix(newdata[, object$vars, drop = FALSE])
  storage.mode(X) <- "double"
  Xs <- scale(X, center = object$center, scale = object$scale)
  scores <- Xs %*% object$rotation
  colnames(scores) <- colnames(object$scores)
  scores
}

#' @export
#' @method print ePCfit
print.ePCfit <- function(x, ...) {
  cat("<ePCfit>  ", ncol(x$scores), " components on ", length(x$vars),
      " predictors\n", sep = "")
  cat("  Variance explained: ",
      paste(sprintf("%.1f%%", 100 * x$var_explained), collapse = ", "),
      "\n", sep = "")
  cat("  Cumulative: ", sprintf("%.1f%%", 100 * sum(x$var_explained)), "\n", sep = "")
  invisible(x)
}
