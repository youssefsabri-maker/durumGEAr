#' Per-Target and Per-Target-Mean Regression Metrics
#'
#' Computes R2 and RMSE for a multi-output regression, reporting each target
#' separately and returning the \strong{per-target-mean R2} as the canonical
#' summary. This function deliberately does \emph{not} return a pooled R2 over
#' the flattened target matrix, because a pooled score is dominated by the
#' highest-variance target (here gPC1) and systematically overstates skill.
#'
#' \strong{Intuition:} averaging the five per-axis R2 values weights each genetic
#' axis equally; pooling first concatenates all targets and lets the easiest,
#' most variable axis dominate the single number.
#'
#' @param observed A numeric matrix or data frame of observed target values
#'   (rows = units, columns = targets).
#' @param predicted A numeric matrix or data frame of predictions, same shape.
#'
#' @return A list of class \code{"geaMetrics"} with \code{per_target} (a data
#'   frame with R2 and RMSE per target), \code{mean_R2} (per-target-mean R2),
#'   and \code{mean_RMSE}.
#'
#' @examples
#' set.seed(1)
#' obs <- matrix(rnorm(50), 10, 5)
#' pred <- obs + matrix(rnorm(50, sd = 0.3), 10, 5)
#' getMetrics(obs, pred)
#'
#' @seealso \code{\link{spatialBlockCV}}, \code{\link{locoCV}}
#' @export
getMetrics <- function(observed, predicted) {
  observed <- as.matrix(observed); predicted <- as.matrix(predicted)
  stopifnot(all(dim(observed) == dim(predicted)))
  nt <- ncol(observed)
  tnames <- colnames(observed)
  if (is.null(tnames)) tnames <- paste0("target", seq_len(nt))
  r2 <- numeric(nt); rmse <- numeric(nt)
  for (j in seq_len(nt)) {
    o <- observed[, j]; p <- predicted[, j]
    ok <- is.finite(o) & is.finite(p)
    o <- o[ok]; p <- p[ok]
    ss_res <- sum((o - p)^2)
    ss_tot <- sum((o - mean(o))^2)
    r2[j] <- if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_
    rmse[j] <- sqrt(mean((o - p)^2))
  }
  structure(list(
    per_target = data.frame(target = tnames, R2 = r2, RMSE = rmse,
                            stringsAsFactors = FALSE),
    mean_R2 = mean(r2, na.rm = TRUE),
    mean_RMSE = mean(rmse, na.rm = TRUE)),
    class = "geaMetrics")
}

#' @export
#' @method print geaMetrics
print.geaMetrics <- function(x, ...) {
  cat("<geaMetrics>\n")
  pt <- x$per_target
  for (i in seq_len(nrow(pt)))
    cat(sprintf("  %-8s R2 = %+.4f   RMSE = %.4f\n", pt$target[i], pt$R2[i], pt$RMSE[i]))
  cat(sprintf("  %-8s R2 = %+.4f   (canonical)\n", "MEAN", x$mean_R2))
  invisible(x)
}
