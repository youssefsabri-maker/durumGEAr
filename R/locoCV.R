#' Leave-One-Country-Out Cross-Validation (Extrapolation Skill)
#'
#' Evaluates multi-output regression skill under leave-one-group-out (LOGO)
#' cross-validation, holding out one whole country (or other coarse spatial
#' grouping) at a time. This is the harsh, honest test of \emph{extrapolation}
#' to unseen geography, as opposed to the interpolation estimate from
#' \code{\link{spatialBlockCV}}.
#'
#' \strong{Intuition:} interpolation asks "can you predict a new site in a
#' country you already know?"; LOCO asks the harder "can you predict an entirely
#' new country you have never seen?". A model that memorizes country identity
#' scores well on the former and collapses on the latter.
#'
#' A large gap between \code{\link{spatialBlockCV}} and \code{locoCV} is the
#' signature of geographic-identity confounding (see the package vignette and
#' \code{\link{residualize}}).
#'
#' @param data A confident modelling frame (e.g. \code{durumUnits}).
#' @param predictors,targets Predictor and target columns (see \code{\link{spatialBlockCV}}).
#' @param group Coarse group to hold out. Default \code{"Country"}.
#' @param seed Random seed for the ExtraTrees fitter. Default 42.
#' @param floor Lower clamp for a country's R2 (small held-out countries can
#'   produce extreme negative values). Default \code{-1}.
#' @param num.trees,min.node.size Passed to the internal ExtraTrees fitter.
#'
#' @return An object of class \code{"geaCV"} with \code{oof}, \code{metrics}
#'   (per-target-mean R2 across all held-out units), \code{per_group}
#'   (a data frame of per-country mean R2 and whether it beat the global-mean
#'   baseline), \code{per_group_target} (a countries x targets matrix of the
#'   individual per-target R2 behind that mean, for country-level breakdowns),
#'   and \code{n_beat_baseline}.
#'
#' @examples
#' \donttest{
#' data(durumUnits)
#' lo <- locoCV(durumUnits)
#' lo$metrics$mean_R2
#' lo$n_beat_baseline           # how many countries beat the no-skill baseline
#' }
#' @seealso \code{\link{spatialBlockCV}}, \code{\link{residualize}}
#' @export
locoCV <- function(data,
                   predictors = c(paste0("BIO", 1:19), "Altitude",
                                  "Latitude", "Longitude"),
                   targets = paste0("gPC", 1:5),
                   group = "Country", seed = 42, floor = -1,
                   num.trees = 600, min.node.size = 3) {
  predictors <- intersect(predictors, names(data))
  Y <- as.matrix(data[, targets, drop = FALSE])
  grp <- as.character(data[[group]])
  ug <- unique(grp)
  oof <- matrix(NA_real_, nrow(data), length(targets), dimnames = list(NULL, targets))
  per_group <- data.frame(group = ug, n = NA_integer_, R2 = NA_real_,
                          beats_baseline = NA, stringsAsFactors = FALSE)
  # per-country x per-target R2 breakdown (issue #16 / spec §7): per_group$R2 is
  # the per-target MEAN per country; this matrix keeps individual targets so a
  # caller can see, e.g., a country extrapolating for gPC1 but not gPC3.
  per_group_target <- matrix(NA_real_, length(ug), length(targets),
                             dimnames = list(ug, targets))
  for (i in seq_along(ug)) {
    g <- ug[i]
    te <- which(grp == g); tr <- which(grp != g)
    if (length(tr) < 5 || length(te) < 1) next
    pred <- .fitPredictMO(data[tr, ], data[te, ], predictors, targets,
                          num.trees = num.trees, min.node.size = min.node.size,
                          seed = seed)
    oof[te, ] <- pred
    # per-country per-target-mean R2 vs the global-train-mean baseline
    r2s <- numeric(length(targets)); base_r2s <- numeric(length(targets))
    for (j in seq_along(targets)) {
      o <- Y[te, j]; p <- pred[, j]; gm <- mean(Y[tr, j])
      ss_tot <- sum((o - mean(o))^2)
      r2s[j] <- if (ss_tot > 0) max(floor, 1 - sum((o - p)^2) / ss_tot) else NA_real_
      # baseline = predict global train mean
      base_r2s[j] <- if (ss_tot > 0) max(floor, 1 - sum((o - gm)^2) / ss_tot) else NA_real_
    }
    per_group_target[i, ] <- r2s
    per_group$n[i] <- length(te)
    per_group$R2[i] <- mean(r2s, na.rm = TRUE)
    per_group$beats_baseline[i] <- mean(r2s, na.rm = TRUE) > mean(base_r2s, na.rm = TRUE)
  }
  metrics <- getMetrics(Y, oof)
  keep <- !is.na(per_group$R2)
  per_group <- per_group[keep, ]
  per_group_target <- per_group_target[keep, , drop = FALSE]
  structure(list(oof = oof, metrics = metrics, per_group = per_group,
                 per_group_target = per_group_target,
                 n_beat_baseline = sum(per_group$beats_baseline, na.rm = TRUE),
                 n_groups = nrow(per_group),
                 median_group_R2 = stats::median(per_group$R2, na.rm = TRUE),
                 type = sprintf("leave-one-%s-out - extrapolation", group)),
            class = "geaCV")
}
