#' Collapse Pseudoreplicated Rows to Independent Effective Units
#'
#' Collapses accession-level rows that share a collection site to independent
#' \code{(SiteCode x Cluster)} effective units, and returns a full cleaning
#' audit. This is the single most important preprocessing step in a GEA analysis
#' of genebank data: because climate predictors are derived from coordinates,
#' all accessions at one site carry an \emph{identical} predictor vector, so
#' treating them as independent observations inflates the effective sample size
#' and biases every downstream cross-validation statistic (Hurlbert, 1984).
#'
#' \strong{Intuition:} measuring the same site 100 times does not give 100
#' independent facts about climate-genetics coupling; it gives one. Collapsing
#' to unit means restores statistical independence before any modelling.
#'
#' @param data A data frame of accession rows. Must contain \code{group_cols},
#'   the \code{predictors}, the \code{targets}, and (optionally) the flag columns.
#' @param group_cols Character vector of columns defining an effective unit.
#'   Default \code{c("SiteCode", "Cluster")}.
#' @param predictors Character vector of predictor columns that are constant
#'   within a group (verified, then taken from the first row). Default: BIO1-BIO19,
#'   Altitude, Latitude, Longitude.
#' @param targets Character vector of target columns to average over non-cloned
#'   rows. Default \code{c("gPC1","gPC2","gPC3","gPC4","gPC5")}.
#' @param cloned_flag Name of the logical column flagging cloned-target rows, or
#'   \code{NULL} to skip. Default \code{"is_cloned_gpc"}.
#' @param confident_flag Name of the logical column flagging admixture-confident
#'   rows, or \code{NULL}. Default \code{"admix_confident"}.
#' @param constancy_tol Numeric tolerance for the within-group predictor
#'   constancy check. Default \code{1e-6}.
#'
#' @return A list of class \code{"unitCollapse"} with elements:
#'   \describe{
#'     \item{units}{Data frame of all collapsed effective units.}
#'     \item{confident}{Data frame of confident units (clean gPC mean AND
#'       admixture-confident) - the canonical modelling frame.}
#'     \item{audit}{Named list of cleaning counts (raw rows, exact duplicates,
#'       cloned rows, IG duplicates, min membership confidence, incomplete-admixture
#'       rows, max within-group predictor SD, unit counts).}
#'   }
#'
#' @examples
#' data(durumRaw)
#' cu <- collapseUnits(durumRaw)
#' cu$audit$n_raw          # 4545
#' nrow(cu$confident)      # 1060
#' print(cu)
#'
#' @references
#' Hurlbert, S.H. (1984). Pseudoreplication and the design of ecological field
#'   experiments. \emph{Ecological Monographs}, 54(2), 187-211.
#' @seealso \code{\link{computeEPC}}, \code{\link{spatialBlockCV}}
#' @export
collapseUnits <- function(data,
                          group_cols = c("SiteCode", "Cluster"),
                          predictors = c(paste0("BIO", 1:19), "Altitude",
                                         "Latitude", "Longitude"),
                          targets = paste0("gPC", 1:5),
                          cloned_flag = "is_cloned_gpc",
                          confident_flag = "admix_confident",
                          constancy_tol = 1e-6) {

  stopifnot(is.data.frame(data))
  need <- c(group_cols, intersect(predictors, names(data)), targets)
  missing_cols <- setdiff(c(group_cols, targets), names(data))
  if (length(missing_cols))
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  predictors <- intersect(predictors, names(data))

  audit <- list()
  audit$n_raw <- nrow(data)
  audit$n_exact_duplicates <- sum(duplicated(data))

  # cloned-target rows
  if (!is.null(cloned_flag) && cloned_flag %in% names(data)) {
    cloned <- as.logical(data[[cloned_flag]])
    cloned[is.na(cloned)] <- FALSE
  } else {
    cloned <- rep(FALSE, nrow(data))
  }
  audit$n_cloned_gpc <- sum(cloned)

  if ("is_ig_duplicate" %in% names(data)) {
    igdup <- as.logical(data[["is_ig_duplicate"]]); igdup[is.na(igdup)] <- FALSE
    audit$n_ig_duplicate <- sum(igdup)
  } else audit$n_ig_duplicate <- NA_integer_

  if ("membership_confidence" %in% names(data))
    audit$min_membership_confidence <- min(data[["membership_confidence"]], na.rm = TRUE)
  else audit$min_membership_confidence <- NA_real_

  if ("incomplete_admixture" %in% names(data)) {
    inc <- as.logical(data[["incomplete_admixture"]]); inc[is.na(inc)] <- FALSE
    audit$n_incomplete_admixture <- sum(inc)
  } else audit$n_incomplete_admixture <- NA_integer_

  gkey <- interaction(data[group_cols], drop = TRUE, sep = "::")

  # within-group predictor constancy check
  max_sd <- 0
  for (p in predictors) {
    sds <- tapply(data[[p]], gkey, function(v) stats::sd(v, na.rm = TRUE))
    max_sd <- max(max_sd, max(sds, na.rm = TRUE), na.rm = TRUE)
  }
  audit$max_within_group_predictor_sd <- max_sd
  audit$constancy_pass <- max_sd < constancy_tol

  # collapse
  idx_by_group <- split(seq_len(nrow(data)), gkey)
  rows <- lapply(names(idx_by_group), function(g) {
    ix <- idx_by_group[[g]]
    first <- ix[1]
    out <- data[first, c(group_cols, predictors), drop = FALSE]
    non_cloned <- ix[!cloned[ix]]
    n_used <- length(non_cloned)
    if (n_used == 0) {
      tvals <- rep(NA_real_, length(targets)); miss_all <- 1L
    } else {
      tvals <- sapply(targets, function(t) mean(data[[t]][non_cloned], na.rm = TRUE))
      miss_all <- 0L
    }
    names(tvals) <- targets
    out$n_group <- length(ix)
    out$n_used_for_gpc <- n_used
    out$gpc_missing_all_cloned <- miss_all
    if (!is.null(confident_flag) && confident_flag %in% names(data))
      out$admix_confident <- as.integer(any(as.logical(data[[confident_flag]][ix]) %in% TRUE))
    else out$admix_confident <- 1L
    if ("Country" %in% names(data)) out$Country <- data[["Country"]][first]
    cbind(out, as.data.frame(as.list(tvals)))
  })
  units <- do.call(rbind, rows)
  rownames(units) <- NULL

  audit$n_units <- nrow(units)
  audit$n_gpc_missing_all_cloned <- sum(units$gpc_missing_all_cloned == 1L)

  confident <- units[units$gpc_missing_all_cloned == 0L & units$admix_confident == 1L, ]
  rownames(confident) <- NULL
  audit$n_confident <- nrow(confident)
  audit$n_countries_confident <- if ("Country" %in% names(confident))
    length(unique(confident$Country)) else NA_integer_

  structure(list(units = units, confident = confident, audit = audit),
            class = "unitCollapse")
}

#' @export
#' @method print unitCollapse
print.unitCollapse <- function(x, ...) {
  a <- x$audit
  cat("<unitCollapse>\n")
  cat(sprintf("  Raw rows ................... %d\n", a$n_raw))
  cat(sprintf("  Exact duplicate rows ....... %d\n", a$n_exact_duplicates))
  cat(sprintf("  Cloned-gPC rows ............ %d\n", a$n_cloned_gpc))
  cat(sprintf("  IG-duplicate rows .......... %s\n", a$n_ig_duplicate))
  cat(sprintf("  Min membership confidence .. %.4f\n", a$min_membership_confidence))
  cat(sprintf("  Incomplete-admixture rows .. %s\n", a$n_incomplete_admixture))
  cat(sprintf("  Max within-group pred. SD .. %.2e  (%s)\n",
              a$max_within_group_predictor_sd,
              if (isTRUE(a$constancy_pass)) "PASS - predictors constant" else "FAIL"))
  cat(sprintf("  Effective units ............ %d\n", a$n_units))
  cat(sprintf("  gpc_missing_all_cloned ..... %d (excluded)\n", a$n_gpc_missing_all_cloned))
  cat(sprintf("  Confident modelling units .. %d  (%s countries)\n",
              a$n_confident, a$n_countries_confident))
  invisible(x)
}
