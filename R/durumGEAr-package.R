#' durumGEAr: Confound-Aware Genotype-Environment Association Modelling
#'
#' The \pkg{durumGEAr} package provides a defensible, reproducible workflow for
#' genotype-environment association (GEA) modelling of genebank accession data.
#' It was developed on a durum wheat collection in which the predictors are
#' bioclimatic variables (BIO1-BIO19) plus altitude and coordinates, and the
#' targets are five genetic principal-component scores (gPC1-gPC5).
#'
#' @section The central methodological problem:
#' Genebank accessions collected at the same site share \emph{identical}
#' coordinate-derived climate predictors and near-duplicate genetic targets.
#' Treating each accession row as an independent observation inflates the
#' effective sample size (pseudoreplication) and, more subtly, lets a model
#' achieve a high cross-validated score by memorizing \emph{country identity}
#' rather than learning a transferable climate-to-genetics relationship. This
#' package supplies the tools to (a) remove the pseudoreplication, (b) measure
#' predictive skill honestly, and (c) separate genuine within-country climate
#' signal from geographic-identity confounding.
#'
#' @section Core functions:
#' \describe{
#'   \item{\code{\link{collapseUnits}}}{Collapse pseudoreplicated site x cluster
#'     rows to independent effective units, with a full cleaning audit.}
#'   \item{\code{\link{computeEPC}}}{Fit an environmental PCA on a confident
#'     modelling frame (leakage-safe: scaler and rotation from training data).}
#'   \item{\code{\link{spatialBlockCV}}}{Spatial-block cross-validation grouped
#'     by whole site (interpolation skill).}
#'   \item{\code{\link{locoCV}}}{Leave-one-country-out cross-validation
#'     (extrapolation skill).}
#'   \item{\code{\link{scoreThenCluster}}}{Leakage-free Stage-2 assignment of
#'     accessions to discrete genetic clusters from predicted genetic scores,
#'     with a perfect-score ceiling, entropy calibration, and off-manifold
#'     typicality flags.}
#'   \item{\code{\link{residualize}}}{Fold-safe country-mean residualization to
#'     strip between-country variance without leakage.}
#'   \item{\code{\link{driverAnalysis}}}{Univariate and permutation-importance
#'     driver ranking on a residual target, with a within-country permutation
#'     significance test.}
#'   \item{\code{\link{getMetrics}}}{Per-target and per-target-mean R2 / RMSE.}
#'   \item{\code{\link{mapAccessions}}}{Plot accession or unit locations.}
#' }
#'
#' @section Canonical evaluation rule:
#' Throughout this package, multi-output regression skill is reported as the
#' \strong{per-target-mean R2} - the mean of the five individual per-gPC R2
#' values - never the pooled R2 over the flattened target matrix (which is
#' dominated by the highest-variance target and overstates skill).
#'
#' @references
#' Hurlbert, S.H. (1984). Pseudoreplication and the design of ecological field
#'   experiments. \emph{Ecological Monographs}, 54(2), 187-211.
#'
#' Roberts, D.R. et al. (2017). Cross-validation strategies for data with
#'   temporal, spatial, hierarchical, or phylogenetic structure.
#'   \emph{Ecography}, 40(8), 913-929.
#'
#' Geurts, P., Ernst, D. & Wehenkel, L. (2006). Extremely randomized trees.
#'   \emph{Machine Learning}, 63(1), 3-42.
#'
#' @keywords internal
"_PACKAGE"
