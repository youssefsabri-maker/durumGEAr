#' Raw durum wheat accession export (flagged)
#'
#' The accession-level source export for the durum wheat GEA project, with
#' quality-control flags added. Each row is one genotyped accession. Climate
#' predictors (BIO1-BIO19, Altitude) are derived from the collection-site
#' coordinates and are therefore \emph{identical} for accessions sharing a site
#' - the source of the pseudoreplication addressed by \code{\link{collapseUnits}}.
#'
#' @format A data frame with 4,545 rows and 43 columns:
#' \describe{
#'   \item{individual}{Accession identifier.}
#'   \item{SiteCode}{Collection-site code (many accessions per site).}
#'   \item{ICARDA_IG}{Genebank accession (IG) number.}
#'   \item{Country}{Country of origin (43 levels).}
#'   \item{Latitude, Longitude, Altitude}{Georeference and elevation (m).}
#'   \item{gPC1, gPC2, gPC3, gPC4, gPC5}{Genetic principal-component scores (targets).}
#'   \item{Cluster}{K7-admixture-derived genetic cluster (5 levels). Never used
#'     as a predictor - it is a genetic label, not an environmental one.}
#'   \item{K7_Q1, K7_Q2, K7_Q3, K7_Q4, K7_Q5}{K=7 admixture proportions (collapsed to 5 retained).}
#'   \item{BIO1, BIO2, BIO3, BIO4, BIO5, BIO6, BIO7, BIO8, BIO9, BIO10, BIO11, BIO12, BIO13, BIO14, BIO15, BIO16, BIO17, BIO18, BIO19}{The 19 standard WorldClim bioclimatic variables.}
#'   \item{is_cloned_gpc}{TRUE if this row's gPC vector is a clone within its group (58 rows).}
#'   \item{is_ig_duplicate}{TRUE if the accession shares an IG number (2 rows).}
#'   \item{membership_confidence}{max(K7_Q1..Q5); minimum observed is 0.500.}
#'   \item{q_sum}{Sum of retained admixture proportions.}
#'   \item{incomplete_admixture}{TRUE if q_sum < 0.95 (957 rows; retained).}
#'   \item{admix_confident}{TRUE if the unit passes the admixture-confidence rule.}
#' }
#' @source Durum wheat genebank collection, integrated climate/altitude export.
"durumRaw"

#' Confident durum wheat modelling units
#'
#' The canonical analysis frame: pseudoreplicated accession rows collapsed to
#' independent (SiteCode x Cluster) effective units, restricted to units with a
#' clean genetic mean and confident admixture assignment. This is the output of
#' \code{\link{collapseUnits}} on \code{\link{durumRaw}} and is the frame on
#' which every modelling result in the package vignette is computed.
#'
#' @format A data frame with 1,060 rows and 30 columns:
#' \describe{
#'   \item{SiteCode}{Collection-site code.}
#'   \item{Cluster}{Genetic cluster (Cluster_1=104, _2=98, _3=83, _4=167, _5=608).}
#'   \item{Country}{Country of origin (43 levels).}
#'   \item{Latitude, Longitude, Altitude}{Georeference and elevation.}
#'   \item{BIO1, BIO2, BIO3, BIO4, BIO5, BIO6, BIO7, BIO8, BIO9, BIO10, BIO11, BIO12, BIO13, BIO14, BIO15, BIO16, BIO17, BIO18, BIO19}{Bioclimatic predictors (constant within a site x cluster group).}
#'   \item{gPC1, gPC2, gPC3, gPC4, gPC5}{Group-mean genetic PC scores (targets).}
#' }
#' @source Derived from \code{\link{durumRaw}} via \code{\link{collapseUnits}}.
"durumUnits"

#' Cached residualization validation result
#'
#' A worked \code{\link{residualize}} result on \code{\link{durumUnits}}, shipped
#' so the deployment reliability table and vignette can be shown without
#' recomputing. Produced at deliberately reduced, laptop-runnable settings
#' (\code{seeds = 42:44}, \code{n_perm = 30}, \code{num.trees = 150}, \code{k = 5});
#' rerun \code{residualize()} at production settings (\code{n_perm = 100},
#' \code{num.trees = 600}) for a publication-grade significance verdict. The gPC2
#' residual-scale R2 (approx 0.18) and gPC5 (approx 0.06) are the signal-carrying
#' axes; gPC1/gPC4 are flagged \code{"fragile"} (country-identity artifacts).
#'
#' @format A \code{geaResidual} object (a list) with fields including
#'   \code{per_target_R2}, \code{per_target_R2_CI}, \code{per_seed_R2},
#'   \code{obs_R2_test}, \code{p_value}, \code{p_value_CI}, \code{verdict},
#'   \code{verdict_fragility}, \code{null_mean}, \code{null_sd}, and
#'   \code{null_mat}. See \code{\link{residualize}} for field definitions.
#' @seealso \code{\link{residualize}}, \code{\link{fitGeneticScoreModel}}
#' @source \code{residualize(durumUnits, seeds = 42:44, n_perm = 30, num.trees = 150, k = 5)}
"durum_residualize_results"

#' Cached leave-one-country-out validation result
#'
#' A worked \code{\link{locoCV}} result on \code{\link{durumUnits}}, shipped so the
#' country-level transferability breakdown can be shown without recomputing.
#' Produced at reduced settings (\code{num.trees = 150}); rerun \code{locoCV()} at
#' production settings for a final figure. Only \code{n_beat_baseline}
#' (12 of 34 held-out countries) improve on their own mean-prediction baseline -
#' the quantitative basis for the package's headline finding that transferable
#' climate-to-genetics skill is near-zero out of country.
#'
#' @format A \code{geaCV} object (a list) with fields \code{oof}, \code{metrics},
#'   \code{per_group} (per-country n, R2, and \code{beats_baseline}),
#'   \code{per_group_target} (per-country x per-target R2), \code{n_beat_baseline},
#'   \code{n_groups}, and \code{median_group_R2}. See \code{\link{locoCV}}.
#' @seealso \code{\link{locoCV}}, \code{\link{checkExtrapolationRisk}}
#' @source \code{locoCV(durumUnits, num.trees = 150)}
"durum_loco_results"
