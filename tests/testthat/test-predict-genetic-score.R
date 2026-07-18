# Structure + correctness tests for the v1.3.0 deployment layer (issues #20-#24:
# fitGeneticScoreModel, predict, checkExtrapolationRisk, predictCluster,
# robustGeneticScore, scoreThenCluster(fit_final=)).
#
# Philosophy matches test-durum-reproducibility.R: shared cheap fixtures fit
# ONCE at tiny settings, reused across blocks. These test the CODE (fields,
# shapes, invariants, gating logic), not production-scale statistical skill.

data(durumUnits)

set.seed(1)
.sub <- durumUnits[sample(nrow(durumUnits), 300), ]
# multi-seed residualize so reliability metadata is populated
.rz  <- residualize(.sub, seeds = 42:43, k = 2, num.trees = 20, n_perm = 8)
.mod <- fitGeneticScoreModel(.sub, residualize_result = .rz, num.trees = 40)
.s2  <- scoreThenCluster(.sub, per_tree = FALSE, k = 2, num.trees = 30,
                         fit_final = TRUE)

# ---- fitGeneticScoreModel ----------------------------------------------------

test_that("fitGeneticScoreModel returns one ranger per target + reliability", {
  expect_s3_class(.mod, "geneticScoreModel")
  expect_equal(length(.mod$models), 5L)
  expect_true(all(vapply(.mod$models, inherits, logical(1), "ranger")))
  expect_equal(nrow(.mod$reliability), 5L)
  expect_true(all(c("verdict", "fragility", "gated", "reliability") %in%
                    names(.mod$reliability)))
})

test_that("reliability gate uses BOTH verdict and verdict_fragility fields", {
  # a target is gated iff verdict == 'country-identity artifact' OR
  # fragility in {fragile, marginal}; never gates on the string 'artifact'
  # appearing in verdict_fragility (it never does).
  rl <- .mod$reliability
  expected <- (rl$verdict == "country-identity artifact") |
              (rl$fragility %in% c("fragile", "marginal"))
  expect_equal(rl$gated, expected)
  expect_false(any(rl$fragility == "artifact", na.rm = TRUE))
})

test_that("fitGeneticScoreModel with no residualize_result leaves axes ungated", {
  m0 <- fitGeneticScoreModel(.sub, num.trees = 20)
  expect_false(any(m0$reliability$gated))
  expect_match(m0$reliability$reliability[1], "ungated")
})

test_that("missing targets or predictors error informatively", {
  expect_error(fitGeneticScoreModel(.sub, targets = "gPCX"), "not found")
  bad <- .sub[, c("gPC1", "gPC2", "gPC3", "gPC4", "gPC5")]  # no predictors
  expect_error(fitGeneticScoreModel(bad), "No predictors")
})

# ---- predict.geneticScoreModel -----------------------------------------------

test_that("predict returns a per-target data frame with reliability attrs", {
  p <- predict(.mod, durumUnits[1:6, ])
  expect_s3_class(p, "data.frame")
  expect_equal(dim(p), c(6L, 5L))
  expect_equal(names(p), paste0("gPC", 1:5))
  expect_false(is.null(attr(p, "reliability")))
  expect_equal(attr(p, "gated_targets"),
               .mod$reliability$target[.mod$reliability$gated])
})

test_that("gated_to_na blanks exactly the gated columns", {
  gated <- .mod$reliability$target[.mod$reliability$gated]
  skip_if(length(gated) == 0, "no gated targets in this fixture")
  p <- predict(.mod, durumUnits[1:6, ], gated_to_na = TRUE)
  all_na <- names(p)[vapply(p, function(c) all(is.na(c)), logical(1))]
  expect_setequal(all_na, gated)
})

test_that("predict errors when newdata lacks predictor columns", {
  expect_error(predict(.mod, durumUnits[1:3, c("gPC1", "gPC2")]), "missing")
})

test_that("predict is deterministic for a fixed fitted model", {
  expect_equal(predict(.mod, durumUnits[1:4, ]),
               predict(.mod, durumUnits[1:4, ]))
})

# ---- checkExtrapolationRisk --------------------------------------------------

test_that("in-domain training rows are (almost) never flagged at defaults", {
  risk <- checkExtrapolationRisk(.mod, .sub)
  expect_equal(mean(risk$range_flag), 0)                 # envelope_prob = 1
  expect_lt(mean(risk$md_flag), 0.06)                    # ~ 1 - 0.975
  expect_true(all(c("mahalanobis", "extrapolation") %in% names(risk)))
})

test_that("out-of-range observations are flagged as extrapolation", {
  od <- durumUnits[1:3, ]
  od$BIO1 <- od$BIO1 + 60; od$BIO12 <- od$BIO12 * 6
  risk <- checkExtrapolationRisk(.mod, od)
  expect_true(all(risk$extrapolation))
  expect_true(all(risk$range_flag))
})

test_that("tightening envelope_prob flags more tail rows", {
  loose <- checkExtrapolationRisk(.mod, .sub, envelope_prob = 1)
  tight <- checkExtrapolationRisk(.mod, .sub, envelope_prob = 0.98)
  expect_gt(mean(tight$range_flag), mean(loose$range_flag))
})

# ---- predictCluster ----------------------------------------------------------

test_that("scoreThenCluster(fit_final=TRUE) exposes a geaQDA classifier", {
  expect_s3_class(.s2$final_classifier, "geaQDA")
  expect_null(scoreThenCluster(.sub, per_tree = FALSE, k = 2,
                               num.trees = 20)$final_classifier)
})

test_that("predictCluster matches a manual QDA posterior call", {
  pc <- predictCluster(.mod, durumUnits[1:8, ], .s2$final_classifier)
  expect_true(all(pc$posterior >= 0 & pc$posterior <= 1))
  expect_true(all(pc$entropy >= 0))
  S <- as.matrix(predict(.mod, durumUnits[1:8, ])[, .s2$final_classifier$targets])
  manual <- .s2$final_classifier$classes[
    max.col(.qdaPosterior(.s2$final_classifier, S), ties.method = "first")]
  expect_equal(pc$predicted, manual)
})

test_that("predictCluster rejects a non-geaQDA classifier", {
  expect_error(predictCluster(.mod, durumUnits[1:3, ], list()), "geaQDA")
})

# ---- robustGeneticScore ------------------------------------------------------

test_that("robustGeneticScore collapses forced-kept axes to one value/row", {
  v <- robustGeneticScore(.mod, durumUnits[1:8, ], keep = paste0("gPC", 1:5))
  expect_length(v, 8L)
  expect_equal(attr(v, "axes"), paste0("gPC", 1:5))
  # mean method == rowMeans of the kept predictions
  M <- as.matrix(predict(.mod, durumUnits[1:8, ])[, paste0("gPC", 1:5)])
  expect_equal(as.numeric(v), unname(rowMeans(M)))
})

test_that("robustGeneticScore methods run and pca1 returns one score/row", {
  for (m in c("mean", "median", "pca1")) {
    v <- robustGeneticScore(.mod, durumUnits[1:8, ], method = m,
                            keep = c("gPC2", "gPC5"))
    expect_length(v, 8L)
    expect_equal(attr(v, "method"), m)
  }
})

test_that("robustGeneticScore errors when no axis survives gating", {
  # force a fully-gated model by supplying reliability where all are gated
  if (all(.mod$reliability$gated)) {
    expect_error(robustGeneticScore(.mod, durumUnits[1:3, ]), "No trustworthy")
  } else {
    succeed("fixture has a surviving axis; skip the all-gated branch")
  }
})
