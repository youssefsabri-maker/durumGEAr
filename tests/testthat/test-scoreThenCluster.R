test_that("scoreThenCluster runs leakage-free and returns a coherent geaStage2", {
  skip_if_not_installed("ranger")
  data(durumUnits)
  # small, fast configuration
  s2 <- scoreThenCluster(durumUnits, per_tree = FALSE, num.trees = 30)

  expect_s3_class(s2, "geaStage2")
  # accuracy and ceiling are valid probabilities
  expect_true(s2$accuracy >= 0 && s2$accuracy <= 1)
  expect_true(s2$ceiling  >= 0 && s2$ceiling  <= 1)
  # the ceiling (perfect scores) must not be below the honest OOF accuracy
  expect_gte(s2$ceiling + 1e-8, s2$accuracy)
  # confusion matrix is square over the clusters, and totals to n
  expect_equal(sum(s2$confusion), nrow(durumUnits))
  expect_equal(nrow(s2$confusion), ncol(s2$confusion))
  # posteriors are a proper distribution per row
  expect_equal(unname(rowSums(s2$posterior)), rep(1, nrow(durumUnits)), tolerance = 1e-6)
  # entropy calibration has the expected columns
  expect_true(all(c("quartile", "accuracy", "n") %in% names(s2$entropy_calibration)))
})

test_that("cv_ceiling = TRUE (issue #7) gives an out-of-fold ceiling that is <= the in-sample one", {
  skip_if_not_installed("ranger")
  data(durumUnits)
  s2_insample <- scoreThenCluster(durumUnits, per_tree = FALSE, num.trees = 30,
                                  cv_ceiling = FALSE)
  s2_oof <- scoreThenCluster(durumUnits, per_tree = FALSE, num.trees = 30,
                             cv_ceiling = TRUE)
  expect_equal(s2_insample$ceiling_type,
               "in-sample (optimistic upper bound - see Caveat in ?scoreThenCluster)")
  expect_equal(s2_oof$ceiling_type, "out-of-fold (spatial-block, leakage-free)")
  # the in-sample ceiling is an optimistic upper bound: it should not be
  # exceeded by the honest out-of-fold ceiling on the same data/folds
  expect_gte(s2_insample$ceiling + 1e-8, s2_oof$ceiling)
  # both ceilings must still dominate their own honest accuracy
  expect_gte(s2_oof$ceiling + 1e-8, s2_oof$accuracy)
})

test_that(".qdaFit/.qdaPosterior recover an easy 2-class problem", {
  set.seed(1)
  X <- rbind(matrix(rnorm(200, -3), 100, 2), matrix(rnorm(200, 3), 100, 2))
  y <- rep(c("A", "B"), each = 100)
  fit <- durumGEAr:::.qdaFit(X, y)
  post <- durumGEAr:::.qdaPosterior(fit, X)
  pred <- fit$classes[max.col(post)]
  expect_gt(mean(pred == y), 0.95)
})
