# Tests the max-statistic permutation null (issue #4) and the matched
# observed-vs-null hyperparameters (issue #5).

test_that("driverAnalysis's null is a max over ALL candidate predictors, not just the top one", {
  data(durumUnits)
  dr <- driverAnalysis(durumUnits, "gPC2", seeds = 42, n_perm = 8,
                       num.trees = 20, min.node.size = 5)
  expect_s3_class(dr, "geaDrivers")
  expect_equal(length(dr$null_max), 8)
  expect_equal(dr$n_predictors_tested, nrow(dr$univariate))
  expect_true(dr$n_predictors_tested > 1)
  expect_true(dr$p_value > 0 && dr$p_value <= 1)
  # observed_R2 must equal the top row of the univariate table
  expect_equal(dr$observed_R2, dr$univariate$R2[1], tolerance = 1e-9)
})

test_that("driverAnalysis's num.trees is threaded into BOTH observed and null fits (issue #5 fix)", {
  data(durumUnits)
  # num.trees is a hyperparameter of the ranger fit used for BOTH the observed
  # top-predictor statistic and every null permutation (issue #5: previously
  # the observed fit used num.trees=300 while the null loop was hardcoded to
  # num.trees=150, biasing the null noisier than the signal it was compared
  # against). Under the pre-fix bug, changing the num.trees ARGUMENT would
  # have moved observed_R2 (which reads the argument) but NOT null_max (which
  # ignored it and always used the hardcoded value) -- so the regression test
  # is that null_max actually moves when num.trees changes, not that
  # observed_R2 stays fixed (it legitimately depends on num.trees too, since
  # it comes from the same kind of ranger fit).
  dr_lo <- driverAnalysis(durumUnits, "gPC2", seeds = 42, n_perm = 5, num.trees = 20)
  dr_hi <- driverAnalysis(durumUnits, "gPC2", seeds = 42, n_perm = 5, num.trees = 80)
  # Both the observed statistic and the null should respond to num.trees --
  # if either were hardcoded/decoupled from the argument (the issue #5 bug),
  # one of these would stay fixed while the other moved.
  expect_false(isTRUE(all.equal(dr_lo$observed_R2, dr_hi$observed_R2, tolerance = 1e-9)))
  expect_false(isTRUE(all.equal(dr_lo$null_max, dr_hi$null_max, tolerance = 0.01)))
})
