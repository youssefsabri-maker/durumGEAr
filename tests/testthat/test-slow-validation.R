# SKIPPED by default. `test-cv-residualize.R` covers correctness (fold
# grouping, verdict logic, structure) at cheap settings in seconds; this file
# re-runs the same functions at (closer to) production settings to sanity-check
# that the permutation null still looks reasonable at scale. This is a
# validation/sanity check, not a correctness test, so it does not belong in
# the default `devtools::test()` path -- it is gated behind an environment
# variable and takes on the order of many minutes.
#
# Run explicitly with:
#   Sys.setenv(RUN_VALIDATION_TESTS = "true"); devtools::test()

skip_if_not(identical(Sys.getenv("RUN_VALIDATION_TESTS"), "true"),
            "set RUN_VALIDATION_TESTS=true to run the slow production-settings validation checks")

test_that("residualize() permutation null is well-behaved near production settings [SLOW]", {
  data(durumUnits)
  rz <- residualize(durumUnits, seeds = 42:43, k = 5, num.trees = 100, n_perm = 100)
  expect_s3_class(rz, "geaResidual")
  expect_true(all(rz$p_value > 0 & rz$p_value <= 1))
  # p-value resolution should not be sitting at the n_perm=100 floor for every
  # target if there's real heterogeneity in the underlying effect sizes
  floor_p <- 1 / (100 + 1)
  expect_false(all(abs(rz$p_value - floor_p) < 1e-9))
})

test_that("confoundGapTest() null distribution is well-behaved near production settings [SLOW]", {
  data(durumUnits)
  gt <- confoundGapTest(durumUnits, seed = 42, n_perm = 100, num.trees = 100)
  expect_s3_class(gt, "geaGapTest")
  expect_true(gt$p_value > 0 && gt$p_value <= 1)
  expect_equal(length(gt$null_gap), 100)
})
