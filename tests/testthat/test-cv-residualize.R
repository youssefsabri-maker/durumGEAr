# Design note (post-review): unit tests here check that the CV/residualization
# *code* is correct (fold grouping, structure, verdict logic), not that the
# underlying permutation null is statistically well-calibrated -- that is a
# validation-study question, answered once in the vignette/thesis analysis at
# production settings (n_perm=100, num.trees=600), not on every `devtools::test()`.
# Accordingly:
#   - `.residualVerdict()` (the actual logic issue #8 broke) is tested directly
#     on synthetic p_value/obs_R2_test vectors -- zero ranger fits, instant.
#   - CV/residualize structural tests use k = 2 folds and tiny num.trees (fold
#     count and tree count only change CV granularity/model capacity, not the
#     code path being tested) and a SHARED fixture fit once and reused, instead
#     of each test_that() re-fitting from scratch.
#   - Exactly ONE small end-to-end permutation smoke test (tiny n_perm, tiny
#     num.trees, k = 2) proves the permutation *plumbing* isn't broken, without
#     trying to be a real calibration check.
# This keeps `devtools::test()` fast (seconds, not tens of minutes) while still
# exercising every code path a regression could break.

# ---- Fast, zero-fit test of the verdict logic itself (issue #8) -----------

test_that(".residualVerdict requires BOTH p < alpha AND obs_R2_test > 0 (issue #8)", {
  p <- c(0.01, 0.01, 0.5, 0.5, 0.01)
  r2 <- c(0.10, -0.10, 0.10, -0.10, 0)
  alpha <- 0.05
  v <- .residualVerdict(p, r2, alpha)
  expect_equal(v, c("real within-country signal",   # sig p,  positive R2
                    "country-identity artifact",     # sig p,  negative R2 <- issue #8
                    "country-identity artifact",      # ns  p,  positive R2
                    "country-identity artifact",      # ns  p,  negative R2
                    "country-identity artifact"))     # sig p,  zero R2 (not > 0)
})

test_that(".residualVerdict never assigns 'real signal' to negative R2, at any p-value (issue #8)", {
  # Regression test: sweep p-values from far-significant to far-nonsignificant
  # with a fixed negative R2 -- the verdict must be "country-identity artifact"
  # in every case, regardless of how small p gets.
  p_grid <- c(1e-6, 1e-3, 0.01, 0.049, 0.05, 0.5, 0.99)
  v <- .residualVerdict(p_grid, rep(-0.05, length(p_grid)), alpha = 0.05)
  expect_true(all(v == "country-identity artifact"))
})

# ---- Shared cheap fixtures (k = 2, tiny trees, fit once) -------------------

data(durumUnits)
.rz_fixture <- residualize(durumUnits, seeds = 42, k = 2, num.trees = 20, n_perm = 5)
.sb_fixture <- spatialBlockCV(durumUnits, seeds = 42, k = 2, num.trees = 20)
.lo_fixture <- locoCV(durumUnits, seed = 42, num.trees = 20)
.nv_fixture <- naiveCV(durumUnits, seeds = 42, k = 2, num.trees = 20)
.gt_fixture <- confoundGapTest(durumUnits, seed = 42, n_perm = 3,
                               num.trees = 15, min.node.size = 5)

test_that("spatialBlockCV keeps whole groups together and returns geaCV", {
  cv <- .sb_fixture
  expect_s3_class(cv, "geaCV")
  expect_equal(nrow(cv$metrics$per_target), 5)
  expect_true(is.finite(cv$metrics$mean_R2))
})

test_that("locoCV returns per-group table and is harsher than interpolation", {
  expect_s3_class(.lo_fixture, "geaCV")
  expect_true(all(c("group", "R2", "beats_baseline") %in% names(.lo_fixture$per_group)))
  # the honest extrapolation estimate must be below interpolation (the confound gap)
  expect_lt(.lo_fixture$metrics$mean_R2, .sb_fixture$metrics$mean_R2)
})

test_that("residualize returns well-formed per-target structure", {
  rz <- .rz_fixture
  expect_s3_class(rz, "geaResidual")
  expect_equal(length(rz$per_target_R2), 5)
  expect_true(all(rz$verdict %in%
                  c("real within-country signal", "country-identity artifact")))
})

test_that("residualize's permutation test (issue #6) returns a coherent p-value per target", {
  rz <- .rz_fixture
  expect_equal(length(rz$p_value), 5)
  expect_true(all(rz$p_value > 0 & rz$p_value <= 1))
  expect_equal(names(rz$p_value), paste0("gPC", 1:5))
  # the object's own verdict must match what .residualVerdict() would compute
  # from its own p_value/obs_R2_test -- i.e. residualize() must actually be
  # calling the shared verdict logic, not some inline duplicate of it.
  expect_equal(unname(rz$verdict),
               unname(.residualVerdict(rz$p_value, rz$obs_R2_test, rz$alpha)))
})

test_that("residualize never assigns 'real signal' to negative R2 on real output (issue #8)", {
  # End-to-end confirmation (on top of the synthetic sweep above) that this
  # invariant holds on actual residualize() output, not just the isolated fn.
  rz <- .rz_fixture
  for (t in names(rz$obs_R2_test)) {
    if (rz$obs_R2_test[t] < 0) {
      expect_equal(unname(rz$verdict[t]), "country-identity artifact",
                   info = sprintf("Target %s has obs_R2_test=%+.4f (negative) but verdict=%s",
                                  t, rz$obs_R2_test[t], rz$verdict[t]))
    }
  }
})

test_that("naiveCV (issue #2, leaky ungrouped baseline) is not grouped by site", {
  nv <- .nv_fixture
  expect_s3_class(nv, "geaCV")
  expect_true(is.finite(nv$metrics$mean_R2))
  expect_false(grepl("SiteCode", nv$type, fixed = TRUE))
  # naiveCV and spatialBlockCV differ only in fold *grouping*, not in target
  # variance or model family; guard against a wiring bug (e.g. accidentally
  # grouped, or a broken fold assignment) rather than asserting a direction
  # (a directional claim is not reliable at n = 1 seed / tiny num.trees; see
  # the vignette's multi-seed production comparison for the qualitative
  # pattern).
  expect_lt(abs(nv$metrics$mean_R2 - .sb_fixture$metrics$mean_R2), 0.4)
})

test_that("confoundGapTest (issue #3) returns a coherent gap and null distribution", {
  gt <- .gt_fixture
  expect_s3_class(gt, "geaGapTest")
  expect_equal(gt$observed_gap, gt$interp_R2 - gt$loco_R2, tolerance = 1e-9)
  expect_equal(length(gt$null_gap), 3)
  expect_true(gt$p_value > 0 && gt$p_value <= 1)
})
