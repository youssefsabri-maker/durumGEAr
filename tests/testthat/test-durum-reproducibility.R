# Reproducibility + structure tests for the v1.2.0 additions (spec issues #13
# per_seed_R2/CIs, #14 fragility, #15 driver sign/correlation).
#
# Philosophy matches test-cv-residualize.R: shared cheap fixtures fit ONCE at
# tiny settings (k = 2, num.trees = 20, small n_perm) and reused across
# test_that() blocks. These test that the new *code* is correct (fields present,
# shapes right, invariants hold, determinism under fixed seed) -- NOT that the
# permutation null is well-calibrated at production settings (that is the
# RUN_VALIDATION_TESTS path and the vignette validation study).

data(durumUnits)

# Shared fixtures: multi-seed so CIs are computable; tiny everything else.
.rz2 <- residualize(durumUnits, seeds = 42:43, k = 2, num.trees = 20, n_perm = 8)
.gt2 <- confoundGapTest(durumUnits, seed = 42, n_perm = 5, num.trees = 15,
                        min.node.size = 5)
.dr2 <- driverAnalysis(durumUnits, "gPC2", seeds = 42:43, n_perm = 5,
                       num.trees = 20)

# ---- §6: per_seed_R2 exposure -------------------------------------------------

test_that("residualize() exposes per_seed_R2 as a seeds x targets matrix", {
  ps <- .rz2$per_seed_R2
  expect_true(is.matrix(ps))
  expect_equal(dim(ps), c(2L, 5L))               # 2 seeds x 5 targets
  expect_equal(colnames(ps), paste0("gPC", 1:5))
  # the headline per_target_R2 must be the column means of per_seed_R2 -- i.e.
  # the exposed matrix is genuinely what feeds the summary, not a separate object
  expect_equal(unname(colMeans(ps)), unname(.rz2$per_target_R2), tolerance = 1e-9)
})

# ---- §1: confidence intervals -------------------------------------------------

test_that("residualize() returns per-seed R2 CI and bootstrap p-value CI", {
  ci <- .rz2$per_target_R2_CI
  expect_equal(dim(ci), c(2L, 5L))
  expect_equal(rownames(ci), c("2.5%", "97.5%"))
  # lower bound <= upper bound for every target
  expect_true(all(ci["2.5%", ] <= ci["97.5%", ]))
  # the per-seed CI must bracket the per-seed values it was computed from
  ps <- .rz2$per_seed_R2
  for (j in 1:5) {
    expect_gte(max(ps[, j]) + 1e-9, ci["97.5%", j])
    expect_lte(min(ps[, j]) - 1e-9, ci["2.5%", j])
  }
  pci <- .rz2$p_value_CI
  expect_equal(dim(pci), c(2L, 5L))
  expect_true(all(pci >= 0 & pci <= 1))
  expect_true(all(pci["2.5%", ] <= pci["97.5%", ]))
})

test_that("confoundGapTest() returns bootstrap CIs for gap and p-value", {
  expect_length(.gt2$observed_gap_CI, 2L)
  expect_true(.gt2$observed_gap_CI[1] <= .gt2$observed_gap_CI[2])
  # the observed gap sits inside its own CI (CI is centred on it by construction)
  expect_gte(.gt2$observed_gap, .gt2$observed_gap_CI[1] - 1e-9)
  expect_lte(.gt2$observed_gap, .gt2$observed_gap_CI[2] + 1e-9)
  expect_length(.gt2$p_value_CI, 2L)
  expect_true(all(.gt2$p_value_CI >= 0 & .gt2$p_value_CI <= 1))
})

test_that("driverAnalysis() returns per-driver CI columns and top-driver CI", {
  u <- .dr2$univariate
  expect_true(all(c("R2_CI_low", "R2_CI_high", "correlation", "sign") %in% names(u)))
  ok <- is.finite(u$R2_CI_low) & is.finite(u$R2_CI_high)
  expect_true(all(u$R2_CI_low[ok] <= u$R2_CI_high[ok]))
  expect_length(.dr2$observed_R2_CI, 2L)
  expect_true(is.matrix(.dr2$per_seed_uni))
  expect_equal(nrow(.dr2$per_seed_uni), 2L)       # one row per seed
})

# ---- §1 resolution note: p-value CI cannot be finer than the null sample ------

test_that("residualize() p-value CI resolution is capped by n_perm", {
  # With n_perm = 8, the finest achievable p-value step is 1/(n_perm+1). The
  # bootstrap CI is drawn from a null of only n_perm rows, so its width cannot
  # imply a resolution finer than that grid -- guard against a future change that
  # silently makes the CI look artificially precise.
  step <- 1 / (.rz2$n_perm + 1)
  pci <- .rz2$p_value_CI
  widths <- pci["97.5%", ] - pci["2.5%", ]
  # Each bootstrap p-value is (count + 1)/(n_perm + 1), so both CI endpoints must
  # land exactly on that discrete grid -- a smaller effective step would mean the
  # CI is being reported at a finer resolution than the n_perm-row null supports.
  on_grid <- function(p) abs(p * (.rz2$n_perm + 1) - round(p * (.rz2$n_perm + 1))) < 1e-6
  expect_true(all(on_grid(as.numeric(pci))))
  # and any non-zero width is at least one grid step (never a sub-step artefact)
  nz <- widths[widths > 1e-9]
  if (length(nz)) expect_true(all(nz >= step - 1e-9))
  expect_true(all(widths >= 0))   # always at least one assertion runs
})

# ---- §2: fragility scoring ----------------------------------------------------

test_that(".fragility_score() grades the three shakiness conditions correctly", {
  # robust: healthy effect, mid p, signs agree
  expect_equal(.fragility_score(0.08, 0.02, 0.07, n_perm = 100), "robust")
  # marginal: tiny effect (|R2| < 0.01)
  expect_equal(.fragility_score(0.005, 0.02, 0.005, n_perm = 100), "marginal")
  # marginal: floor p-value
  expect_equal(.fragility_score(0.08, 1/101, 0.07, n_perm = 100), "marginal")
  # fragile: sign disagreement between single-seed and multi-seed R2 (overrides)
  expect_equal(.fragility_score(0.08, 0.02, -0.05, n_perm = 100), "fragile")
  # fragile wins even when a marginal condition also holds
  expect_equal(.fragility_score(0.005, 1/101, -0.05, n_perm = 100), "fragile")
})

test_that("residualize() attaches a per-target fragility grade from that logic", {
  vf <- .rz2$verdict_fragility
  expect_length(vf, 5L)
  expect_true(all(vf %in% c("robust", "marginal", "fragile")))
  expect_equal(names(vf), paste0("gPC", 1:5))
  # the attached grade must equal recomputing .fragility_score from the object's
  # own fields -- i.e. residualize() calls the shared helper, no inline drift.
  for (j in 1:5) {
    expect_equal(unname(vf[j]),
                 .fragility_score(.rz2$obs_R2_test[j], .rz2$p_value[j],
                                  .rz2$per_target_R2[j], .rz2$n_perm))
  }
})

test_that("driverAnalysis() attaches a fragility grade", {
  expect_true(.dr2$verdict_fragility %in% c("robust", "marginal", "fragile"))
})

# ---- §8: driver sign / correlation --------------------------------------------

test_that("driverAnalysis() sign column matches the correlation column", {
  u <- .dr2$univariate
  ok <- is.finite(u$correlation)
  expect_equal(u$sign[ok & u$correlation >= 0], rep("positive", sum(ok & u$correlation >= 0)))
  expect_equal(u$sign[ok & u$correlation < 0],  rep("negative", sum(ok & u$correlation < 0)))
  # $correlation / $sign named vectors cover every predictor tested
  expect_length(.dr2$correlation, .dr2$n_predictors_tested)
  expect_length(.dr2$sign, .dr2$n_predictors_tested)
})

# ---- Reproducibility: identical seeds -> identical output ----------------------

test_that("residualize() is deterministic under fixed seeds", {
  a <- residualize(durumUnits, seeds = 42:43, k = 2, num.trees = 20, n_perm = 8)
  expect_equal(a$per_target_R2, .rz2$per_target_R2, tolerance = 1e-9)
  expect_equal(a$p_value, .rz2$p_value, tolerance = 1e-9)
  expect_equal(a$per_target_R2_CI, .rz2$per_target_R2_CI, tolerance = 1e-9)
  expect_equal(a$p_value_CI, .rz2$p_value_CI, tolerance = 1e-9)  # boot seed fixed
})

test_that("confoundGapTest() bootstrap CI is deterministic under fixed seed", {
  b <- confoundGapTest(durumUnits, seed = 42, n_perm = 5, num.trees = 15,
                       min.node.size = 5)
  expect_equal(b$observed_gap_CI, .gt2$observed_gap_CI, tolerance = 1e-9)
  expect_equal(b$p_value_CI, .gt2$p_value_CI, tolerance = 1e-9)
})

# ---- Backward compatibility: v1.1.1 fields still present -----------------------

test_that("v1.1.1 geaResidual fields survive the v1.2.0 additions", {
  for (f in c("metrics", "per_target_R2", "obs_R2_test", "p_value",
              "null_mean", "null_sd", "n_perm", "alpha", "verdict", "shrink"))
    expect_true(f %in% names(.rz2), info = f)
})
