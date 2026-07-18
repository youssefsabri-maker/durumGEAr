test_that("collapseUnits reproduces the canonical cleaning audit", {
  data(durumRaw)
  cu <- collapseUnits(durumRaw)
  expect_s3_class(cu, "unitCollapse")
  a <- cu$audit
  expect_equal(a$n_raw, 4545L)
  expect_equal(a$n_exact_duplicates, 0L)
  expect_equal(a$n_cloned_gpc, 58L)
  expect_equal(a$n_confident, 1060L)
  expect_equal(a$n_countries_confident, 43L)
  # predictors must be constant within each (SiteCode x Cluster) group
  expect_true(a$constancy_pass)
  expect_lt(a$max_within_group_predictor_sd, 1e-6)
})

test_that("confident frame is a subset of all units and has no all-cloned units", {
  data(durumRaw)
  cu <- collapseUnits(durumRaw)
  expect_lte(nrow(cu$confident), nrow(cu$units))
  expect_true(all(cu$confident$gpc_missing_all_cloned == 0L))
})

test_that("missing required columns raise an informative error", {
  bad <- data.frame(SiteCode = 1:3)
  expect_error(collapseUnits(bad), "Missing required columns")
})
