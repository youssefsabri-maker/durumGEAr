test_that("getMetrics returns per-target-mean R2, not pooled", {
  set.seed(1)
  obs <- matrix(rnorm(500), 100, 5); colnames(obs) <- paste0("gPC", 1:5)
  pred <- obs + matrix(rnorm(500, sd = 0.5), 100, 5)
  m <- getMetrics(obs, pred)
  expect_s3_class(m, "geaMetrics")
  expect_equal(nrow(m$per_target), 5)
  expect_equal(m$mean_R2, mean(m$per_target$R2), tolerance = 1e-12)
  # perfect prediction -> R2 == 1 on every target
  mp <- getMetrics(obs, obs)
  expect_equal(mp$mean_R2, 1)
})

test_that("getMetrics R2 differs from a pooled R2 when target variances differ", {
  set.seed(2)
  y1 <- rnorm(200, sd = 10); y2 <- rnorm(200, sd = 0.1)
  obs <- cbind(y1, y2)
  pred <- cbind(y1 + rnorm(200, sd = 5), y2 + rnorm(200, sd = 0.09))
  m <- getMetrics(obs, pred)
  pooled <- 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
  expect_false(isTRUE(all.equal(m$mean_R2, pooled)))
})

test_that("computeEPC is leakage-safe: predict on training rows reproduces scores", {
  data(durumUnits)
  ep <- computeEPC(durumUnits, n_comp = 5)
  expect_s3_class(ep, "ePCfit")
  expect_equal(length(ep$var_explained), 5)
  re <- predict(ep, durumUnits)
  # projecting the same data back must match the fitted scores (up to sign/scale of prcomp)
  expect_equal(abs(stats::cor(re[, 1], ep$scores[, 1])), 1, tolerance = 1e-6)
  expect_true(sum(ep$var_explained) > 0.85)   # 5 comps explain >85% of climate variance
})
