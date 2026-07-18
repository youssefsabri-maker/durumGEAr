# Internal: fit one ExtraTrees (ranger, splitrule="extratrees") per target and
# predict on newdata. Returns a matrix (nrow(newdata) x n_targets).
# Not exported.
.fitPredictMO <- function(train, test, predictors, targets,
                          num.trees = 600, min.node.size = 3,
                          num.random.splits = 1, seed = 42) {
  preds <- matrix(NA_real_, nrow = nrow(test), ncol = length(targets),
                  dimnames = list(NULL, targets))
  for (t in targets) {
    df <- data.frame(train[, predictors, drop = FALSE], .y = train[[t]])
    rf <- ranger::ranger(
      dependent.variable.name = ".y",
      data = df,
      num.trees = num.trees,
      min.node.size = min.node.size,
      splitrule = "extratrees",
      num.random.splits = num.random.splits,
      mtry = length(predictors),          # max_features = 1.0
      respect.unordered.factors = "order",
      seed = seed,
      num.threads = 1
    )
    preds[, t] <- stats::predict(rf, data = test[, predictors, drop = FALSE])$predictions
  }
  preds
}

# Internal: build folds grouping whole levels of `group` together.
.groupFolds <- function(group, k = 5, seed = 42) {
  ug <- unique(group)
  set.seed(seed)
  ug <- sample(ug)
  fold_of_group <- stats::setNames(rep(seq_len(k), length.out = length(ug)), ug)
  fold_of_group[as.character(group)]
}
