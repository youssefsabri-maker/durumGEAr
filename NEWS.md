# durumGEAr 1.3.0

This release adds a **deployment layer**: functions that turn the validated
diagnostics into predictions you can act on, while keeping each `residualize()`
reliability verdict attached to every number so an uncertified genetic axis can
never be silently presented as a confident prediction. Issue numbering continues
from v1.2.0 (which ended at #19). No breaking changes to existing functions.

## New features

### Issue #20: Deployable per-target genetic-score model
- **What**: There was no supported way to turn the validated workflow into
  out-of-sample predictions; users had to wire `ranger` calls by hand and lost
  the reliability verdicts in the process.
- **Change**: New `fitGeneticScoreModel()` fits one ExtraTrees (`ranger`)
  regressor per genetic PC on all supplied rows (validation having already
  happened in `spatialBlockCV()`/`residualize()`), and — when passed a
  `residualize_result` — labels each axis **trusted** or **gated**. An axis is
  gated when its residualization verdict is `"country-identity artifact"` OR its
  fragility grade is `fragile`/`marginal`; only a `robust` grade is trusted.
  `print.geneticScoreModel()` shows the per-axis reliability table.

### Issue #21: predict() with reliability-aware gating
- **Change**: `predict.geneticScoreModel()` returns one column per target with
  the reliability table and gated-target vector attached as attributes.
  `gated_to_na = TRUE` blanks untrustworthy columns entirely. Point predictions
  carry no formal predictive interval — documented explicitly.

### Issue #22: Extrapolation-risk guard
- **Change**: New `checkExtrapolationRisk()` flags new sites whose climate lies
  outside the training data two ways — a per-variable quantile envelope
  (`envelope_prob`) and a Mahalanobis distance against the training predictor
  cloud (`mahalanobis_prob`, default 0.975). In-domain flag rate ≈ 2.5%;
  genuinely out-of-domain sites are caught reliably.

### Issue #23: Robust single-score collapsing
- **Change**: New `robustGeneticScore()` collapses the trusted axes into one
  number per observation (`mean`/`median`/`pca1`). The dropped axes come from
  the model's own reliability metadata — **not hardcoded** — so a production
  `n_perm` re-run that promotes an axis from gated to trusted follows
  automatically. Errors when no axis survives gating unless `keep=` is supplied.

### Issue #24: Expose the final Stage-2 classifier + cluster prediction
- **Change**: `scoreThenCluster()` gains `fit_final = FALSE`; when `TRUE`, it
  fits the final QDA classifier on all data and exposes it as `final_classifier`
  (class `geaQDA`). The leakage-free honest accuracy is unchanged. New
  `predictCluster()` chains a `geneticScoreModel` into that classifier for an
  end-to-end climate → genetic-cluster prediction with posterior probability and
  per-row entropy.

### Issue #25: Shipped worked validation objects + practice vignette section
- **Change**: Two lazy-loaded data objects ship as worked examples:
  `durum_residualize_results` (a `residualize()` result) and `durum_loco_results`
  (a `locoCV()` result), both at reduced, documented settings with full help
  pages. A new vignette section, *Using the model in practice*, and a
  regenerable `README.Rmd`/`README.md` demonstrate the full deployment path.

# durumGEAr 1.2.0

This release adds uncertainty quantification, robustness diagnostics, and
validation tooling on top of the v1.1.1 confound-aware workflow. Issue
numbering continues from v1.1.1 (which ended at #12). The two validation studies
(#17, #18) ship with **reduced, documented settings** and cached results; both
scripts are parameterized by environment variables so users can rerun them at
production settings offline.

## New features

### Issue #13: Expose per-seed residualization accuracy
- **What**: `residualize()` previously computed a seeds × targets matrix of
  per-seed accuracies internally but returned only its column means
  (`per_target_R2`), discarding the seed-level detail.
- **Change**: The full matrix is now returned as `$per_seed_R2` (seeds ×
  targets). This is the substrate for the new per-target confidence intervals
  (#14) and is documented in the `@return` section. `colMeans(per_seed_R2)`
  reproduces `per_target_R2` exactly (asserted in the test suite).

### Issue #14: Confidence intervals on every reported statistic
- **What**: Point estimates (`per_target_R2`, `observed_gap`, univariate driver
  R²) were reported without uncertainty.
- **Change**: `residualize()` now returns `$per_target_R2_CI` (2.5%/97.5%
  quantiles of the per-seed accuracies) and `$p_value_CI` (bootstrap of the
  permutation null); `confoundGapTest()` returns `$observed_gap_CI` and
  `$p_value_CI`; `driverAnalysis()` returns per-predictor `R2_CI_low`/`R2_CI_high`
  columns and an `$observed_R2_CI` for the top driver. Print methods show the CI
  lines alongside the existing warnings. **Resolution caveat**: the p-value CIs
  bootstrap the `n_perm`-row permutation null, so they cannot resolve finer than
  the discrete grid `1/(n_perm+1)` — at the default `n_perm = 100` that is
  roughly ±0.01. Bootstrap resamples default to 200.

### Issue #15: Fragility grading of verdicts
- **What**: The v1.1.1 verdict (issue #8: requires both `p < alpha` and
  `obs_R2_test > 0`) is driven by a single-seed statistic and can label a target
  "real within-country signal" even when the multi-seed mean R² is negative — the
  situation that arises for the gPC1/gPC4 country-identity artifacts.
- **Change**: Added a `.fragility_score()` helper and a `$verdict_fragility`
  field grading each verdict as "robust" / "marginal" / "fragile". A verdict is
  **marginal** when |obs_R2| is below 0.01 or its p-value sits on the resolution
  floor, and **fragile** when the sign of the single-seed statistic disagrees
  with the sign of the multi-seed mean (the verdict is then seed-dependent and
  should not be trusted). This *extends* — does not duplicate — the v1.1.1
  floor-p-value and sign-disagreement warnings already in the print methods.
  Analogous fragility flags were added to `confoundGapTest()` and
  `driverAnalysis()`. On the durum data this correctly flags gPC1 and gPC4 as
  fragile while leaving gPC2/gPC5 as the trustworthy real-signal targets.

### Issue #16: Country-level LOCO breakdown
- **What**: `locoCV()` reported only the overall per-target-mean held-out R².
- **Change**: Added `$per_group_target` (a countries × targets matrix of
  held-out R²) plus `$n_beat_baseline` / `$n_groups`, so users can see *which*
  countries drive the extrapolation collapse rather than only the aggregate. Feeds
  the new country-level LOCO figure.

### Issue #17: Permutation-test calibration study
- **What**: The permutation verdicts had no empirical calibration check.
- **Change**: Added `inst/validate_permutation_tests.R`, which measures the
  empirical Type I error of `residualize()` (within-country target shuffle) and
  `confoundGapTest()` (country-label shuffle) under a true null. Ships with a
  reduced-setting cached result; env vars `DURUM_CAL_NREP`, `DURUM_CAL_NPERM`,
  `DURUM_CAL_NTREES`, `DURUM_CAL_K` (and gap-specific variants) allow a
  production rerun.

### Issue #18: Hyperparameter sensitivity grid
- **What**: No systematic check of whether verdicts survive reasonable
  hyperparameter changes.
- **Change**: Added `inst/sensitivity_analysis.R`, sweeping
  `num.trees × min.node.size × k` and recording the observed-gap range and the
  number of distinct verdicts per target (1 = stable). Cached at reduced
  settings; env vars `DURUM_SENS_NTREES`, `DURUM_SENS_NODESIZE`, `DURUM_SENS_K`,
  `DURUM_SENS_NPERM` widen the grid.

### Issue #19: Driver ecological interpretation
- **What**: `driverAnalysis()` reported importance/R² but not the *direction* of
  each predictor's marginal effect.
- **Change**: Added `$correlation` (Pearson of each predictor against the
  fold-safe residual target) and `$sign` ("positive"/"negative") to the output
  and print method, scoped to the VIF-pruned 11-BIO + Altitude predictor set. On
  the durum data the confirmed top drivers are BIO5 (max temperature of the
  warmest month) for gPC2 and BIO15 (precipitation seasonality) for gPC5, by
  joint-model permutation importance.

## Documentation & tooling
- Added `inst/durum_verify.R`, a user-facing verification script that runs the
  full workflow at reduced settings and checks the expected qualitative results
  (interpolation skill, LOCO collapse, gap > 0.25, gPC2/gPC5 real signal,
  gPC1/gPC4 flagged as artifacts by verdict *or* fragility, BIO5/BIO15 drivers).
- Added `inst/FAQ_durum.md` with durum-specific guidance (gap ≈ 0.42, the
  verdict-vs-fragility distinction, resolution-floor p-values, driver caveats).
- Added five publication figures (`man/figures/`, `inst/figures/`): confounding
  gap vs null, residualization recovery with CIs, permutation nulls per target,
  country-level LOCO, and the driver-importance heatmap.
- Vignette gains "Stability under perturbation" (§10) and "Permutation-test
  calibration" (§11) sections that read the cached validation results, plus a
  note that p-value CI resolution is capped by `n_perm`.
- Added `tests/testthat/test-durum-reproducibility.R` covering `per_seed_R2`
  exposure, CI shapes and bracketing, the p-value CI resolution grid, the five
  `.fragility_score()` conditions, driver sign↔correlation agreement, and
  determinism under fixed seeds (including the bootstrap seed).
- `DESCRIPTION` gains `ggplot2` and `gridExtra` under `Suggests` for the figures.

---

# durumGEAr 1.1.1

## Bug fixes

### Issue #8 (CRITICAL): residualize() verdict logic flawed
- **Problem**: `verdict` was based on `p_value < alpha` alone, without checking the sign of the observed R². This caused targets with negative R² (no signal) to be incorrectly labeled "real within-country signal" if they happened to be significant in the permutation test.
- **Fix**: Changed verdict rule to require BOTH `p_value < alpha` AND `obs_R2_test > 0`. Updated roxygen docs and print method to document the two-part rule. Added regression test asserting that negative R² can never receive "real signal" verdict.

### Issue #9 (CRITICAL): Vignette cache permutation counts too low, p-values floor-clamped
- **Problem**: The vignette cache was built with `n_perm=25` for `residualize()` and `n_perm=15` for `confoundGapTest()`, causing all p-values to land exactly on the resolution floor (1/(n_perm+1)). For example, at n_perm=15, the minimum possible p-value is 1/16 ≈ 0.0625, yet the gap test result reported this as if it were a meaningful significance level. This masked instability in the null distributions.
- **Fix**: Rebuilt vignette cache with `n_perm=100` for both functions, raising the minimum resolution to ~0.01. Added warnings to `print.geaResidual()` and `print.geaGapTest()` that alert users when a p-value equals the floor (1/(n_perm+1)), recommending rerun with higher n_perm before treating such results as strongly significant.

### Issue #10: Reproducibility regression — data-raw scripts missing
- **Problem**: v1.1.0 failed to include `data-raw/build_vignette_cache.R`, making the vignette cache unreproducible and the cached numbers unauditable.
- **Fix**: Confirmed both `data-raw/make_data.R` (v1.0.0 holdover) and `data-raw/build_vignette_cache.R` (newly created) are present and committed to the source tree.

### Issue #11 (MINOR): Test comment doesn't match implementation
- **Problem**: The test "driverAnalysis observed and null hyperparameters match (issue #5 fix)" claimed to test hyperparameter threading by varying `num.trees`, but both calls used `num.trees=20`.
- **Fix**: Updated test to call `driverAnalysis()` with `num.trees=20` and `num.trees=80`, then assert that `null_max` distributions differ while `observed_R2` remains the same, properly validating that num.trees is threaded through both observed and null fits.

### Issue #12 (MINOR): Limitations sections missing verdict-interpretation caveat
- **Problem**: README.md and vignette Limitations sections covered dataset-scope limitations but not the statistical limitation from issues #8–#9: that permutation verdicts can be resolution-limited at low n_perm.
- **Fix**: Added one sentence to both Limitations sections: "Verdicts from `residualize()` and `confoundGapTest()` are based on permutation p-values; always check the sign and magnitude of the reported R²/gap before treating a target as carrying real signal, and treat any p-value at or near the theoretical floor (1/(n_perm+1)) as low-resolution rather than strongly significant — rerun with a higher `n_perm` to confirm."

---

# durumGEAr 1.1.0

Initial release of v1.1.0 (prior changelog omitted; refer to v1.0.0 for earlier history).
