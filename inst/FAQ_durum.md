# Frequently Asked Questions: durumGEAr Durum Wheat Analysis

These answers are specific to the bundled durum wheat genebank dataset
(`durumUnits`) and the confound-aware workflow in durumGEAr. Numbers quoted are
from real runs of the package's own functions; where a value depends on
settings (`num.trees`, `n_perm`, seeds) that is stated.

## Q: The confounding gap is around 0.4. Is that big?

Yes. `confoundGapTest()` computes `observed_gap = interpolation_R2 - LOCO_R2`.
On durum wheat, per-target-mean interpolation R² is ≈ 0.35 and leave-one-country-out
R² is ≈ −0.07, so the gap is ≈ **0.42** at production settings (it lands a little
lower at reduced `num.trees`). A gap that large means most of the apparent
predictive power comes from *which country* an accession is from, not *what
climate* it experiences. For breeding this is the key warning: you cannot
reliably predict genetic values in a new region from climate alone without first
residualizing out country identity.

## Q: gPC1 and gPC4 are artifacts. Can I ignore them for climate work?

From a climate-adaptation standpoint, yes — their residual within-country signal
is not distinguishable from zero once country identity is removed. They may still
carry other biological meaning (population structure, lineage), so don't discard
them wholesale; just don't use them for climate-based selection.

## Q: The verdict says "real within-country signal" for gPC1/gPC4 but they're supposed to be artifacts. Which is right?

Read the **fragility grade** next to the verdict, not the verdict alone. The
verdict is driven by the single-seed `obs_R2_test` statistic, which can be
positive for a target even when the *multi-seed mean* `per_target_R2` is
negative. When that sign disagreement happens, `residualize()` grades the
verdict **"fragile"** — which is exactly what it does for gPC1 and gPC4. A fragile
verdict means "the sign itself is seed-dependent; do not trust it." So the honest
reading is: gPC2 and gPC5 are robust real signal; gPC1 and gPC4 are fragile and
should be treated as country-identity artifacts. This is why v1.2.0 added
`$verdict_fragility` — the verdict on its own is not enough for the borderline
targets.

## Q: Why do gPC2 and gPC5 have signal but the others don't?

gPC2 is the strongest within-country climate axis (residual R² ≈ 0.18 at
`shrink = FALSE`), gPC5 the second (≈ 0.06). gPC3 sits in between — a weak
residual signal that may or may not survive a stricter permutation test (rerun
with higher `n_perm` to sharpen it). gPC1 and gPC4 are consistently
country-identity artifacts. To distinguish borderline cases, rerun
`residualize()` with more seeds and higher `n_perm`, or collect more
within-country data.

## Q: Should I use num.trees = 600 on my own durum data?

That is the current default, but run the sensitivity grid
(`inst/sensitivity_analysis.R`, §4) on your own data first. The shipped grid
shows how stable the gap and the per-target verdicts are across
`num.trees × min.node.size × k` on *this* dataset; do not assume the same
settings transfer unchanged to a different collection size or country count.

## Q: What is a resolution-floor p-value, and why does the package warn about it?

At `n_perm` permutations the smallest possible p-value is `1 / (n_perm + 1)` —
e.g. `1/101 ≈ 0.0099` at `n_perm = 100`. If a reported p-value sits exactly at
that floor, it means the observed statistic beat *every* permutation: suggestive,
but not finely resolved. The print methods warn when a p-value is at the floor
(NEWS.md issue #9), and v1.2.0's `$verdict_fragility` grades such a verdict at
least "marginal". Rerun with `n_perm ≥ 500` for a sharper estimate before
treating a floor p-value as strong evidence.

## Q: Why does the p-value confidence interval look coarse?

The §1 p-value CI is built by bootstrapping the permutation null matrix, which
only has `n_perm` rows (default 100). Each bootstrap p-value is therefore
`(count + 1)/(n_perm + 1)` — a value on a discrete grid whose spacing is
`1/(n_perm + 1)`. The CI cannot be finer than that grid, so at `n_perm = 100`
it resolves to roughly ±0.01. Raise `n_perm` if you need a finer interval.

## Q: Can I use these climate drivers to breed better durum?

In principle. BIO5 (max temperature of the warmest month) is the top driver of
gPC2, and BIO15 (precipitation seasonality) the top driver of gPC5, by
joint-model permutation importance — both on targets with real within-country
signal. Note this is observational correlation, not a causal or field-validated
result: confirm with known genotypes and targeted stress treatments before
acting on it. Also note that *univariate* residual R² per predictor is noisy at
low settings; the headline driver metric is permutation importance from the
joint model, not univariate R².

## Q: Why is residualize() slow?

It runs `n_perm × n_targets × k` ranger fits at `num.trees = 600` by default —
of order a couple of hours at full defaults. To speed up exploration: use
`n_perm = 50` (quick) / `100` (standard) / `200+` (production), and
`num.trees = 100–150` while iterating. Lock production settings only for the run
whose numbers you report.

## Q: My sensitivity analysis shows a verdict flips across some settings. Should I trust it?

Compare the observed flip count against the number of grid cells the sensitivity
script reports (`verdict_distinct_per_target`: 1 = perfectly stable). A target
whose verdict is stable across the grid is safe to report; one that flips is
genuinely borderline — increase seeds and `n_perm`, or treat it as weak. Report
the real observed flip pattern for your data rather than assuming it matches the
durum case.

## Q: The verdict says "real signal" but the CI overlaps zero. Which do I believe?

The CI. If the 95% interval on `per_target_R2` includes zero, treat the target as
ambiguous regardless of the verdict label — this is the situation
`$verdict_fragility` is designed to catch. Add such a target to the "don't use
for climate work" list.
