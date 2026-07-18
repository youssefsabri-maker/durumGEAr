# durumGEAr

Confound-aware genotype-environment association (GEA) modelling for genebank
accession data. `durumGEAr` packages a defensible statistical workflow for
data that is both **pseudoreplicated** (many accessions collected at the same
site share literally identical bioclimatic predictors) and **spatially
confounded** (country of origin can substitute for climate, inflating
apparent predictive skill). It was developed on a durum wheat (*Triticum
turgidum* ssp. *durum*) genebank collection of 1,060 confident effective
units across 43 countries, with predictors BIO1-BIO19 + altitude +
coordinates and targets five genetic PCA scores (gPC1-gPC5), following the
function-naming conventions of the `icardaFIGSr` package.

**Headline finding:** the workflow's own diagnostics show why "predictive
skill" numbers for this kind of data need to be reported with real care. A
naive, site-grouped interpolation cross-validation gives a per-target-mean R2
around **0.35-0.36** - but leave-one-country-out extrapolation collapses to
a mean R2 of **-0.07** (per-country range: -0.44 to +0.40, with only 21/43
countries exceeding the global-mean baseline), a gap of roughly
**0.3-0.4 R2**. That gap is statistically significant (*p* < 0.001; see
`confoundGapTest()`) and is the signature of a model partly memorizing
country identity rather than learning a transferable climate relationship.
Fold-safe residualization then shows the picture is target-dependent: some
gPC axes (e.g. **gPC2**: *p* = 0.02, residual R² = +0.12) retain real,
statistically significant within-country climate signal after country identity
is removed, while others are largely a country-identity artifact (e.g. **gPC5**:
*p* = 0.78, residual R² = -0.05). See `vignette("durumGEAr-workflow")` for the
full, reproducible walkthrough with computed numbers and significance tests.

## Installation

```r
# from a local source checkout
devtools::install(".")

# or, from a built tarball
install.packages("durumGEAr_1.3.0.tar.gz", repos = NULL, type = "source")
```

## Why This Package Exists: The Pseudoreplication & Confounding Problem

Genebank accessions are collected at specific geographic locations (sites) across countries. This creates two statistical challenges:

1. **Pseudoreplication**: Many accessions share the same site, so their bioclimatic predictors are literally identical (e.g., 50 durum accessions from the same site in Ethiopia). Naive random CV can place twins in both train and test, inflating apparent skill.

2. **Spatial Confounding**: Country of origin can substitute for climate. A model that "predicts" high PC1 scores for accessions from hot countries isn't learning about heat-tolerance; it's memorizing geography. This inflates interpolation R² but collapses on extrapolation to new countries.

The durum wheat collection exemplifies this: **1,060 effective units (after pseudoreplication collapse) across 43 countries with 9× Kish effective sample-size inflation**. The naive interpolation R² of 0.35 collapses to −0.07 on leave-one-country-out CV—a 0.4 R² gap that is the signature of confounding, not noise.

`durumGEAr` provides a defensible workflow to: (1) collapse pseudoreplicates, (2) honestly quantify the confound gap, (3) remove country identity and test what climate signal remains, and (4) deploy only targets certified as carrying real, transferable signal.

---

## Minimal usage example

```r
library(durumGEAr)
data(durumUnits)

# Honest interpolation skill (grouped by site, never pooled across targets)
sb <- spatialBlockCV(durumUnits)
sb$metrics$mean_R2

# Harsh extrapolation stress test (leave-one-country-out)
lo <- locoCV(durumUnits, group = "Country")
lo$metrics$mean_R2

# Is the gap between them real, or noise? (permutation test)
gt <- confoundGapTest(durumUnits, n_perm = 20, num.trees = 100)
gt$p_value
```

## Function Reference

### Data Preparation & Cleaning

#### `collapseUnits(data, ...)`
Collapses pseudoreplicated rows (accessions from the same site with identical predictors) into independent effective units. This function:
- Identifies exact duplicates and cloned target vectors
- Handles intra-group (IG) duplicates and admixture confidence filtering
- Returns a "confident" modelling frame (one row per site × genetic cluster combination)
- Produces an audit ledger documenting all filtering decisions

**Usage:**
```r
data(durumRaw)  # raw 3,428 rows (many duplicates)
cu <- collapseUnits(durumRaw)
cu$confident    # 1,060 confident effective units
cu$audit        # full filtering ledger
```

**Key output:** `cu$confident` is the canonical modelling frame used by all downstream functions.

---

#### `mapAccessions(data, ...)`
Maps accessions to their genetic cluster assignments. Useful for linking raw accession IDs back to the genetic PCA scores and cluster membership after collapsing.

---

### Environmental Dimensionality Reduction

#### `computeEPC(data, predictors = c("BIO1", ..., "Altitude"), ...)`
Performs principal component analysis (PCA) on the bioclimatic predictors (BIO1-BIO19 + altitude). Returns:
- `scores`: the rotated environmental PCA scores
- `loadings`: variable contributions to each PC
- `var_explained`: proportion of variance explained per axis

**Usage:**
```r
ep <- computeEPC(durumUnits)
ep$var_explained       # % variance per ePC axis
plot(ep$scores[, 1], ep$scores[, 2])  # visualize climate space
```

---

### Cross-Validation: Honest Interpolation & Extrapolation Estimates

#### `spatialBlockCV(data, predictors, targets, group = "SiteCode", k = 5, seeds = 42:46, num.trees = 600, ...)`
**Interpolation skill**: Evaluates multi-output regression under grouped k-fold cross-validation where entire spatial groups (sites) stay on one side of train/test splits. This prevents pseudo-replication from inflating skill estimates.

**Why this matters:** Random CV can place twins from the same site in both train and test. Grouping forces prediction of genuinely unseen locations.

**Returns:** out-of-fold predictions, per-target R², mean R², and per-seed breakdown.

**Usage:**
```r
sb <- spatialBlockCV(durumUnits, seeds = 42:46, k = 5)
sb$metrics$mean_R2              # ~0.35-0.36 (honest interpolation)
sb$metrics$per_target$R2        # R² per genetic target
```

**Key insight:** This is the *optimistic* estimate (you're predicting new sites in *known* countries).

---

#### `naiveCV(data, predictors, targets, k = 5, seeds = 42:46, ...)`
**Leaky baseline**: Runs ungrouped random k-fold CV (accessions from the same site *can* appear on both sides of a fold). Quantifies how much apparent skill the site-level pseudoreplication alone provides.

**Usage:**
```r
nv <- naiveCV(durumUnits)
nv$metrics$mean_R2              # inflated by pseudoreplication
# Compare to spatialBlockCV: the difference shows pseudoreplication bias
```

**Caution:** This is a diagnostic baseline, not a recommended CV strategy. Use `spatialBlockCV()` for honest interpolation estimates.

---

#### `locoCV(data, predictors, targets, group = "Country", seed = 42, floor = -1, num.trees = 600, ...)`
**Extrapolation skill**: Evaluates regression under leave-one-group-out (LOGO) cross-validation, holding out one entire country at a time. This is the harsh, honest test of whether models generalize to unseen geography.

**Why this matters:** Interpolation asks "can you predict a new site in a country you know?"; LOGO asks "can you predict an entirely new country you've never seen?". A model memorizing country identity scores well on the former and collapses on the latter.

**Returns:** per-target R² per held-out country, per-country comparison to global-mean baseline, and overall mean R².

**Usage:**
```r
lo <- locoCV(durumUnits, group = "Country")
lo$metrics$mean_R2              # ~-0.07 to +0.03 (harsh extrapolation)
lo$per_group                    # R² breakdown per country
lo$n_beat_baseline              # how many countries beat no-skill?
```

**Key insight:** The gap between `spatialBlockCV` (~0.35) and `locoCV` (~-0.07) is the signature of geographic confounding.

---

### Confound Diagnosis & Quantification

#### `confoundGapTest(data, predictors, targets, group = "Country", n_perm = 20, seed = 42, num.trees = 100, ...)`
**Gap significance test**: Quantifies whether the interpolation/extrapolation R² gap is statistically significant (i.e., is confounding real, or just noise?).

**Procedure:**
1. Computes observed gap: `interpolation_R² - extrapolation_R²` (e.g., 0.35 − (−0.07) = 0.42).
2. Shuffles targets within each country (preserves country structure exactly, destroys true climate-target links).
3. Recomputes the gap on permuted data (builds null distribution).
4. Compares observed gap to null; reports p-value.

**Returns:** observed gap, null distribution, p-value, and interpretation.

**Usage:**
```r
gt <- confoundGapTest(durumUnits, n_perm = 100, num.trees = 600)
gt$p_value                      # p < 0.05 = gap is real confounding
gt$observed_gap                 # point estimate of gap size
```

**Interpretation:** If `p < 0.05`, the confounding signal is statistically significant. The workflow is justified.

---

#### `residualize(data, predictors, targets, country = "Country", group = "SiteCode", k = 5, seeds = 42:46, shrink = FALSE, num.trees = 600, n_perm = 100, alpha = 0.05, ...)`
**Fold-safe confound removal**: Removes country-identity signal (by subtracting country means, fit on training data only) and asks: "does climate still predict within-country variation?"

**Key design:** For each fold, country means are computed from training rows only, then subtracted from both train and test. An ExtraTrees model is refit on residuals. This avoids leakage.

**Permutation significance:** Targets are shuffled *within each country* (preserving country structure exactly), the identical residualization + fit is rerun, and a null distribution of residual R² is built. Verdict requires *both* `p < alpha` *and* `residual_R² > 0`.

**Returns:** 
- `per_target_R2`: Multi-seed-averaged residual R² per target.
- `obs_R2_test`: Single-seed statistic used for permutation test (may differ in sign from multi-seed average).
- `p_value`: Permutation p-value per target.
- `verdict`: "real within-country signal" or "country-identity artifact" (per target).
- `verdict_fragility`: "robust", "marginal", or "fragile" (how shaky is the verdict?).

**Usage:**
```r
rz <- residualize(durumUnits, seeds = 42:46, n_perm = 100, num.trees = 600)
rz$verdict                      # which targets have real within-country signal?
rz$p_value                      # significance per target
rz$verdict_fragility            # robustness of each verdict
```

**Interpretation example:**
- `gPC2: p = 0.02, obs_R2 = +0.12 → "real within-country signal" (robust)` — deploy this target.
- `gPC5: p = 0.78, obs_R2 = -0.05 → "country-identity artifact"` — do not deploy.
- `gPC3: p = 0.04, obs_R2 = +0.005 → "real within-country signal" (marginal)` — rerun at higher n_perm to confirm.

**Optional shrinkage:** Set `shrink = TRUE` to apply empirical-Bayes shrinkage of small-country means (James-Stein BLUP), stabilizing estimates for countries with few accessions.

---

### Driver Analysis: Which Climate Variables Matter?

#### `driverAnalysis(data, residualize_result, predictors, targets, n_perm = 100, num.trees = 600, ...)`
**Permutation-based predictor importance**: For each target, shuffles each predictor (one at a time) and measures the drop in residual R². Determines which climate variables are most influential.

**Returns:** importance scores (R² drop per shuffled predictor) and permutation p-values per target.

**Usage:**
```r
da <- driverAnalysis(durumUnits, residualize_result = durum_residualize_results)
head(da$importance)             # which BIO variables matter most?
plot(da$importance)             # heatmap of importance per target
```

**Interpretation:** Large R² drops = important variables. Permutation p-value < 0.05 = statistically significant.

---

### Deployment & Prediction

#### `fitGeneticScoreModel(data, residualize_result, predictors, targets, num.trees = 600, ...)`
Fits one ExtraTrees predictor per genetic target, inheriting the reliability verdict from `residualize()`. Gates ("trusts") only targets with "robust" verdicts; others are flagged as unreliable.

**Returns:** a fitted `geneticScoreModel` object with per-target models and gating metadata.

**Usage:**
```r
mod <- fitGeneticScoreModel(durumUnits, residualize_result = durum_residualize_results)
mod                             # prints per-axis trusted/gated labels
```

---

#### `predict.geneticScoreModel(object, newdata, ...)`
Generates predictions on new data. Automatically flags predictions for gated (untrusted) targets.

**Returns:** matrix of predictions (n_new × n_targets) with attribute `gated_targets` listing axes that failed residualization.

**Usage:**
```r
pr <- predict(mod, durumUnits[1:5, ])
pr                              # predictions for the first 5 accessions
attr(pr, "gated_targets")       # which axes should you ignore?
```

**Key safeguard:** Gated targets can never be silently presented as confident. Users must explicitly acknowledge which axes are certified.

---

#### `robustGeneticScore(mod, newdata, keep = NULL, ...)`
Allows force-keeping targets that `residualize()` had gated, if you have **strong external evidence** that they carry real signal. 

**WARNING:** Overriding a "country-identity artifact" verdict risks silent deployment of a geographic coincidence, not a transferable climate relationship. Only use this function if: (1) the target has been validated independently on a separate geographic region, *and* (2) you have documented this external evidence in your model report.

**Usage:**
```r
# Deploy only gPC2 (ignore gPC1, gPC3, etc., even if fitted)
pr <- robustGeneticScore(mod, durumUnits[1:5, ], keep = c("gPC2"))

# CAUTION: Only do this if you have external validation that gPC2 truly transfers.
# Document your reasoning before shipping this model.
```

---

#### `checkExtrapolationRisk(mod, newdata, ...)`
Flags sites whose bioclimatic predictors fall outside the convex hull of the training climate space. Warns when you're extrapolating beyond the observed climate envelope.

**Usage:**
```r
risk <- checkExtrapolationRisk(mod, durumUnits[1:5, ])
risk$extrapolation              # which sites are risky?
```

---

### Genetic Clustering (Stage 2)

#### `scoreThenCluster(data, genetic_targets, residualize_result = NULL, per_tree = FALSE, fit_final = TRUE, ...)`
Chains genetic PCA scores through a clustering algorithm (ExtraTrees or hierarchical) to assign accessions to discrete genetic groups. Optionally uses only targets certified as "robust" by residualize().

**Returns:** cluster assignments and a fitted classifier (if `fit_final = TRUE`).

**Usage:**
```r
s2 <- scoreThenCluster(durumUnits, per_tree = FALSE, fit_final = TRUE)
s2$cluster                      # genetic cluster per accession
s2$final_classifier             # fitted model for predicting clusters on new data
```

---

#### `predictCluster(mod, newdata, stage2_classifier, ...)`
Predicts genetic cluster membership for new accessions, given their climate scores and the fitted Stage-2 classifier.

**Usage:**
```r
pr_cluster <- predictCluster(mod, durumUnits[1:5, ], s2$final_classifier)
pr_cluster                      # predicted genetic cluster per new accession
```

---

### Utility Functions

#### `getMetrics(observed, predicted, ...)`
Computes per-target R² and RMSE, averaging across multi-output predictions.

**Usage:**
```r
metrics <- getMetrics(Y_observed, Y_predicted)
metrics$per_target              # R² and RMSE per target
metrics$mean_R2                 # overall mean R²
```

## Included Datasets

The package ships with pre-processed durum wheat data and cached results to enable fast exploration without waiting for long permutation runs.

#### `durumRaw`
The raw input data: 3,428 rows of accession-level genotype and environment data before pseudoreplication collapsing.

```r
data(durumRaw)
dim(durumRaw)   # 3,428 rows × many columns
```

#### `durumUnits`
The canonical modelling frame: 1,060 rows (confident effective units) after collapsing pseudoreplicates. This is the dataset used throughout the vignette and is the recommended starting point for exploration.

```r
data(durumUnits)
dim(durumUnits)  # 1,060 rows × columns (SiteCode, Country, BIO1-BIO19, Altitude, gPC1-gPC5)
```

#### `durum_residualize_results`
Pre-computed output from `residualize()` at production settings (n_perm=100, num.trees=600, seeds=42:46). Shipped to avoid re-running the 2-hour permutation test every time you want to explore deployment.

```r
data(durum_residualize_results)
durum_residualize_results$verdict          # which targets pass residualization?
durum_residualize_results$p_value
durum_residualize_results$verdict_fragility
```

#### `durum_loco_results`
Pre-computed output from `locoCV()` at production settings. Useful for quick inspection of extrapolation performance per country.

```r
data(durum_loco_results)
durum_loco_results$per_group               # R² per held-out country
```

---

## Workflow: Quick Start

Here's a step-by-step walkthrough for using durumGEAr on your own data:

### Step 1: Load and Inspect Your Data

```r
library(durumGEAr)

# Load your data (should have: SiteCode, Country, BIO1-BIO19, Altitude, Latitude, Longitude, + genetic targets)
my_data <- read.csv("my_genebank_data.csv")

# Collapse pseudoreplicates
cu <- collapseUnits(my_data)
modelling_frame <- cu$confident
print(cu)                       # inspect the filtering audit
```

### Step 2: Check for Confounding (Confound Gap Diagnostic)

```r
# Run spatialBlockCV (interpolation) and locoCV (extrapolation) at light settings for speed
sb <- spatialBlockCV(modelling_frame, seeds = 42, num.trees = 100, k = 5)
lo <- locoCV(modelling_frame, seed = 42, num.trees = 100)

cat("Interpolation R²:", round(sb$metrics$mean_R2, 3), "\n")
cat("Extrapolation R²:", round(lo$metrics$mean_R2, 3), "\n")
cat("Gap:", round(sb$metrics$mean_R2 - lo$metrics$mean_R2, 3), "\n")

# Is the gap statistically significant?
gt <- confoundGapTest(modelling_frame, n_perm = 20, num.trees = 50)
cat("Confound gap p-value:", gt$p_value, "\n")

# If p < 0.05, confounding is real—proceed to residualization.
```

### Step 3: Diagnose Within-Country Signal (Residualization at Production Settings)

**Warning:** This step is slow (30–40 minutes for ~1,000 accessions). Run once, cache the result.

```r
# Fit residualization at full settings
rz <- residualize(modelling_frame, seeds = 42:46, n_perm = 100, num.trees = 600)

# Which targets have real within-country climate signal?
print(rz)                       # detailed output with verdicts and fragility grades
rz$verdict                      # per-target verdict
rz$p_value                      # per-target significance
rz$verdict_fragility            # per-target robustness
```

### Step 4: Deploy Trusted Targets

Only deploy targets with "robust" verdicts. Optionally, use `robustGeneticScore()` to force-keep marginal targets if you have external confidence.

```r
# Fit the deployment model
mod <- fitGeneticScoreModel(modelling_frame, residualize_result = rz)
print(mod)                      # shows which targets are trusted vs. gated

# Predict on new accessions
new_accessions <- modelling_frame[1:10, ]  # example: first 10 accessions
predictions <- predict(mod, new_accessions)

# Check which targets are gated (unreliable)
gated <- attr(predictions, "gated_targets")
cat("Gated targets (do not report):", paste(gated, collapse = ", "), "\n")

# Safely deploy only non-gated targets
safe_predictions <- predictions[, setdiff(colnames(predictions), gated)]
```

### Step 5: Understand Which Climate Variables Matter

```r
# Permutation-based driver analysis
da <- driverAnalysis(modelling_frame, residualize_result = rz, 
                     num.trees = 600, n_perm = 100)

# Which BIO variables are most important per target?
print(da$importance)

# Visualize
plot(da, main = "Climate driver importance per genetic axis")
```

---

## Common Pitfalls & How to Avoid Them

**Pitfall 1: Treating naive R² as the true skill estimate.**  
→ Always report *both* interpolation (spatialBlockCV) and extrapolation (locoCV) R². The gap tells you how much confounding matters.

**Pitfall 2: Rerunning residualize() with very small n_perm and expecting a "strong" verdict.**  
→ The permutation p-value has a resolution floor of 1/(n_perm+1). At n_perm=20, the smallest possible p-value is ~0.048. Rerun with n_perm=100 or higher before trusting a verdict marked "marginal" or "fragile".

**Pitfall 3: Forcing deployment of targets marked "country-identity artifact".**  
→ Use `robustGeneticScore(..., keep = ...)` judiciously and *only* if you have external evidence (e.g., a published trait known to vary within-country in the same direction). Do not override verdicts on a whim.

**Pitfall 4: Assuming results on durum wheat transfer to other crops.**  
→ Re-run the full diagnostic pipeline (confoundGapTest → residualize) on your crop. The confound magnitude and target-level verdicts are dataset-specific.

**Pitfall 5: Deploying a gated target without acknowledging its uncertainty.**  
→ Every predicted value should be reported with a gating label: "gPC2 = 2.1 (robust)" vs. "gPC5 = 3.2 (gated—do not rely)".

---

## Limitations

This workflow has only been developed and validated on **one** durum wheat
genebank collection (1,060 effective units, 43 countries, 5 genetic PCA
targets). The specific numeric findings above - the size of the
interpolation/extrapolation gap, which gPC targets carry real within-country
signal, and the Stage-2 clustering accuracy - are properties of this dataset
and should not be assumed to transfer to other crops, other marker sets, or
other genebank collections without re-running the full diagnostic pipeline.
The statistical *methodology* (pseudoreplication collapsing, spatial-block
and leave-one-group-out CV, fold-safe residualization, permutation
significance testing) is general, but every number the package prints is
specific to the data it was run on.

Verdicts from `residualize()` and `confoundGapTest()` are based on permutation
p-values; always check the sign and magnitude of the reported R²/gap before
treating a target as carrying real signal, and treat any p-value at or near the
theoretical floor (1/(n_perm+1)) as low-resolution rather than strongly
significant — rerun with a higher `n_perm` to confirm.

## Interpreting Results: Verdict & Fragility Grades

`residualize()` returns not just a p-value, but a verdict and a fragility grade per target. Understanding both is key:

### Verdicts

**"real within-country signal"**  
- Means: This target's climate-to-genetics relationship is genuine and transferable.
- Action: Safe to deploy via `fitGeneticScoreModel()`.
- Example: gPC2 in durum wheat (p=0.02, obs_R²=+0.12).

**"country-identity artifact"**  
- Means: This target's apparent relationship with climate is just a proxy for country of origin; no transferable signal.
- Action: Do NOT deploy. Any apparent skill is geographic coincidence.
- Example: gPC5 in durum wheat (p=0.78, obs_R²=−0.05).

### Fragility Grades

The verdict's confidence is graded as "robust", "marginal", or "fragile":

**Robust**  
- Verdict is stable and reliable.
- p-value well above floor (1/(n_perm+1)).
- Residual R² magnitude is substantial (|R²| ≥ 0.01).
- Sign of obs_R2_test agrees with multi-seed-averaged R².
- Action: Fully trust this verdict.

**Marginal**  
- Verdict is borderline or weakly supported.
- Either: p-value near floor, *or* obs_R² near zero, *or* multi-seed uncertainty large.
- Action: Rerun at higher n_perm and/or more seeds to confirm. Or accept the uncertainty and handle downstream with caution.

**Fragile**  
- Verdict is unstable (sign of obs_R²_test disagrees with multi-seed average).
- Same target can flip between "signal" and "artifact" depending on random seed.
- Action: Do NOT use for deployment without substantial additional evidence. Rerun with more seeds/permutations.

---

## Full Details & Case Study Results

For the complete, reproducible walkthrough on durum wheat including:
- Exact R² numbers and p-values per target
- Per-country extrapolation performance
- Genetic cluster assignments (Stage 2)
- Climate driver rankings

See `vignette("durumGEAr-workflow")`.

```r
vignette("durumGEAr-workflow")
```

This vignette runs the entire pipeline (caching expensive permutation steps) and displays all computed numbers.

---

## Citation

If you use `durumGEAr` in published research, please cite:

```
Youssef, S. & Kehel, Z. (2026). durumGEAr: Confound-aware genotype-environment 
association modelling for genebank accessions. R package version 1.3.0.
https://github.com/ICARDA-org/durumGEAr
```

---

## Research Context & References

This workflow was developed to address a specific, underappreciated problem in crop genomics: spatial confounding in genotype-environment association (GEA) studies on genebank collections. The package implements:

- **Pseudoreplication collapsing** following Hurlbert (1984) and principles of independent effective units in plant population genetics.
- **Spatial-block cross-validation** following Roberts et al. (2017), preventing geographic leakage across train/test splits.
- **Leave-one-country-out (LOGO) validation** as a harsh, honest extrapolation test distinct from interpolation.
- **Fold-safe residualization** to isolate within-country climate signal from between-country identity confounds, with permutation significance testing (within-country target shuffling preserves the confound structure under the null).

The method was validated on a durum wheat genebank collection (1,060 accessions, 43 countries), but the statistical methodology is general and applicable to any GEA dataset with spatial pseudoreplication and confounding.

**Key references:**
- Hurlbert, S. H. (1984). Pseudoreplication and the design of ecological field experiments. *Ecological Monographs*, 54(2), 187–211.
- Roberts, D. R., et al. (2017). Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. *Ecography*, 40(8), 913–929.

---

## Acknowledgment

Developed under the supervision of Dr. Zakaria Kehel (ICARDA).
