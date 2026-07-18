
<!-- README.md is generated from README.Rmd. Please edit that file, then run
     rmarkdown::render("README.Rmd") to regenerate README.md. -->

# durumGEAr

Confound-aware genotype-environment association (GEA) modelling for
genebank accession data. `durumGEAr` packages a defensible statistical
workflow for data that is both **pseudoreplicated** (many accessions
collected at the same site share literally identical bioclimatic
predictors) and **spatially confounded** (country of origin can
substitute for climate, inflating apparent predictive skill). It was
developed on a durum wheat (*Triticum turgidum* ssp. *durum*) genebank
collection of 1,060 confident effective units across 43 countries, with
predictors BIO1-BIO19 + altitude + coordinates and targets five genetic
PCA scores (gPC1-gPC5), following the function-naming conventions of the
`icardaFIGSr` package.

**Headline finding:** the workflow’s own diagnostics show why
“predictive skill” numbers for this kind of data need to be reported
with real care. A naive, site-grouped interpolation cross-validation
gives a per-target-mean R2 around **0.35-0.36** - but
leave-one-country-out extrapolation collapses to around **-0.07 to
+0.03** (near or below no-skill), a gap of roughly **0.3-0.4 R2**. That
gap is the signature of a model partly memorizing country identity
rather than learning a transferable climate relationship. Fold-safe
residualization then shows the picture is target-dependent: some gPC
axes (e.g. gPC2) retain real, statistically significant within-country
climate signal after country identity is removed, while others are
largely a country-identity artifact. See
`vignette("durumGEAr-workflow")` for the full, reproducible walkthrough
with computed numbers and significance tests.

## Installation

``` r
# from a local source checkout
devtools::install(".")

# or, from a built tarball
install.packages("durumGEAr_1.3.0.tar.gz", repos = NULL, type = "source")
```

## Minimal usage example

``` r
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

## Using the fitted model

Once the diagnostics have told you which genetic axes carry transferable
signal, the deployment layer turns that judgement into predictions —
with the reliability verdict kept attached to every number, so an
uncertified axis can never be silently presented as confident.

``` r
library(durumGEAr)
data(durumUnits)
data(durum_residualize_results)     # shipped worked residualize() result

# Fit one predictor per genetic axis; inherit the reliability verdict
mod <- fitGeneticScoreModel(durumUnits,
                            residualize_result = durum_residualize_results)
mod                                 # prints per-axis trusted/gated labels

# Predict, and see which axes the validation does not certify
pr <- predict(mod, durumUnits[1:5, ])
attr(pr, "gated_targets")

# Flag sites whose climate extrapolates beyond the training envelope
checkExtrapolationRisk(mod, durumUnits[1:5, ])$extrapolation

# Chain to genetic-cluster assignment (climate -> scores -> cluster)
s2 <- scoreThenCluster(durumUnits, per_tree = FALSE, fit_final = TRUE)
predictCluster(mod, durumUnits[1:5, ], s2$final_classifier)
```

Gating is conservative: an axis is trusted only when `residualize()`
grades it `robust`, which requires a production-scale permutation run
(`n_perm = 100`, `num.trees = 600`). The shipped
`durum_residualize_results` uses reduced, laptop-runnable settings, so
it gates every axis by default — force-keep the axes you have confirmed
with `robustGeneticScore(mod, newx, keep = c("gPC2", "gPC5"))`. See
`vignette("durumGEAr-workflow")`, section *Using the model in practice*.

## Limitations

This workflow has only been developed and validated on **one** durum
wheat genebank collection (1,060 effective units, 43 countries, 5
genetic PCA targets). The specific numeric findings above - the size of
the interpolation/extrapolation gap, which gPC targets carry real
within-country signal, and the Stage-2 clustering accuracy - are
properties of this dataset and should not be assumed to transfer to
other crops, other marker sets, or other genebank collections without
re-running the full diagnostic pipeline. The statistical *methodology*
(pseudoreplication collapsing, spatial-block and leave-one-group-out CV,
fold-safe residualization, permutation significance testing) is general,
but every number the package prints is specific to the data it was run
on.

Verdicts from `residualize()` and `confoundGapTest()` are based on
permutation p-values; always check the sign and magnitude of the
reported R²/gap before treating a target as carrying real signal, and
treat any p-value at or near the theoretical floor (1/(n_perm+1)) as
low-resolution rather than strongly significant — rerun with a higher
`n_perm` to confirm.

## Acknowledgment

Developed under the supervision of Dr. Zakaria Kehel (ICARDA).
