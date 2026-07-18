#!/usr/bin/env Rscript
#
# durumGEAr v1.2.0: Verify Durum Wheat Analysis Reproducibility
#
# Verifies that durumGEAr produces the expected qualitative results on the
# bundled durum wheat dataset. If all checks pass, the package is correctly
# installed and behaving as documented.
#
# Usage:  Rscript inst/durum_verify.R
#    or:  source(system.file("durum_verify.R", package = "durumGEAr"))
#
# Runtime: a few minutes at the reduced settings used here (num.trees and n_perm
# are deliberately modest so the script is a smoke test, not a production run).
# For production-grade numbers raise NUM_TREES / N_PERM below.
#
# Output: prints a report and writes durum_verification_report.txt

suppressMessages(library(durumGEAr))
data(durumUnits)

# ---- reduced-but-adjustable settings ----------------------------------------
NUM_TREES <- 200L    # production: 600
N_PERM    <- 30L     # production: 100-200

con <- file("durum_verification_report.txt", open = "wt")
emit <- function(...) { s <- sprintf(...); cat(s); writeLines(s, con) }
pass <- character(0)
check <- function(label, ok) {
  pass[[label]] <<- if (isTRUE(ok)) "PASS" else "FAIL"
  emit("  %-45s %s\n", label, pass[[label]])
}

emit("durumGEAr v%s Verification Report\n", as.character(utils::packageVersion("durumGEAr")))
emit("==========================================\n\n")

## Test 1: dataset structure ---------------------------------------------------
emit("Test 1: Dataset structure\n")
n_units     <- nrow(durumUnits)
n_countries <- length(unique(durumUnits$Country))
n_targets   <- sum(grepl("^gPC[0-9]+$", names(durumUnits)))
emit("  N units      = %d\n", n_units)
emit("  N countries  = %d\n", n_countries)
emit("  N targets    = %d\n", n_targets)
check("dataset has >= 1000 units", n_units >= 1000)
check("dataset has 5 gPC targets", n_targets == 5L)

## Test 2: interpolation skill (spatialBlockCV) --------------------------------
emit("\nTest 2: Interpolation skill (spatialBlockCV)\n")
sb <- spatialBlockCV(durumUnits, seeds = 42, num.trees = NUM_TREES)
emit("  Mean R2 = %.4f\n", sb$metrics$mean_R2)
check("interpolation R2 is clearly positive (> 0.20)", sb$metrics$mean_R2 > 0.20)

## Test 3: extrapolation skill (locoCV) ----------------------------------------
emit("\nTest 3: Extrapolation skill (locoCV)\n")
lo <- locoCV(durumUnits, num.trees = NUM_TREES)
emit("  Mean R2 = %.4f\n", lo$metrics$mean_R2)
emit("  Countries beating baseline: %d / %d\n", lo$n_beat_baseline, lo$n_groups)
check("extrapolation R2 collapses toward/below zero (< 0.10)",
      lo$metrics$mean_R2 < 0.10)

## Test 4: confounding gap (confoundGapTest) -----------------------------------
emit("\nTest 4: Confounding gap (confoundGapTest)\n")
gt <- confoundGapTest(durumUnits, seed = 42, n_perm = N_PERM, num.trees = NUM_TREES)
emit("  Gap = %.4f  (interp R2 - LOCO R2)\n", gt$observed_gap)
emit("  p-value = %.4f\n", gt$p_value)
# The gap is large and positive (~0.3-0.5 depending on settings); at reduced
# num.trees it lands a little lower than the ~0.42 production value.
check("confounding gap is large and positive (> 0.25)", gt$observed_gap > 0.25)
check("gap p-value is significant (< 0.05)", gt$p_value < 0.05)

## Test 5: residual signal (residualize) ---------------------------------------
emit("\nTest 5: Residual within-country signal (residualize)\n")
rz <- residualize(durumUnits, seeds = 42:43, n_perm = N_PERM, num.trees = NUM_TREES)
emit("  %-6s %-28s %-10s %s\n", "target", "verdict", "fragility", "per_target_R2")
for (t in names(rz$verdict))
  emit("  %-6s %-28s %-10s %+.4f\n", t, rz$verdict[t],
       rz$verdict_fragility[t], rz$per_target_R2[t])
# gPC2 and gPC5 are the confirmed real-signal targets: positive multi-seed R2.
check("gPC2 shows positive within-country signal", rz$per_target_R2["gPC2"] > 0)
check("gPC5 shows positive within-country signal", rz$per_target_R2["gPC5"] > 0)
# gPC1 and gPC4 are country-identity artifacts. The package flags them EITHER by
# a 'country-identity artifact' verdict OR (when the single-seed statistic is
# positive but the multi-seed mean flips sign) by a 'fragile' fragility grade.
artifact_flagged <- function(t)
  rz$verdict[t] == "country-identity artifact" ||
  rz$verdict_fragility[t] == "fragile"
check("gPC1 flagged as artifact (verdict or fragility)", artifact_flagged("gPC1"))
check("gPC4 flagged as artifact (verdict or fragility)", artifact_flagged("gPC4"))

## Test 6: climate drivers (driverAnalysis) ------------------------------------
emit("\nTest 6: Climate drivers (driverAnalysis)\n")
dr2 <- driverAnalysis(durumUnits, "gPC2", seeds = 42, n_perm = N_PERM, num.trees = NUM_TREES)
dr5 <- driverAnalysis(durumUnits, "gPC5", seeds = 42, n_perm = N_PERM, num.trees = NUM_TREES)
top2 <- names(dr2$importance)[1]; top5 <- names(dr5$importance)[1]
emit("  gPC2 top driver (permutation importance) = %s\n", top2)
emit("  gPC5 top driver (permutation importance) = %s\n", top5)
# Permutation importance from the joint model is the headline driver metric;
# univariate residual R2 is noisy at reduced settings.
check("gPC2 top driver is BIO5", top2 == "BIO5")
check("gPC5 top driver is BIO15", top5 == "BIO15")

## Summary ---------------------------------------------------------------------
emit("\n==========================================\n")
n_pass <- sum(unlist(pass) == "PASS"); n_tot <- length(pass)
emit("Verification: %d / %d checks passed.\n", n_pass, n_tot)
if (n_pass == n_tot)
  emit("All checks passed - durumGEAr is correctly installed and functioning.\n") else
  emit("Some checks did not pass - see individual lines above. Note reduced\n  settings (NUM_TREES=%d, N_PERM=%d) can shift borderline numbers; rerun\n  with production settings before concluding a genuine failure.\n", NUM_TREES, N_PERM)
close(con)
invisible(NULL)
