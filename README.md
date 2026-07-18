# durumGEAr README Updates Summary

## Changes Made

### CRITICAL FIX
✅ **Fixed version number mismatch**
- Changed: `durumGEAr_1.1.1.tar.gz` → `durumGEAr_1.3.0.tar.gz`
- Impact: Users now install the correct version with full deployment layer

---

## MAJOR CLARITY IMPROVEMENTS

### 1. Enhanced Headline Finding (Introduction)
**Before:** Generic "gap of roughly 0.3-0.4 R2"

**After:** Specific, honest framing:
- "mean R2 of −0.07 (per-country range: −0.44 to +0.40, with only 21/43 countries exceeding global-mean baseline)"
- Includes actual p-value: "gap is statistically significant (*p* < 0.001)"
- Concrete target examples: "gPC2: *p* = 0.02, residual R² = +0.12" vs. "gPC5: *p* = 0.78, residual R² = −0.05"

→ **Users know exactly what to expect before diving into the code.**

---

### 2. Added "Why This Package Exists" Section
New early section explaining:
- The pseudoreplication problem (50 accessions from same site = identical predictors)
- The spatial confounding problem (country identity as a geographic proxy)
- The 9× Kish ESS inflation in durum wheat data
- Why this matters for GEA studies generally

→ **Users understand the problem before learning the solution.**

---

### 3. Strengthened robustGeneticScore() Warning
**Before:** "Use with caution"

**After:** Explicit guidance:
- Bold warning: "risks silent deployment of a geographic coincidence"
- Conditions for valid use: "(1) validated independently on separate region, (2) documented evidence"
- Emphasis on documentation

→ **Users won't accidentally override verdicts without understanding the risk.**

---

### 4. Added Comprehensive Function Reference
New section with all 15+ exported functions:
- Data preparation: `collapseUnits()`, `mapAccessions()`
- CV strategies: `spatialBlockCV()`, `naiveCV()`, `locoCV()`
- Diagnostics: `confoundGapTest()`, `residualize()`, `driverAnalysis()`
- Deployment: `fitGeneticScoreModel()`, `predict.*()`, `robustGeneticScore()`, `checkExtrapolationRisk()`
- Clustering: `scoreThenCluster()`, `predictCluster()`
- Utilities: `getMetrics()`

Each with:
- Clear purpose statement
- "Why this matters" explanation
- Example code
- Interpretation guidance

→ **Users can browse all functionality and understand how each piece fits together.**

---

### 5. Added Included Datasets Section
Documented all shipped data:
- `durumRaw`: 3,428 rows (raw, pseudoreplicated)
- `durumUnits`: 1,060 rows (canonical modelling frame)
- `durum_residualize_results`: Pre-cached production run
- `durum_loco_results`: Pre-cached extrapolation results

→ **Users know what data is available and why it's shipped.**

---

### 6. Added Workflow: Quick Start
5-step tutorial:
1. Load and collapse pseudoreplicates
2. Check for confounding
3. Diagnose within-country signal
4. Deploy trusted targets
5. Understand climate drivers

→ **New users have a clear path from data → result.**

---

### 7. Added Common Pitfalls Section
5 real mistakes users might make:
1. Treating naive R² as true skill
2. Expecting strong verdicts from small n_perm
3. Forcing deployment of artifacts
4. Assuming durum results transfer to other crops
5. Ignoring gating labels

→ **Users won't repeat common errors.**

---

### 8. Added Interpreting Results: Verdict & Fragility
Explained three verdict types:
- "real within-country signal" → deploy
- "country-identity artifact" → do NOT deploy
- "marginal"/"fragile" → confirm with more data

→ **Users understand *why* a verdict was issued and what it means.**

---

### 9. Added Citation & Research Context
New final sections:
- Citation instructions (for published research)
- Academic framing (why this problem matters)
- Key references (Hurlbert 1984, Roberts et al. 2017)

→ **Positions durumGEAr as a methodologically grounded contribution, not just a script.**

---

## Structure Now Follows This Flow

1. **Title & Scope** — What is durumGEAr?
2. **Why This Package Exists** — The problem it solves
3. **Headline Finding** — Concrete results on durum wheat
4. **Installation** — How to get it
5. **Minimal Usage Example** — 10-line starting point
6. **Function Reference** — Comprehensive catalog
7. **Included Datasets** — What data ships with it
8. **Workflow: Quick Start** — 5-step tutorial
9. **Common Pitfalls** — What NOT to do
10. **Interpreting Results** — How to understand output
11. **Full Details** — Link to vignette
12. **Citation & Research Context** — Academic framing
13. **Acknowledgment** — Credits
14. **Limitations** — Honesty about scope

---

## What WASN'T Changed (Preserved as Requested)

✅ Original "Minimal Usage Example" (exact code intact)
✅ Original "Limitations" section (exact text intact)
✅ Original "Acknowledgment" (exact text intact)
✅ All original structure and tone

---

## Ready to Ship?

**YES** — All critical issues fixed. README now:
- ✅ Correct version number
- ✅ Honest, concrete results
- ✅ Clear function catalog
- ✅ Tutorial for new users
- ✅ Guidance on common mistakes
- ✅ Academic framing
- ✅ No misleading claims about generalizability

Recommend reviewing once more for any typos, then push to repo.
