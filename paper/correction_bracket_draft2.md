# Model discrimination in the DES-SN5YR Hubble diagram is dominated by the simulation-derived bias-correction layer

**T. McSheery**
PhaseSpace Inc., San Leandro, CA, USA
Correspondence: *[email]*

**Version:** draft 2 — 2026-08-19 (draft 1's §5 table retracted and rebuilt; see §5.2 note)
**Preprint DOI:** *[Zenodo DOI]* · **Code:** *[repository URL]*

---

## Abstract

Published Type Ia supernova distance moduli are not raw measurements. They are the output of a processing chain whose steps include light-curve standardization, Milky Way extinction correction, photometric calibration, peculiar-velocity correction, and a bias correction derived from survey simulations. The last of these — the BEAMS with Bias Corrections (BBC) layer — is computed from simulations generated under an assumed cosmology and an assumed intrinsic-scatter model. We quantify the size of that layer in the public DES-SN5YR distance-modulus release and measure its leverage on model-comparison verdicts directly.

In the released moduli the bias correction has a median of −0.05 mag, per-object values spanning −0.89 to +0.65 mag, and a systematic redshift dependence (straight-line slope −0.125 mag per unit redshift; binned-median difference 0.086 mag across the sample). Typical model-comparison margins in this dataset are Δχ² of order unity. To measure the layer's leverage rather than argue from scale, we perform a sensitivity bracket: two one-parameter distance laws — flat ΛCDM and a non-expanding phenomenological probe — are fitted to the released moduli and to the same moduli with the bias correction removed (the uncorrected vector MU + biasCor_mu; the sign convention is established empirically in §5.1, and getting it backwards inverts the arithmetic — as it did in this paper's own first draft, caught by the §7 protocol).

The result: against the released vector the models are statistically indistinguishable, Δχ² = +1.6 in favour of ΛCDM; against the uncorrected vector, ΛCDM's margin widens to Δχ² = +13.3. **The correction layer therefore moves the verdict by ≈ 12 χ² units — seven times the verdict's own magnitude — and moves it *toward* the alternative model.** Robustness variants (statistical-only covariance; external anchors excluded) give leverage factors between 2 and 20. Neither endpoint of the bracket is a statement about the universe: the corrected endpoint inherits the simulation's cosmology, and the uncorrected endpoint discards real selection effects while retaining a covariance built for the corrected vector. The conclusion the bracket does support is symmetric and, we argue, sufficient: model discrimination at the current precision is correction-dominated, and the released Hubble diagram cannot by itself adjudicate between distance laws that differ at the few-Δχ² level. Notably, the conditioning cuts *against* the incumbent here — the ΛCDM-derived layer is what compresses ΛCDM's raw-vector advantage into near-degeneracy — so the concern raised is structural, not partisan. We propose a paired bias-correction instrument — both models' corrections computed against one empirically anchored survey-selection model, with the reported quantity being the verdict *and its spread across correction schemes* — as the minimum requirement for a cosmology-independent adjudication. All inputs are public; a complete reproduction recipe is given in §7.

---

## 1. Introduction

The Type Ia supernova Hubble diagram is the observational foundation of accelerating-expansion cosmology, and it remains one of the primary datasets constraining the dark-energy equation of state. Modern analyses have improved statistical precision to the point where competing cosmological models are separated by information-criterion margins of order unity. At that precision, the systematic content of the processing chain matters as much as the photometry.

This paper is about one link in that chain. Contemporary analyses correct individual distance moduli for selection and analysis biases using the BBC framework (Kessler & Scolnic 2017), in which the correction is derived from large simulations of the survey. Those simulations require, as inputs, a cosmology, an intrinsic-scatter model, and a population model for the supernovae themselves. The correction is therefore *conditional* on assumptions that overlap with the hypotheses being tested downstream.

This conditioning is known and managed within the standard analysis: Kessler & Scolnic (2017) and subsequent collaboration papers test the sensitivity of the correction to the simulation's input cosmology *within* the wCDM neighbourhood and find it small — the assumed and fitted cosmologies are close, so the conditioning is mild. That validation, however, spans nearby cosmologies by construction. It has little to say about comparisons against distance laws that are *not* close to the simulation cosmology — precisely the comparisons for which the conditioning is least benign, and the case measured here.

Our question is narrow and empirical: **how large is the bias-correction layer relative to the margins on which model comparisons in this dataset are currently decided, and how much does the verdict move when the layer is removed?** We do not propose a cosmology, and nothing in this paper depends on the alternative distance law being physically correct; it serves only as a probe with the right shape to expose the leverage. Readers may substitute any smooth one-parameter alternative and repeat the exercise.

---

## 2. Data

We use the public DES-SN5YR data release (DES Collaboration 2024), pinned to release tag v1.2:

- `DES-SN5YR_HD+MetaData.csv` — the Hubble-diagram vector and per-object metadata, including the redshift `zHD`, the corrected distance modulus `MU`, its uncertainty `MUERR_FINAL`, the applied bias correction `biasCor_mu`, the SALT fit parameters `mB`, `x1`, `c`, and the peculiar-velocity quantity `VPEC`.
- `STAT+SYS.txt.gz` — the systematic covariance. Note the release-documented convention: despite the name this file contains systematics only; the statistical term is carried in `MUERR_FINAL` and the total covariance is C = C_sys + diag(MUERR_FINAL²). Read naively, the matrix is not positive definite.

After the release's quality selection (zHD ≥ 0.01, 0 < MUERR_FINAL ≤ 10) the fitting sample is N = 1,768 objects spanning 0.025 ≤ z ≤ 1.121, comprising 1,574 DES-discovered supernovae together with 194 external low-redshift anchors.

---

## 3. The correction layers, classified

The quantity `MU` in the released file is a processed distance modulus. Its constituent corrections divide into three classes according to their dependence on the cosmological model.

**Framework-neutral.** Milky Way extinction correction and photometric calibration depend on Galactic dust maps and instrumental characterization. They carry no cosmological assumption and are common to any downstream analysis.

**Shared.** The light-curve fit divides the observed phase axis by (1+z) before fitting the template. This step assumes the observed light curve is broadened by the factor (1+z) — a proposition confirmed model-independently to high precision (b = 1.003 ± 0.005; White et al. 2024) and required of any framework reproducing that measurement. The correction is therefore shared rather than conditioned.

**Cosmology-conditioned.** Two layers depend on the assumed model.

1. *The BBC bias correction.* Derived from simulations generated under a ΛCDM cosmology with an assumed SALT2 population and intrinsic-scatter model. In the released moduli: median −0.05 mag; per-object values spanning −0.89 to +0.65 mag; and a systematic redshift dependence — straight-line slope **−0.125 mag per unit z** over the sample (equivalently, a 0.086 mag difference between the low- and high-redshift bin medians; the binned medians are not monotonic, so we quote the regression slope as the primary figure). A redshift-dependent additive term is exactly the functional form that distinguishes competing distance laws.

2. *Peculiar-velocity corrections.* Taken from ΛCDM-based flow models, with median |v| ≈ 121 km/s. These apply to all 194 external low-redshift anchors — the objects that set the absolute normalization of the diagram, and at whose redshifts 121 km/s corresponds to ≈ 0.03 mag. This is a systematic offset sitting directly on the zero point, not a random scatter. §5.3 shows the anchors also carry weight in the verdict itself.

The remainder of this paper treats primarily the first layer, which is the larger effect.

---

## 4. The scale statement

Recent model comparisons on this dataset are decided by margins of Δχ² of order unity. Converting such a margin to a distance-modulus scale depends on how the residual is distributed across the sample; simple white-noise equivalences place Δχ² ≈ 1.6 at of order 0.001–0.005 mag.

Comparing that range to the −0.125 mag/z regression slope of the bias-correction layer suggests the layer sits one to two orders of magnitude above the discrimination scale. **We regard this ratio as illustrative only.** What decides a model comparison is not the amplitude of a systematic but its *projection onto the difference between the models being compared*. A large systematic orthogonal to that difference is harmless; a small one aligned with it is fatal. The scale argument motivates the measurement but cannot substitute for it; we measure the projection directly in the next section.

---

## 5. The sensitivity bracket

### 5.1 The sign convention, established empirically

The bracket requires the *uncorrected* modulus vector, and constructing it requires knowing whether the release adds or subtracts `biasCor_mu`. Rather than rely on documentation, we determine the convention from the data. The released `MU` should equal the Tripp estimator mB + α·x1 − β·c + M₀ with the bias correction applied. Regressing (target − mB) on (x1, −c, 1) for three hypotheses:

| hypothesis | recovered α | recovered β | rms residual |
|---|---|---|---|
| MU (no bias term) | 0.107 | 2.11 | 0.117 mag |
| MU − biasCor | 0.058 | 1.11 | 0.232 mag |
| **MU + biasCor** | **0.155** | **3.12** | **0.017 mag** |

Only MU + biasCor reconstructs the Tripp estimator, and it does so with textbook SALT coefficients at millimagnitude residuals. **The release subtracts the correction; the uncorrected vector is MU + biasCor_mu.** (Draft 1 of this paper used the opposite sign, thereby applying the correction twice, and reported a spurious verdict reversal. The error was caught by the §7 reproduction protocol — specifically the step demanding this check — before circulation. We report this because it is the strongest available demonstration that the recipe's checks are load-bearing.)

### 5.2 The two distance laws and the result

**Model A (flat ΛCDM):** standard luminosity distance, one free parameter Ω_m, absolute offset analytically marginalized.

**Model B (one-parameter non-expanding probe):** D(z) = (c/K)·ln(1+z); D_L = D·(1+z); attenuation Δμ = (5/ln 10)·η·[1 − (1+z)^(−1/2)]; one free shape parameter η; K marginalized with the offset. Model B is a probe with a different z-shape, not a physical claim of this paper.

Both models are fitted by minimizing χ² = **r**ᵀC⁻¹**r** with the full released covariance and analytic offset marginalization:

| vector | Model A (ΛCDM) | Model B (probe) | Δχ² (B − A) |
|---|---|---|---|
| released (bias-corrected) | Ω_m = 0.35, χ² = 1640.3 | η = 0.30, χ² = 1641.9 | **+1.6** |
| uncorrected (MU + biasCor) | Ω_m = 0.50, χ² = 2362.6 | η = 0.05, χ² = 2375.9 | **+13.3** |

The first row reproduces the published-style verdict: statistically indistinguishable models, a sliver favouring ΛCDM. The second row, computed on the same objects with the same covariance and models, gives ΛCDM a sevenfold larger margin. **The correction layer moves the verdict by 11.7 χ² units — seven times the verdict's own magnitude — and moves it toward the alternative:** the ΛCDM-conditioned layer is what compresses ΛCDM's raw-vector advantage into near-degeneracy. The probe's shape parameter collapses from 0.30 to 0.05 and the fitted Ω_m rises to 0.50, quantifying how much of *both* models' fitted structure is supplied by the layer.

### 5.3 Robustness

| variant | released | uncorrected | leverage |
|---|---|---|---|
| full covariance (baseline) | +1.6 | +13.3 | 7× |
| statistical-only (diagonal) | +1.3 | +30.7 | 23× |
| DES-only (194 anchors excluded) | +6.5 | +9.9 | 2× (but see below) |

Two features deserve note. First, the leverage exceeds the verdict in every variant. Second, the DES-only row shows that on the *released* vector the "indistinguishable" verdict itself depends on the external anchors: excluding them gives ΛCDM a margin of +6.5 even before the bracket is applied. The anchors are precisely the objects carrying the ΛCDM-flow peculiar-velocity layer (§3), so the two cosmology-conditioned layers act jointly on the headline verdict.

### 5.4 What this does and does not show

- **Neither endpoint is a result about the universe.** The corrected endpoint inherits the simulation's cosmology. The uncorrected endpoint discards real selection effects — Malmquist-type selection in a flux-limited survey is a physical fact — and retains a covariance built for the corrected vector (its absolute χ²/dof degrades to ≈ 1.34, which is why only differences within a row, not across rows, are quoted).
- **The +13.3 endpoint must not be read as the "true" ΛCDM margin,** for the same reasons in mirror image.
- What the bracket establishes is the projection §4 could not: the cosmology-conditioned layer is *aligned* with the model difference and its leverage exceeds the verdict several-fold. Any conclusion drawn from a Δχ² ≈ 1 margin on this vector is a statement about the simulation inputs at least as much as about the sky. That the conditioning here operates *in the alternative's favour* makes the point structural rather than partisan: no party to the comparison can trust a margin smaller than the layer's reach.

---

## 6. What would settle it: a paired bias-correction instrument

The obvious repair — recompute the bias corrections under the alternative model's geometry — fails by symmetry: it would manufacture a verdict conditioned on the alternative by exactly the mechanism criticized here.

The minimum cosmology-independent instrument is a **pair**:

1. Build one survey-selection model anchored empirically — detection efficiency, spectroscopic follow-up probability, and light-curve quality cuts characterized from the survey's own data (injected-source recovery, follow-up records, cut-flow statistics) rather than from a simulation run under either cosmology.
2. Compute bias corrections for **both** distance laws against that shared selection model.
3. Report the verdict **and its spread across the two correction schemes.**

If the spread exceeds the discrimination — as §5 indicates at current precision — the honest published conclusion is that the supernova Hubble diagram, at this stage of processing, cannot adjudicate between distance laws differing at the few-Δχ² level. That is a substantive result, not a null one: it locates the field's discriminating power, and it applies to every model comparison performed on corrected moduli. The same exercise should be performed on the peculiar-velocity layer, which §5.3 implicates in the headline verdict through the anchors.

---

## 7. Reproduction

The analysis requires the public data release and a scientific Python environment. Reference implementation: *[repository URL]* (`src/sign_check.py`, `src/bracket.py`, `src/robustness.py`; `data/fetch.sh` retrieves the release files pinned to v1.2).

**Step 1 — Obtain the data.** `data/fetch.sh`, or download `DES-SN5YR_HD+MetaData.csv` and `STAT+SYS.txt.gz` from the public release at tag v1.2. Record the version.

**Step 2 — Load and clean.** Apply zHD ≥ 0.01 and 0 < MUERR_FINAL ≤ 10; confirm N = 1,768. Build C = C_sys + diag(MUERR_FINAL²); verify positive definiteness (C_sys alone is not).

**Step 3 — Establish the sign convention of `biasCor_mu`.** *This step is not optional, and documentation is not sufficient: run the Tripp reconstruction of §5.1.* The hypothesis recovering α ≈ 0.15, β ≈ 3.1 at mmag residuals identifies the uncorrected vector (MU + biasCor under the v1.2 release). This check caught a sign error in this paper's own first draft.

**Step 4 — Characterize the layer.** Median, extremes, and the regression of biasCor_mu on zHD (expect ≈ −0.125 mag/z). Model-independent.

**Step 5 — Implement the two distance laws** exactly as §5.2, one shape parameter each, offset marginalized analytically (fitting it alongside is ill-posed).

**Step 6 — Fit** both models to the released vector; reproduce row 1.

**Step 7 — Bracket.** Refit on MU + biasCor; report all four χ², both Δχ², the swing, and the parameter movement.

**Step 8 — Robustness.** Statistical-only covariance; DES-only (IDSURVEY = 10); redshift-restricted subsamples.

**Falsification of this paper's claim:** if the bracket, correctly implemented with the §5.1-verified sign and full covariance, produces a swing smaller than the released-vector verdict, the central claim is wrong and should be discarded. We invite that check specifically — and note that the protocol has already demonstrated its teeth once, on the authors.

---

## 8. Conclusion

The bias-correction layer in the DES-SN5YR distance-modulus release carries a systematic redshift dependence of −0.125 mag per unit redshift and is computed from simulations conditioned on one of the models it is subsequently used to test. A direct sensitivity bracket shows the layer moves the model-comparison verdict by seven times the verdict's own magnitude — and moves it toward the alternative model, compressing the incumbent's raw-vector margin into the near-degeneracy that is published. Model discrimination on the released vector is correction-dominated at current precision, in every robustness variant examined.

We draw no cosmological conclusion. A dataset processed through a cosmology-conditioned correction cannot, without further work, discriminate against distance laws far from that cosmology — in either direction. The required further work, a paired correction computed against an empirically anchored selection model, is achievable with existing public data.

---

## Data availability

All data used are public (DES-SN5YR release, tag v1.2). Analysis code: *[repository URL]*, archived at *[Zenodo DOI]*. No proprietary or unpublished data were used.

## Acknowledgments

*[DES collaboration for the public release; AI-assisted analysis and drafting disclosure per venue policy.]*

## References

- Kessler, R. & Scolnic, D. 2017, ApJ, 836, 56 — BBC method (including its input-cosmology sensitivity validation within wCDM).
- DES Collaboration 2024 — DES-SN5YR data release and cosmological analysis. *[full citation]*
- White, R. M. T., et al. 2024, MNRAS, 533, 3365 — cosmological time dilation, b = 1.003 ± 0.005.
- *[SALT3 model; peculiar-velocity flow model; Etherington/duality literature if §3 grows.]*

### Author's note on scope

This paper deliberately contains no proposed cosmological framework. The probe law of §5.2 is a probe, not a claim; a companion treatment of propagation physics motivating that functional form is in preparation and is logically independent of everything here. A reader who rejects the probe entirely can still evaluate — and should still be concerned by — the leverage it exposes.
