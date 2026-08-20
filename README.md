# Model Discrimination in DES-SN5YR

> **ADDENDUM (2026-08-20) — the headline number is release-specific.**
> This paper's L_BBC = −11.7 ("eleven times the verdict") is measured on
> DES-SN5YR **v1.2**, which DES has since superseded with the
> recalibrated **DES-Dovekie** release (Popovic et al. 2026). Our
> pre-registered rerun on Dovekie (`src/dovekie_rerun.py`, registration
> + result [here](https://github.com/tracyphasespace/Static-Universe-DES-SN5YR/blob/master/v8/DOVEKIE_RERUN_PREREG.md))
> finds verdict **+1.8**, pre-correction **+4.2**, **L_BBC = −2.4** — a
> five-fold reduction. The qualitative concern survives (leverage still
> exceeds the verdict, 1.3×); the quantitative headline does not
> transfer to the current release. Accurate framing: the concern was
> real for v1.2, and DES's own recalibration and SALT3 retraining
> materially reduced it — this paper is a diagnostic that the field's
> own improvements answered. Published number-pairs for the record:
> v1.2 → (+1.1, +12.8, −11.7); Dovekie → (+1.8, +4.2, −2.4). A Zenodo
> addendum version is pending (author's action).


**How sensitive is supernova-cosmology model discrimination to the
simulation-derived bias-correction layer?** Companion repository for
*"Model Discrimination in the DES-SN5YR Hubble Diagram Is Highly Sensitive
to the Simulation-Derived Bias-Correction Layer"* (T. McSheery; PDF
in `paper/`; published at DOI [10.5281/zenodo.22022089](https://doi.org/10.5281/zenodo.22022089)).

Everything here runs from **public DES data** on a laptop. The two release
files the headline needs total ~17 MB; the full raw light curves for all
8,293 DES transient candidates are included in `data/` (22 MB compressed).

## The finding, in one table

Two one-parameter distance laws — flat ΛCDM (free Ω_m) and a non-expanding
probe law (free η) — fitted to the public DES-SN5YR v1.2 moduli with the
full STAT+SYS covariance (N = 1768, offset marginalized):

| distance-modulus vector | Δχ² (probe − ΛCDM) |
|---|---|
| released `MU` (BBC-corrected) | **+1.1** (near-degenerate) |
| reconstructed pre-BBC (`MU + biasCor_mu`) | **+12.8** |

(Official DES redshift convention — `zHEL` in the external factor, `zHD` in
the cosmological integral; identical at quoted precision retaining all 1829
released entries.) The BBC step — computed from simulations generated under
a **reference** cosmology — carries a pairwise discrimination leverage of
**L_BBC = released − preBBC = −11.7 χ² units, roughly eleven times the
released verdict**, and moves this comparison *toward* the probe. The
sensitivity is not unique to this probe: flat wCDM (w = −0.5) shows
L = −5.7, and an Einstein–de Sitter control L = +365.3 (opposite sign — no
consistent partisanship). DES has directly tested the reference-cosmology
dependence (Camilleri et al. 2024, MNRAS 533, 2615) with simulation reruns
spanning ~0.15 mag of distance-modulus variation and found it small; this
probe's released-best-fit shape departs from ΛCDM by at most **0.055 mag**
offset-removed — well below that tested amplitude (`src/envelope.py`; an
amplitude comparison — Camilleri's formal criterion is a shape envelope,
and the pointwise containment test is the named upgrade). The untested
regime begins at larger departures (EdS: 0.276 mag) and with non-FLRW
functional coverage.

**What this does NOT show:** neither endpoint measures the universe
(stripping the correction restores real survey selection — χ²/dof degrades
0.93 → 1.34 — and the released covariance belongs to the corrected vector;
the reconstruction also retains the globally fitted α, β, γ, making it a
hybrid diagnostic). The bracket measures the *total* step leverage; no single
component is bounded by it. The quadrature sensitivity of the verdict to
twelve released BBC-related systematic modes is 0.71 χ² units
(`src/bbc_systematic_modes.py`) — below the margin — but none of those modes
is a reference-cosmology variation; that question is addressed by the
Camilleri simulation reruns and the envelope placement above.

## Replicate it

```bash
pip install -r requirements.txt
cd data && ./fetch.sh && cd ..     # pulls the two v1.2 release files

python src/sign_check.py           # STEP ZERO. Establishes the biasCor sign
                                   # empirically (Tripp reconstruction). This
                                   # check caught a sign error in our own
                                   # draft 1 — run it before anything else.
python src/bracket.py              # the headline table + layer characterization
python src/robustness.py           # stat-only covariance; anchors excluded
python src/projection.py           # rho_C = -0.068; L_fixed/adapt = -10.1/-1.6
python src/probes.py               # EdS and wCDM second probes; Tripp errors
python src/envelope.py             # Camilleri amplitude comparison
cd data && ./fetch_systematics.sh && cd ..   # 12 single-sys blocks, then:
python src/bbc_systematic_modes.py # quadrature sensitivity to those modes
```

Every number in the paper regenerates from these seven scripts. Expected
runtimes: seconds to ~3 minutes each.

## Falsify it

The central claim fails if a correct implementation — verified sign
convention (step zero), full covariance C = C_sys + diag(MUERR_FINAL²),
offset marginalized analytically — produces a correction-step swing that is
small relative to the released discrimination. We invite exactly that check.
The protocol has already demonstrated its teeth once, on the authors: draft 1
used the wrong sign and was caught by step zero before circulation.

## Data

| file | contents | provenance |
|---|---|---|
| `data/fetch.sh` | pulls `DES-SN5YR_HD+MetaData.csv` + `STAT+SYS.txt.gz` | DES-SN5YR public release, **pinned at tag v1.2** (`main` has moved to the Dovekie recalibration) |
| `data/des_sn5yr_raw_photometry.csv.gz` | raw griz forced photometry: 8,293 SNIDs, 769,647 rows (`snid, ra, dec, z, mjd, band, flux_nu_jy, flux_nu_jy_err, wavelength_eff_nm, snr, …`) | rebuilt from the DES public release; audit path + schema in `tools/extract_raw_photometry.py`; sha256 in `data/SHA256SUMS` |

Known conventions in the raw CSV: 26.5% of rows carry a sentinel state
(flux error ≤ 0 with snr = 0) — flag them, never divide by snr blindly.
`STAT+SYS.txt.gz` contains systematics ONLY despite the name; the
statistical term lives in `MUERR_FINAL` (release README).

DES-SN5YR data © the Dark Energy Survey Collaboration, used per their public
release. We thank the collaboration for a release detailed enough to make
independent tests like this one possible. Code: MIT.

## Scope note: DES-Dovekie supersession (2026-08-20)

DES has superseded the original SN5YR cosmology products with the
recalibrated **DES-Dovekie** analysis (Popovic et al. 2026,
arXiv:2511.07517; flat-ΛCDM Ωm = 0.330 ± 0.015 after cross-calibration
and SALT3 retraining; the upstream `des-science/DES-SN5YR` repository
HEAD now carries the Dovekie products). Every number in this paper is a
statement about the **frozen v1.2 release** — all fetches here are
pinned at tag v1.2, so the results remain exactly reproducible. A
pre-registered rerun of the model comparison on the Dovekie vector
(frozen estimator, reading bands committed before data inspection) has
been **executed** — `python3 src/dovekie_rerun.py --unblind` reproduces
it (v1.2 validation gate included): **Delta chi2 = +1.8** on N = 1820
(recovered Omega_m = 0.330 matches the published Dovekie value), i.e.
the near-degeneracy is robust to the recalibration per the
pre-committed reading bands. Secondary registered measurement: the
BBC-layer leverage on this comparison fell from **-11.7 (v1.2) to -2.4
(Dovekie)** — released and pre-correction margins are now same-sign and
both small, materially reducing this paper's correction-dominance
concern in the new release. Registration + result:
[companion repository](https://github.com/tracyphasespace/Static-Universe-DES-SN5YR/blob/master/v8/DOVEKIE_RERUN_PREREG.md).

## Companion paper

*A Static-Universe Two-Channel Photon Propagation Model Confronted with
DES-SN5YR* (McSheery 2026, Version 8.0, DOI
[10.5281/zenodo.22032094](https://doi.org/10.5281/zenodo.22032094);
prior records: V7 [10.5281/zenodo.22031627](https://doi.org/10.5281/zenodo.22031627),
Draft 6 [10.5281/zenodo.22031121](https://doi.org/10.5281/zenodo.22031121),
Draft 5 [10.5281/zenodo.22025329](https://doi.org/10.5281/zenodo.22025329)) presents
the propagation physics motivating this paper's probe law; its code lives in
[Static-Universe-DES-SN5YR](https://github.com/tracyphasespace/Static-Universe-DES-SN5YR).
The two works are logically independent — this repository's leverage
measurement makes no cosmological claim and stands regardless of that
interpretation.

## Paper

`paper/main.tex` is the canonical manuscript (Draft 4.1, LaTeX end-to-end;
PDF built by CI on every change). Draft 2.1 (markdown) is retained for
history only. Two disclosed self-corrections: an early draft's
"verdict reverses" claim was a sign-convention error caught by step zero;
draft 4's numbers used a simplified redshift convention, corrected to the
official zHD/zHEL treatment in 4.1 (headline moved +1.6 → +1.1, swing
unchanged at 11.7).

## Acknowledgments

The author would like to thank Anthropic, Google, and OpenAI for their
excellent tools which enabled the scripts and mathematical assistance as
well as collaborative/adversarial interactions to resolve the thousands
of details at the speed of electrons.
