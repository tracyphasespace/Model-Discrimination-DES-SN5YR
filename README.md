# Model Discrimination in DES-SN5YR

**How sensitive is supernova-cosmology model discrimination to the
simulation-derived bias-correction layer?** Companion repository for
*"Model Discrimination in the DES-SN5YR Hubble Diagram Is Highly Sensitive
to the Simulation-Derived Bias-Correction Layer"* (T. McSheery; PDF will be
added to `paper/` on release; reserved archive DOI 10.5281/zenodo.22022089).

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

## Paper

`paper/main.tex` is the canonical manuscript (Draft 4.1, LaTeX end-to-end;
PDF built by CI on every change). Draft 2.1 (markdown) is retained for
history only. Two disclosed self-corrections: an early draft's
"verdict reverses" claim was a sign-convention error caught by step zero;
draft 4's numbers used a simplified redshift convention, corrected to the
official zHD/zHEL treatment in 4.1 (headline moved +1.6 → +1.1, swing
unchanged at 11.7).
