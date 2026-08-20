# DES-SN5YR correction-bracket: is the model-comparison verdict correction-dominated?

Companion repository for *"Model discrimination in the DES-SN5YR Hubble diagram
is dominated by the simulation-derived bias-correction layer"* (draft in
`paper/`). Everything here runs from **public data** in **minutes**.

## The claim, in three numbers

Fitting flat ΛCDM (free Ω_m) and a one-parameter non-expanding probe law to the
public DES-SN5YR distance moduli (full STAT+SYS covariance, N = 1768):

| vector | Δχ² (probe − ΛCDM) |
|---|---|
| released `MU` (BBC-corrected) | **+1.6** (near-degenerate, ΛCDM ahead) |
| uncorrected `MU + biasCor_mu` | **+13.3** (ΛCDM ahead sevenfold more) |

The ΛCDM-simulation-derived bias-correction layer moves the verdict by ~12 χ²
units — seven times the verdict itself — *toward* the alternative. The
comparison is **correction-dominated**: at this precision the Hubble diagram
verdict is a property of the correction layer's simulation inputs as much as
of the sky. (Robustness: leverage 2–23× across covariance/subsample variants.)

## Reproduce it (under an hour, mostly download time)

```bash
pip install -r requirements.txt
cd data && ./fetch.sh && cd ..
python src/sign_check.py    # STEP ZERO: establishes biasCor sign empirically.
                            # (This check caught a sign error in our own
                            #  draft 1 — run it before anything else.)
python src/bracket.py       # the headline table
python src/robustness.py    # stat-only + anchors-excluded variants
```

## Falsify it

If a correct implementation (verified sign, full covariance C = C_sys +
diag(MUERR_FINAL²)) yields a swing smaller than the released-vector verdict,
the paper's claim is wrong. We invite exactly that check.

## Raw light curves (extended replication)

The raw griz forced photometry for all 8,293 DES transient candidates can be
rebuilt from the public 12 GB DES release with `tools/extract_raw_photometry.py`
(schema and sha256 audit path in its docstring), or fetched as a cached CSV
(Zenodo DOI pending). The bracket itself needs only the two small release files.

## Data & licenses

DES-SN5YR release files © the DES Collaboration, public, fetched at tag v1.2
(pinned — `main` has moved to the Dovekie recalibration). Code: MIT.
