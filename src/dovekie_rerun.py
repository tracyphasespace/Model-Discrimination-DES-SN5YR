#!/usr/bin/env python3
"""Preregistered static-law rerun on the DES-Dovekie recalibrated vector.

Registration (committed before any Dovekie data values were inspected):
Static-Universe-DES-SN5YR repo, v8/DOVEKIE_RERUN_PREREG.md. Estimator is
FROZEN: the same two one-parameter laws, same optimizer bounds/tolerance,
same zHD/zHEL convention, same analytic offset marginalization as
`bracket.py` (whose v1.2 result Delta chi2 = +1.1 is CI-verified).

Data: data/dovekie/{DES-Dovekie_HD.csv, STAT+SYS.npz} fetched from
des-science/DES-SN5YR (upstream HEAD; supersedes v1.2). Per the official
DES-Dovekie-SN_Likelihood.py, the npz stores the INVERSE of the full
stat+sys covariance, upper-triangle packed, N = 1820.

Blindness protocol: `--validate` exercises the shared code path on the
frozen v1.2 vector and must reproduce the manifest values. The Dovekie
measurement itself runs only with `--unblind` (the author's call).
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from bracket import load, make_chi2, fit, mu_modelB, mu_lcdm  # noqa: E402

DDIR = os.path.join(os.path.dirname(__file__), "..", "data", "dovekie")


def load_dovekie():
    hd = pd.read_csv(os.path.join(DDIR, "DES-Dovekie_HD.csv"),
                     sep=r"\s+", comment="#")
    hd = hd.drop(columns=["VARNAMES:"])
    hd.columns = ["CID", "IDSURVEY", "zHD", "zHEL", "MU", "MUERR",
                  "MUERR_VPEC", "MUERR_SYS", "PROBIA_BEAMS"]
    d = np.load(os.path.join(DDIR, "STAT+SYS.npz"))
    n = int(np.atleast_1d(d["nsn"])[0])
    assert n == len(hd), f"covariance n={n} != HD rows {len(hd)}"
    inv_cov = np.zeros((n, n))
    inv_cov[np.triu_indices(n)] = d["cov"]
    lo = np.tril_indices(n, -1)
    inv_cov[lo] = inv_cov.T[lo]
    return hd, inv_cov


def make_chi2_from_inv(Ci):
    """Offset-marginalized chi2 built from an INVERSE covariance
    (same algebra as bracket.make_chi2, which starts from C)."""
    u = np.ones(len(Ci))
    uCu = u @ Ci @ u

    def chi2(r):
        r = r - (u @ Ci @ r) / uCu * u
        return r @ Ci @ r
    return chi2


def validate():
    hd, C = load()
    chi2 = make_chi2(C)
    r = fit(chi2, hd.zHD.to_numpy(), hd.zHEL.to_numpy(), hd.MU.to_numpy())
    m = r["chi2_B"] - r["chi2_A"]
    print(f"[validate v1.2] N={len(hd)}  chi2_LCDM={r['chi2_A']:.1f}  "
          f"chi2_static={r['chi2_B']:.1f}  Delta chi2={m:+.1f}  "
          f"(manifest: +1.1)  Om={r['om']:.3f}  eta={r['eta']:.3f}")
    ok = abs(m - 1.1) < 0.15
    print("VALIDATION", "PASS" if ok else "FAIL")
    return ok


def run_dovekie():
    hd, Ci = load_dovekie()
    chi2 = make_chi2_from_inv(Ci)
    r = fit(chi2, hd.zHD.to_numpy(), hd.zHEL.to_numpy(), hd.MU.to_numpy())
    m = r["chi2_B"] - r["chi2_A"]
    print(f"[DOVEKIE, preregistered] N={len(hd)}")
    print(f"  chi2_LCDM   = {r['chi2_A']:.1f}   (Om    = {r['om']:.3f})")
    print(f"  chi2_static = {r['chi2_B']:.1f}   (eta   = {r['eta']:.3f})")
    print(f"  Delta chi2  = {m:+.1f}")
    print("Preregistered reading: |Dchi2|<~5 robust near-degeneracy; "
          ">~+25 adverse for the static law; strongly negative -> caution "
          "+ systematics search; intermediate -> report as measured.")


if __name__ == "__main__":
    if "--unblind" in sys.argv:
        if not validate():
            sys.exit("v1.2 validation failed — Dovekie run aborted.")
        run_dovekie()
    else:
        validate()
        print("\nDovekie measurement is GATED: rerun with --unblind "
              "(author's call; registration v8/DOVEKIE_RERUN_PREREG.md).")
