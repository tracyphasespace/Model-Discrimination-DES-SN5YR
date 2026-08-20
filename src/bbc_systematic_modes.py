#!/usr/bin/env python3
"""Quadrature sensitivity of the pairwise verdict to twelve selected
BBC-related systematic modes (paper Section 6.1).

The release ships single-systematic covariance blocks for BBC-related
inputs (scatter model, dust/population realizations, host library,
selection efficiency, classifier swaps). NONE of these is a
reference-cosmology variation, so this script does NOT isolate the
cosmology-conditioned component; it measures the verdict movement carried
by the released systematic budget represented by these twelve modes.
Reference-cosmology dependence is addressed by the simulation reruns of
Camilleri et al. (2024) — see envelope.py.

Each block is (near) rank-1: C_k ~ d_k d_k^T. Perturb the released moduli
by +/- d_k, refit both models, record the verdict shift; quadrature over
blocks. Expected (v1.2): individual shifts to +/-0.4; quadrature ~0.7 chi2
units against the ~1-unit released margin.

Requires the SingleSYS_CovMatrix/ directory from the release (see data/fetch.sh).
"""
import gzip
import os
import numpy as np

from bracket import load, make_chi2, mu_modelB, mu_lcdm

DATA = os.environ.get("BRACKET_DATA", os.path.join(os.path.dirname(__file__), "..", "data"))
SS = os.path.join(DATA, "SingleSYS_CovMatrix")

BLOCKS = ["BS20", "P21SYS1", "P21SYS2", "P21SYS3", "P21_HOSTCOLOR", "W22_AGE",
          "HOSTEFFshift", "SVAHOSTLIB", "CClikelihood", "SNNtraining",
          "SNIRF", "SCONE"]


def main():
    hd, C = load()
    # rebuild the release-frame selection mask to subset the blocks
    import pandas as pd
    hd_all = pd.read_csv(os.path.join(DATA, "DES-SN5YR_HD+MetaData.csv"))
    sel = ((hd_all.zHD >= 0.01) & (hd_all.MUERR_FINAL > 0)
           & (hd_all.MUERR_FINAL <= 10)).to_numpy()

    chi2 = make_chi2(C)
    z = hd.zHD.to_numpy()
    zhel = hd.zHEL.to_numpy()

    from scipy.optimize import minimize_scalar
    def verdict(mu_obs):
        cq = minimize_scalar(lambda e: chi2(mu_obs - mu_modelB(z, zhel, e)),
                             bounds=(0.0, 2.0), method="bounded",
                             options={"xatol": 1e-3}).fun
        cl = minimize_scalar(lambda o: chi2(mu_obs - mu_lcdm(z, zhel, o)),
                             bounds=(0.01, 1.2), method="bounded",
                             options={"xatol": 1e-3}).fun
        return cq - cl

    base = verdict(hd.MU.to_numpy())
    print(f"baseline verdict (released MU): {base:+.2f}")
    print(f"{'block':>14s} {'rank1':>6s} {'shift+':>8s} {'shift-':>8s}")

    tot2 = 0.0
    for b in BLOCKS:
        with gzip.open(os.path.join(SS, b + ".txt.gz"), "rt") as f:
            n = int(f.readline())
            Ck = np.fromstring(f.read(), sep="\n").reshape(n, n)
        Ck = Ck[np.ix_(sel, sel)]
        w, V = np.linalg.eigh(Ck)
        d = V[:, -1] * np.sqrt(max(w[-1], 0.0))
        sp = verdict(hd.MU.to_numpy() + d) - base
        sm = verdict(hd.MU.to_numpy() - d) - base
        tot2 += max(abs(sp), abs(sm)) ** 2
        print(f"{b:>14s} {w[-1] / max(w.sum(), 1e-30):6.2f} {sp:+8.2f} {sm:+8.2f}")

    print(f"\nquadrature of max shifts: {np.sqrt(tot2):.2f} chi2 units "
          f"(released-vector margin: {base:+.2f})")
    print("Interpretation: the released DES systematic budget represented by")
    print("these twelve modes moves the verdict by less than the released")
    print("margin. This does NOT bound the reference-cosmology component;")
    print("see envelope.py / Camilleri et al. (2024) for that question.")


if __name__ == "__main__":
    main()
