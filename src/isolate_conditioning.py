#!/usr/bin/env python3
"""Isolate the conditioning component within DES's own error budget (paper §5.5).

The release ships single-systematic covariance blocks for the bias-correction
layer's simulation inputs (scatter model, dust/population realizations, host
library, selection efficiency, classifier swaps). Each is (near) rank-1:
C_k ~ d_k d_k^T where d_k is the modulus shift from swapping that input.
Perturb the released moduli by +/- d_k, refit both models, and record the
verdict shift. The quadrature over blocks measures the conditioning-INPUT
sensitivity as DES themselves budget it — i.e. within the LCDM neighborhood.

Expected (v1.2): individual shifts +/-0.01 to +/-0.40; quadrature ~0.7 chi2
units against the 1.6 margin. DES's own budget therefore certifies the
verdict within the neighborhood; the unvalidated territory is the remainder
between this 0.7 and the bracket's ~12-unit total leverage.

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

    def verdict(mu_obs):
        cq = min(chi2(mu_obs - mu_modelB(z, e)) for e in np.arange(0, 1.5, 0.01))
        cl = min(chi2(mu_obs - mu_lcdm(z, o)) for o in np.arange(0.02, 0.7, 0.01))
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
    print("Interpretation: this measures input sensitivity WITHIN the LCDM")
    print("neighborhood — DES's own budget. The unvalidated conditioning against")
    print("distant distance laws lies between this figure and the bracket's")
    print("total-leverage bound (bracket.py).")


if __name__ == "__main__":
    main()
