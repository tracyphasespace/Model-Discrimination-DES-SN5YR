#!/usr/bin/env python3
"""Robustness variants (paper §7 Step 8): statistical-only covariance and
DES-only (194 external anchors excluded). Reuses bracket.py machinery."""
import numpy as np

from bracket import load, make_chi2, fit


def report(tag, hd, C):
    chi2 = make_chi2(C)
    zhd, zhel = hd.zHD.to_numpy(), hd.zHEL.to_numpy()
    b = hd.biasCor_mu.to_numpy()
    out = []
    for label, mu in [("released", hd.MU.to_numpy()),
                      ("uncorrected", hd.MU.to_numpy() + b)]:
        r = fit(chi2, zhd, zhel, mu)
        out.append(f"{label}: D={r['chi2_B'] - r['chi2_A']:+.1f} "
                   f"(eta={r['eta']:.2f}, Om={r['om']:.2f})")
    print(f"{tag:32s} " + " | ".join(out))


def main():
    hd, C = load()
    report("full covariance (baseline)", hd, C)
    report("stat-only (diagonal MUERR)", hd, np.diag(hd.MUERR_FINAL.to_numpy() ** 2))
    des = (hd.IDSURVEY == 10).to_numpy()
    report("DES-only (anchors excluded)", hd[des].reset_index(drop=True),
           C[np.ix_(des, des)])


if __name__ == "__main__":
    main()
