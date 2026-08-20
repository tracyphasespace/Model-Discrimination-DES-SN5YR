#!/usr/bin/env python3
"""The sensitivity bracket (paper §5): fit two one-parameter distance laws to
the released DES-SN5YR moduli and to the uncorrected vector MU + biasCor_mu,
with the full STAT+SYS covariance and analytic marginalization of the
absolute offset. Reports all four chi^2, both Delta chi^2, and the swing.

Convention note: the uncorrected vector is MU + biasCor_mu — DES subtracts
the correction (verify with sign_check.py FIRST; getting this backwards
applies the correction twice and inverts the result).
"""
import gzip
import os
import numpy as np
import pandas as pd

DATA = os.environ.get("BRACKET_DATA", os.path.join(os.path.dirname(__file__), "..", "data"))


def load(all_entries=False):
    """Default sample excludes 61 effectively zero-weight entries
    (MUERR_FINAL > 10) -> N = 1768; this is OUR diagnostic cut, not the DES
    likelihood convention (which uses all 1829 with zHD > 0). Pass
    all_entries=True (or env BRACKET_ALL_1829=1) for the official-sample
    robustness variant."""
    hd_all = pd.read_csv(os.path.join(DATA, "DES-SN5YR_HD+MetaData.csv"))
    with gzip.open(os.path.join(DATA, "STAT+SYS.txt.gz"), "rt") as f:
        n = int(f.readline())
        c_sys = np.fromstring(f.read(), sep="\n").reshape(n, n)
    if n != len(hd_all):
        raise ValueError("covariance dimension mismatch")
    # statistical term is in MUERR_FINAL (release README); C_sys alone is not PD
    C_all = c_sys + np.diag(hd_all.MUERR_FINAL.to_numpy() ** 2)
    if all_entries or os.environ.get("BRACKET_ALL_1829", "0") == "1":
        sel = (hd_all.zHD > 0).to_numpy()
    else:
        sel = ((hd_all.zHD >= 0.01) & (hd_all.MUERR_FINAL > 0)
               & (hd_all.MUERR_FINAL <= 10)).to_numpy()
    hd = hd_all[sel].reset_index(drop=True)
    return hd, C_all[np.ix_(sel, sel)]


def make_chi2(C):
    Ci = np.linalg.inv(C)
    u = np.ones(len(C))
    uCu = u @ Ci @ u

    def chi2(r):
        r = r - (u @ Ci @ r) / uCu * u  # analytic offset marginalization
        return r @ Ci @ r
    return chi2


# Redshift convention (matches the official DES v1.2 likelihood): the
# cosmological integral uses zHD; the external luminosity-distance factor
# uses the HELIOCENTRIC redshift, mu ~ 5 log10[(1+zHEL) * D_M(zHD)].
# Model B adopts the analogous deliberate assignment: observed-frame flux
# factors (photon energy + arrival rate) -> (1+zHEL); the accumulated
# path/attenuation coordinate -> zHD.

def mu_modelB(zhd, zhel, eta):
    """Non-expanding probe: D = (c/K) ln(1+zHD); external factor (1+zHEL);
    attenuation (5/ln10) * eta * [1 - (1+zHD)^(-1/2)]. K absorbed in offset."""
    return (5 * np.log10((1 + zhel) * np.log1p(zhd))
            + (5 / np.log(10)) * eta * (1 - (1 + zhd) ** -0.5))


def mu_lcdm(zhd, zhel, om, n_grid=4000):
    zs = np.linspace(0, zhd.max(), n_grid)
    zm = 0.5 * (zs[1:] + zs[:-1])
    integ = np.concatenate([[0], np.cumsum(np.diff(zs) / np.sqrt(om * (1 + zm) ** 3 + 1 - om))])
    return 5 * np.log10((1 + zhel) * np.interp(zhd, zs, integ))


def fit(chi2, zhd, zhel, mu_obs):
    from scipy.optimize import minimize_scalar
    rB = minimize_scalar(lambda e: chi2(mu_obs - mu_modelB(zhd, zhel, e)),
                         bounds=(0.0, 2.0), method="bounded",
                         options={"xatol": 1e-4})
    rA = minimize_scalar(lambda o: chi2(mu_obs - mu_lcdm(zhd, zhel, o)),
                         bounds=(0.01, 1.2), method="bounded",
                         options={"xatol": 1e-4})
    return dict(eta=rB.x, chi2_B=rB.fun, om=rA.x, chi2_A=rA.fun)


def main():
    all_1829 = os.environ.get("BRACKET_ALL_1829", "0") == "1"
    hd, C = load(all_entries=all_1829)
    chi2 = make_chi2(C)
    zhd = hd.zHD.to_numpy()
    zhel = hd.zHEL.to_numpy()
    print(f"N = {len(hd)}, zHD in [{zhd.min():.3f}, {zhd.max():.3f}]"
          + ("  [ALL-1829 variant]" if all_1829 else "  [default N=1768 cut]"))

    b = hd.biasCor_mu.to_numpy()
    print(f"\nbiasCor_mu: median {np.median(b):+.4f}, range [{b.min():+.3f}, {b.max():+.3f}], "
          f"linear slope {np.polyfit(zhd, b, 1)[0]:+.4f} mag/z")

    results = {}
    for label, mu in [("released", hd.MU.to_numpy()),
                      ("uncorrected (MU + biasCor)", hd.MU.to_numpy() + b)]:
        r = fit(chi2, zhd, zhel, mu)
        d = r["chi2_B"] - r["chi2_A"]
        results[label] = d
        print(f"\n{label}:")
        print(f"  Model A (LCDM):  Om  = {r['om']:.3f}  chi2 = {r['chi2_A']:.1f}")
        print(f"  Model B (probe): eta = {r['eta']:.3f}  chi2 = {r['chi2_B']:.1f}")
        print(f"  Delta chi2 (B - A) = {d:+.1f}  ->  "
              f"{'A (LCDM) ahead' if d > 0 else 'B ahead'}")

    L = results["released"] - results["uncorrected (MU + biasCor)"]
    print(f"\nL_BBC = released - preBBC = {L:+.1f} chi2 units")
    print(f"|L_BBC| = {abs(L):.1f} "
          f"({abs(L / results['released']):.0f}x the released-vector verdict)")


if __name__ == "__main__":
    main()
