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


def load():
    hd_all = pd.read_csv(os.path.join(DATA, "DES-SN5YR_HD+MetaData.csv"))
    with gzip.open(os.path.join(DATA, "STAT+SYS.txt.gz"), "rt") as f:
        n = int(f.readline())
        c_sys = np.fromstring(f.read(), sep="\n").reshape(n, n)
    if n != len(hd_all):
        raise ValueError("covariance dimension mismatch")
    # statistical term is in MUERR_FINAL (release README); C_sys alone is not PD
    C_all = c_sys + np.diag(hd_all.MUERR_FINAL.to_numpy() ** 2)
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


def mu_modelB(z, eta):
    """Non-expanding probe law: D = (c/K) ln(1+z), D_L = D (1+z), plus
    attenuation (5/ln10) * eta * [1 - (1+z)^(-1/2)]. K absorbed in offset."""
    return 5 * np.log10((1 + z) * np.log1p(z)) + (5 / np.log(10)) * eta * (1 - (1 + z) ** -0.5)


def mu_lcdm(z, om, n_grid=400):
    zs = np.linspace(0, z.max(), n_grid)
    zm = 0.5 * (zs[1:] + zs[:-1])
    integ = np.concatenate([[0], np.cumsum(np.diff(zs) / np.sqrt(om * (1 + zm) ** 3 + 1 - om))])
    return 5 * np.log10((1 + z) * np.interp(z, zs, integ))


def fit(chi2, z, mu_obs):
    etas = np.arange(0.0, 1.5001, 0.01)
    cq = [chi2(mu_obs - mu_modelB(z, e)) for e in etas]
    iq = int(np.argmin(cq))
    oms = np.arange(0.02, 0.7001, 0.01)
    cl = [chi2(mu_obs - mu_lcdm(z, o)) for o in oms]
    il = int(np.argmin(cl))
    return dict(eta=etas[iq], chi2_B=cq[iq], om=oms[il], chi2_A=cl[il])


def main():
    hd, C = load()
    chi2 = make_chi2(C)
    z = hd.zHD.to_numpy()
    print(f"N = {len(hd)}, z in [{z.min():.3f}, {z.max():.3f}]")

    b = hd.biasCor_mu.to_numpy()
    print(f"\nbiasCor_mu: median {np.median(b):+.4f}, range [{b.min():+.3f}, {b.max():+.3f}], "
          f"linear slope {np.polyfit(z, b, 1)[0]:+.4f} mag/z")

    results = {}
    for label, mu in [("released", hd.MU.to_numpy()),
                      ("uncorrected (MU + biasCor)", hd.MU.to_numpy() + b)]:
        r = fit(chi2, z, mu)
        d = r["chi2_B"] - r["chi2_A"]
        results[label] = d
        print(f"\n{label}:")
        print(f"  Model A (LCDM):  Om  = {r['om']:.2f}  chi2 = {r['chi2_A']:.1f}")
        print(f"  Model B (probe): eta = {r['eta']:.2f}  chi2 = {r['chi2_B']:.1f}")
        print(f"  Delta chi2 (B - A) = {d:+.1f}  ->  "
              f"{'A (LCDM) ahead' if d > 0 else 'B ahead'}")

    swing = results["uncorrected (MU + biasCor)"] - results["released"]
    print(f"\nSWING = {swing:+.1f} chi2 units "
          f"({abs(swing / results['released']):.0f}x the released-vector verdict)")
    print("The verdict is correction-dominated if |swing| exceeds |verdict|.")


if __name__ == "__main__":
    main()
