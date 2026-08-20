#!/usr/bin/env python3
"""Manifest-to-computation verifier: recompute the headline and diagnostic
results and assert them against frozen_results.json with numeric tolerances.

This closes the loop the reviewer identified: scripts -> computed results
-> frozen_results.json <- README/main.tex (the prose arrows are checked by
check_sync.py). The CI workflow contains NO scientific literals; this file
and the manifest are the only sources of expected truth, and this file
derives its expectations FROM the manifest.
"""
import gzip
import json
import os
import sys

import numpy as np

from bracket import load, make_chi2, mu_modelB, mu_lcdm, fit
from probes import mu_wcdm

ROOT = os.path.join(os.path.dirname(__file__), "..")
F = json.load(open(os.path.join(ROOT, "frozen_results.json")))


def num(key):
    return float(F[key].replace("−", "-"))


FAIL = 0


def check(name, computed, key, tol):
    global FAIL
    expected = num(key)
    ok = abs(computed - expected) <= tol
    print(f"{'OK ' if ok else 'FAIL'} {name}: computed {computed:+.4g}, "
          f"frozen {expected:+.4g}, tol {tol}")
    if not ok:
        FAIL += 1


def main():
    hd, C = load()
    chi2 = make_chi2(C)
    Ci = np.linalg.inv(C)
    u = np.ones(len(hd))
    uCu = u @ Ci @ u
    zhd, zhel = hd.zHD.to_numpy(), hd.zHEL.to_numpy()
    MU = hd.MU.to_numpy()
    b = hd.biasCor_mu.to_numpy()
    pre = MU + b

    def proj(x):
        return x - (u @ Ci @ x) / uCu * u

    # --- bracket ---
    r_rel = fit(chi2, zhd, zhel, MU)
    r_pre = fit(chi2, zhd, zhel, pre)
    d_rel = r_rel["chi2_B"] - r_rel["chi2_A"]
    d_pre = r_pre["chi2_B"] - r_pre["chi2_A"]
    check("Delta released", d_rel, "delta_chi2_released", 0.05)
    check("Delta preBBC", d_pre, "delta_chi2_prebbc", 0.05)
    check("L_BBC", d_rel - d_pre, "L_BBC", 0.05)
    check("Om released", r_rel["om"], "om_released", 0.002)
    check("eta released", r_rel["eta"], "eta_released", 0.002)
    check("Om preBBC", r_pre["om"], "om_prebbc", 0.002)
    check("eta preBBC", r_pre["eta"], "eta_prebbc", 0.002)

    # --- all-1829 variant: L identical at tolerance ---
    hd2, C2 = load(all_entries=True)
    chi2_2 = make_chi2(C2)
    z2, zh2 = hd2.zHD.to_numpy(), hd2.zHEL.to_numpy()
    b2 = hd2.biasCor_mu.to_numpy()
    ra = fit(chi2_2, z2, zh2, hd2.MU.to_numpy())
    rb = fit(chi2_2, z2, zh2, hd2.MU.to_numpy() + b2)
    L1829 = (ra["chi2_B"] - ra["chi2_A"]) - (rb["chi2_B"] - rb["chi2_A"])
    check("L_BBC (all-1829)", L1829, "L_BBC", 0.1)

    # --- biasCor characterization ---
    check("biasCor median", float(np.median(b)), "biascor_median", 0.005)
    check("biasCor slope", float(np.polyfit(zhd, b, 1)[0]), "biascor_slope", 0.005)

    # --- projection & decomposition ---
    # d is defined at the PUBLISHED (rounded) released-vector fits, per the
    # manuscript's Section 5.8; the manifest is the single source of those.
    om0, eta0 = num("om_released"), num("eta_released")
    mA, mB = mu_lcdm(zhd, zhel, om0), mu_modelB(zhd, zhel, eta0)
    bt, dt = proj(b), proj(mB - mA)
    nb, nd = np.sqrt(bt @ Ci @ bt), np.sqrt(dt @ Ci @ dt)
    rho = (bt @ Ci @ dt) / (nb * nd)
    check("rho_C", rho, "rho_C", 0.003)
    check("|b|", nb, "b_norm", 0.05)
    check("|d|", nd, "d_norm", 0.05)

    def q(r):
        r = proj(r)
        return r @ Ci @ r
    L_fixed = (q(MU - mB) - q(MU - mA)) - (q(pre - mB) - q(pre - mA))
    check("L_fixed", L_fixed, "L_fixed", 0.1)
    check("L_adapt", (d_rel - d_pre) - L_fixed, "L_adapt", 0.15)

    mA2 = mu_lcdm(zhd, zhel, num("om_prebbc"))
    mB2 = mu_modelB(zhd, zhel, num("eta_prebbc"))
    dt2 = proj(mB2 - mA2)
    rho2 = (bt @ Ci @ dt2) / (nb * np.sqrt(dt2 @ Ci @ dt2))
    check("rho_C (pre-BBC fits)", rho2, "rho_C_alt", 0.003)

    # --- second probes ---
    from scipy.optimize import minimize_scalar

    def fitA(mu_obs):
        return minimize_scalar(lambda o: chi2(mu_obs - mu_lcdm(zhd, zhel, o)),
                               bounds=(0.01, 1.2), method="bounded",
                               options={"xatol": 1e-4}).fun

    def d_eds(mu_obs):
        return chi2(mu_obs - mu_lcdm(zhd, zhel, 1.0)) - fitA(mu_obs)

    def d_wcdm(mu_obs):
        c = minimize_scalar(lambda o: chi2(mu_obs - mu_wcdm(zhd, zhel, o)),
                            bounds=(0.01, 1.0), method="bounded",
                            options={"xatol": 1e-4}).fun
        return c - fitA(mu_obs)

    check("L_EdS", d_eds(MU) - d_eds(pre), "L_eds", 0.5)
    check("L_wCDM", d_wcdm(MU) - d_wcdm(pre), "L_wcdm", 0.1)

    # --- amplitude comparison ---
    probe_amp = float(np.abs(proj(mu_modelB(zhd, zhel, 0.297)
                                  - mu_lcdm(zhd, zhel, 0.352))).max())
    eds_amp = float(np.abs(proj(mu_lcdm(zhd, zhel, 1.0)
                                - mu_lcdm(zhd, zhel, 0.352))).max())
    check("probe amplitude", probe_amp, "envelope_probe_mag", 0.003)
    check("EdS amplitude", eds_amp, "envelope_eds_mag", 0.005)

    # --- Tripp regression ---
    y = (MU + b - hd.mB.to_numpy())
    X = np.column_stack([hd.x1.to_numpy(), -hd.c.to_numpy(), np.ones(len(hd))])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    check("Tripp alpha", coef[0], "tripp_alpha", 0.001)
    check("Tripp beta", coef[1], "tripp_beta", 0.005)

    # --- systematic modes quadrature (needs SingleSYS_CovMatrix) ---
    import pandas as pd
    hd_all = pd.read_csv(os.path.join(ROOT, "data", "DES-SN5YR_HD+MetaData.csv"))
    sel = ((hd_all.zHD >= 0.01) & (hd_all.MUERR_FINAL > 0)
           & (hd_all.MUERR_FINAL <= 10)).to_numpy()
    SS = os.path.join(ROOT, "data", "SingleSYS_CovMatrix")
    blocks = ["BS20", "P21SYS1", "P21SYS2", "P21SYS3", "P21_HOSTCOLOR",
              "W22_AGE", "HOSTEFFshift", "SVAHOSTLIB", "CClikelihood",
              "SNNtraining", "SNIRF", "SCONE"]

    def verdict(mu_obs):
        cq = minimize_scalar(lambda e: chi2(mu_obs - mu_modelB(zhd, zhel, e)),
                             bounds=(0.0, 2.0), method="bounded",
                             options={"xatol": 1e-3}).fun
        return cq - fitA(mu_obs)

    base = verdict(MU)
    tot2 = 0.0
    for blk in blocks:
        with gzip.open(os.path.join(SS, blk + ".txt.gz"), "rt") as f:
            n = int(f.readline())
            Ck = np.fromstring(f.read(), sep="\n").reshape(n, n)
        Ck = Ck[np.ix_(sel, sel)]
        w, V = np.linalg.eigh(Ck)
        dvec = V[:, -1] * np.sqrt(max(w[-1], 0.0))
        sp = verdict(MU + dvec) - base
        sm = verdict(MU - dvec) - base
        tot2 += max(abs(sp), abs(sm)) ** 2
    check("modes quadrature", np.sqrt(tot2), "modes_quadrature", 0.05)

    print(f"\n{'ALL CHECKS PASS' if FAIL == 0 else f'{FAIL} CHECK(S) FAILED'}")
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
