#!/usr/bin/env python3
"""Projection and decomposition diagnostics (paper Sections 5.8-5.9).

Computes: rho_C (covariance-weighted cosine between the BBC vector and the
model-separation vector, offset mode projected out), the C^-1 norms, and the
fixed-shape vs adaptive leverage decomposition. Identity check:
L_fixed = 2 * rho_C * |b| * |d| (exact, from the quadratic form).

The model-separation vector d is evaluated at the TABULATED THREE-DECIMAL
released-vector best-fit values (read from frozen_results.json), per the
manuscript's Section 5.8; using full-precision optimizer outputs instead
shifts L_fixed by ~0.2. L_total is computed from the two fitted verdicts,
not hard-coded. Expected (v1.2, official zHD/zHEL convention):
rho_C = -0.068, |b| = 29.76, |d| = 2.51, L_fixed = -10.1, L_adapt = -1.6.
"""
import numpy as np
import json
import os

from bracket import load, make_chi2, mu_modelB, mu_lcdm, fit

_F = json.load(open(os.path.join(os.path.dirname(__file__), "..",
                                 "frozen_results.json")))
# Tabulated three-decimal released-vector best-fit values (the evaluation
# point of the projection, per manuscript Section 5.8) — from the manifest.
ETA_REL = float(_F["eta_released"])
OM_REL = float(_F["om_released"])


def main():
    hd, C = load()
    Ci = np.linalg.inv(C)
    u = np.ones(len(hd))
    uCu = u @ Ci @ u
    z = hd.zHD.to_numpy()

    def proj(x):
        return x - (u @ Ci @ x) / uCu * u

    def q(r):
        r = proj(r)
        return r @ Ci @ r

    zhel = hd.zHEL.to_numpy()
    mA, mB = mu_lcdm(z, zhel, OM_REL), mu_modelB(z, zhel, ETA_REL)
    b, d = proj(hd.biasCor_mu.to_numpy()), proj(mB - mA)
    nb, nd = np.sqrt(b @ Ci @ b), np.sqrt(d @ Ci @ d)
    rho = (b @ Ci @ d) / (nb * nd)
    print(f"rho_C = {rho:+.3f}   |b|_Cinv = {nb:.2f}   |d|_Cinv = {nd:.2f}")

    MU = hd.MU.to_numpy()
    pre = MU + hd.biasCor_mu.to_numpy()
    L_fixed = (q(MU - mB) - q(MU - mA)) - (q(pre - mB) - q(pre - mA))
    print(f"L_fixed = {L_fixed:+.2f}   (identity 2*rho*|b|*|d| = {2*rho*nb*nd:+.2f})")

    # L_total and the released margin computed from the two fitted verdicts
    chi2 = make_chi2(C)
    r_rel = fit(chi2, z, zhel, MU)
    r_pre = fit(chi2, z, zhel, pre)
    m = r_rel["chi2_B"] - r_rel["chi2_A"]           # released margin
    L_total = m - (r_pre["chi2_B"] - r_pre["chi2_A"])
    print(f"L_adapt = {L_total - L_fixed:+.2f}   (L_total = {L_total:+.2f}, computed)")
    print(f"\nReading: the vectors are nearly covariance-orthogonal (small cosine),")
    print(f"but |b| is large enough that the residual projection dominates the")
    print(f"released margin. General requirement on the conditioned component:")
    print(f"|rho_cond| * |b_cond|_Cinv < margin/(2*|d|) = {m/(2*nd):.3f}")

if __name__ == "__main__":
    main()
