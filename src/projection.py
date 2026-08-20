#!/usr/bin/env python3
"""Projection and decomposition diagnostics (paper Sections 5.8-5.9).

Computes: rho_C (covariance-weighted cosine between the BBC vector and the
model-separation vector, offset mode projected out), the C^-1 norms, and the
fixed-shape vs adaptive leverage decomposition. Identity check:
L_fixed = 2 * rho_C * |b| * |d| (exact, from the quadratic form).

Expected (v1.2, official zHD/zHEL convention): rho_C = -0.068,
|b| = 29.76, |d| = 2.51, L_fixed = -10.1, L_adapt = -1.6 (L_total = -11.7).
"""
import numpy as np
from bracket import load, mu_modelB, mu_lcdm

ETA_REL, OM_REL = 0.297, 0.352  # released-vector best fits (bracket.py, DES zHD/zHEL convention)
L_TOTAL = -11.7


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
    print(f"L_adapt = {L_TOTAL - L_fixed:+.2f}   (L_total = {L_TOTAL})")
    m = 1.1  # released-vector margin under the DES zHD/zHEL convention
    print(f"\nReading: the vectors are nearly covariance-orthogonal (small cosine),")
    print(f"but |b| is large enough that the residual projection dominates the")
    print(f"released margin. General requirement on the conditioned component:")
    print(f"|rho_cond| * |b_cond|_Cinv < margin/(2*|d|) = {m/(2*nd):.3f}")

if __name__ == "__main__":
    main()
