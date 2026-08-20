#!/usr/bin/env python3
"""Projection and decomposition diagnostics (paper Sections 5.8-5.9).

Computes: rho_C (covariance-weighted cosine between the BBC vector and the
model-separation vector, offset mode projected out), the C^-1 norms, and the
fixed-shape vs adaptive leverage decomposition. Identity check:
L_fixed = 2 * rho_C * |b| * |d| (exact, from the quadratic form).

Expected (v1.2): rho_C = -0.065, |b| = 29.8, |d| = 2.51,
L_fixed = -9.7, L_adapt = -2.0 (L_total = -11.7).
"""
import numpy as np
from bracket import load, mu_modelB, mu_lcdm

ETA_REL, OM_REL = 0.30, 0.35   # released-vector best fits (bracket.py)
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

    mA, mB = mu_lcdm(z, OM_REL), mu_modelB(z, ETA_REL)
    b, d = proj(hd.biasCor_mu.to_numpy()), proj(mB - mA)
    nb, nd = np.sqrt(b @ Ci @ b), np.sqrt(d @ Ci @ d)
    rho = (b @ Ci @ d) / (nb * nd)
    print(f"rho_C = {rho:+.3f}   |b|_Cinv = {nb:.2f}   |d|_Cinv = {nd:.2f}")

    MU = hd.MU.to_numpy()
    pre = MU + hd.biasCor_mu.to_numpy()
    L_fixed = (q(MU - mB) - q(MU - mA)) - (q(pre - mB) - q(pre - mA))
    print(f"L_fixed = {L_fixed:+.2f}   (identity 2*rho*|b|*|d| = {2*rho*nb*nd:+.2f})")
    print(f"L_adapt = {L_TOTAL - L_fixed:+.2f}   (L_total = {L_TOTAL})")
    print("\nReading: the BBC vector is ~94% covariance-orthogonal to the model")
    print("difference — but at |b| ~ 30, the residual 6.5% projection alone")
    print("carries ~6x the released verdict. For the verdict to be safe, the")
    print("CONDITIONED part of b must satisfy 2|b_cond . d|_Cinv < 1.6.")
