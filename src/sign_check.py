#!/usr/bin/env python3
"""Establish the sign convention of biasCor_mu empirically (paper §7 Step 3).

Reconstruct the Tripp estimator mB + alpha*x1 - beta*c + M0 by least squares
against three hypotheses for how biasCor_mu enters the released MU. The
hypothesis that reconstructs with textbook SALT coefficients (alpha ~ 0.15,
beta ~ 3.1) and mmag-level residuals reveals the convention.

Expected result (DES-SN5YR v1.2): MU + biasCor matches at rms ~0.017 mag
with alpha = 0.155, beta = 3.12  =>  DES SUBTRACTS the correction, so the
UNCORRECTED vector is MU + biasCor_mu.

This check caught a sign error in this paper's own first analysis. Run it
before anything else.
"""
import os
import numpy as np
import pandas as pd

DATA = os.environ.get("BRACKET_DATA", os.path.join(os.path.dirname(__file__), "..", "data"))

hd = pd.read_csv(os.path.join(DATA, "DES-SN5YR_HD+MetaData.csv"))
sel = (hd.zHD >= 0.01) & (hd.MUERR_FINAL > 0) & (hd.MUERR_FINAL <= 10)
d = hd[sel]
print(f"N = {sel.sum()} (published analyses use 1768)")

X = np.column_stack([d.x1, -d.c, np.ones(len(d))])
for name, target in [
    ("MU as-is            ", d.MU),
    ("MU - biasCor        ", d.MU - d.biasCor_mu),
    ("MU + biasCor        ", d.MU + d.biasCor_mu),
]:
    y = (target - d.mB).to_numpy()
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    rms = (y - X @ coef).std()
    print(f"{name} alpha={coef[0]:+.3f} beta={coef[1]:+.3f} rms={rms:.4f} mag")

print("\n(The design matrix column is -c, so the printed beta is +beta.)")
print("Smallest rms with textbook (alpha~0.15, beta~3.1) identifies the")
print("Tripp-side vector; that vector is the UNCORRECTED one.")
