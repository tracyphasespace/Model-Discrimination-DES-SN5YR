#!/usr/bin/env python3
"""Camilleri amplitude comparison (Draft 4.1, Section 6).

Camilleri et al. (2024, MNRAS 533, 2615) tested the BBC REFERENCE-cosmology
dependence (the cosmology used to generate the bias-correction simulations,
distinct from the FIDUCIAL cosmology used inside the BBC nuisance fit) with
simulation reruns spanning distance-modulus variations of ~0.15 mag and
found the effect subdominant to DES statistical uncertainties.

SCOPE: this is an AMPLITUDE comparison -- each probe's maximum
offset-removed |dmu| from LCDM against the ~0.15 mag scale of the Camilleri
reference variations. Their formal validation criterion is a SHAPE envelope
(the boundary reference curves of their Fig. 4e); the rigorous pointwise
containment test against those boundary cosmologies is the named upgrade.
Expected: released-best-fit Model B departs by at most ~0.055 mag (well
below the tested amplitude); Einstein-de Sitter by ~0.276 mag (well above),
consistent with its far larger leverage.
"""
import json
import os

import numpy as np
from bracket import load, make_chi2, mu_modelB, mu_lcdm

# External literature constant (Camilleri et al. 2024): tested amplitude scale
CAMILLERI_TESTED_AMPLITUDE_MAG = 0.15
# Fitted coordinates from the frozen-results manifest (tabulated values)
_F = json.load(open(os.path.join(os.path.dirname(__file__), "..",
                                 "frozen_results.json")))
ETA_REL = float(_F["eta_released"])
OM_REL = float(_F["om_released"])
ETA_PRE = float(_F["eta_prebbc"])
OM_PRE = float(_F["om_prebbc"])


def main():
    hd, C = load()
    Ci = np.linalg.inv(C)
    u = np.ones(len(hd))
    uCu = u @ Ci @ u
    zhd, zhel = hd.zHD.to_numpy(), hd.zHEL.to_numpy()

    def proj(x):
        return x - (u @ Ci @ x) / uCu * u

    cases = [
        (f"Model B, released fit (eta={ETA_REL}) vs LCDM Om={OM_REL}",
         mu_modelB(zhd, zhel, ETA_REL) - mu_lcdm(zhd, zhel, OM_REL)),
        (f"Model B, released fit (eta={ETA_REL}) vs Planck Om=0.315",
         mu_modelB(zhd, zhel, ETA_REL) - mu_lcdm(zhd, zhel, 0.315)),
        (f"Model B, pre-BBC fit (eta={ETA_PRE}) vs LCDM Om={OM_PRE}",
         mu_modelB(zhd, zhel, ETA_PRE) - mu_lcdm(zhd, zhel, OM_PRE)),
        (f"Einstein-de Sitter (Om=1) vs LCDM Om={OM_REL}",
         mu_lcdm(zhd, zhel, 1.0) - mu_lcdm(zhd, zhel, OM_REL)),
    ]
    print(f"Camilleri et al. tested reference-cosmology amplitude scale: "
          f"~{CAMILLERI_TESTED_AMPLITUDE_MAG} mag\n")
    for name, d in cases:
        dt = proj(d)
        mx = np.abs(dt).max()
        rel = mx / CAMILLERI_TESTED_AMPLITUDE_MAG
        pos = ("BELOW tested amplitude" if mx < CAMILLERI_TESTED_AMPLITUDE_MAG
               else "ABOVE tested amplitude")
        print(f"{name}:\n  offset-removed max|dmu| = {mx:.3f} mag "
              f"({rel:.2f}x the tested ~0.15 mag scale) -> {pos}\n")
    print("Caveats: amplitude comparison only, not the pointwise Fig-4e shape")
    print("containment test (the named upgrade); the Camilleri envelope used")
    print("FLRW-family reference simulations, so non-FLRW coverage is inferential.")


if __name__ == "__main__":
    main()
