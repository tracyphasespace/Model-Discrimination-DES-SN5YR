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
import numpy as np
from bracket import load, make_chi2, mu_modelB, mu_lcdm

CAMILLERI_ENVELOPE_MAG = 0.15


def main():
    hd, C = load()
    Ci = np.linalg.inv(C)
    u = np.ones(len(hd))
    uCu = u @ Ci @ u
    zhd, zhel = hd.zHD.to_numpy(), hd.zHEL.to_numpy()

    def proj(x):
        return x - (u @ Ci @ x) / uCu * u

    cases = [
        ("Model B, released fit (eta=0.297) vs LCDM Om=0.352",
         mu_modelB(zhd, zhel, 0.297) - mu_lcdm(zhd, zhel, 0.352)),
        ("Model B, released fit (eta=0.297) vs Planck Om=0.315",
         mu_modelB(zhd, zhel, 0.297) - mu_lcdm(zhd, zhel, 0.315)),
        ("Model B, pre-BBC fit (eta=0.049) vs LCDM Om=0.499",
         mu_modelB(zhd, zhel, 0.049) - mu_lcdm(zhd, zhel, 0.499)),
        ("Einstein-de Sitter (Om=1) vs LCDM Om=0.352",
         mu_lcdm(zhd, zhel, 1.0) - mu_lcdm(zhd, zhel, 0.352)),
    ]
    print(f"Camilleri et al. validated reference-cosmology envelope: "
          f"~{CAMILLERI_ENVELOPE_MAG} mag\n")
    for name, d in cases:
        dt = proj(d)
        mx = np.abs(dt).max()
        rel = mx / CAMILLERI_ENVELOPE_MAG
        pos = ("BELOW tested amplitude" if mx < CAMILLERI_ENVELOPE_MAG
               else "ABOVE tested amplitude")
        print(f"{name}:\n  offset-removed max|dmu| = {mx:.3f} mag "
              f"({rel:.2f}x the tested ~0.15 mag scale) -> {pos}\n")
    print("Caveats: amplitude comparison only, not the pointwise Fig-4e shape")
    print("containment test (the named upgrade); the Camilleri envelope used")
    print("FLRW-family reference simulations, so non-FLRW coverage is inferential.")


if __name__ == "__main__":
    main()
