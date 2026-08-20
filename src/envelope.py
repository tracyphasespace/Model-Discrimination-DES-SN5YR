#!/usr/bin/env python3
"""Camilleri envelope comparison (Draft 4.1, Section 6).

Camilleri et al. (2024, MNRAS 533, 2615) tested the BBC REFERENCE-cosmology
dependence (the cosmology used to generate the bias-correction simulations,
distinct from the FIDUCIAL cosmology used inside the BBC nuisance fit) with
simulation reruns spanning distance-modulus variations of ~0.15 mag and
found the effect subdominant to DES statistical uncertainties.

This script measures where each probe's offset-removed shape difference from
LCDM sits relative to that envelope. Expected: released-best-fit Model B
(eta ~ 0.30) lies ~3x INSIDE the envelope (max |dmu| ~ 0.055 mag), so the
near-degenerate released comparison operates within DES's validated
reference-cosmology range; Einstein-de Sitter lies far outside, consistent
with its far larger leverage.
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
        pos = "INSIDE" if mx < CAMILLERI_ENVELOPE_MAG else "OUTSIDE"
        print(f"{name}:\n  offset-removed max|dmu| = {mx:.3f} mag "
              f"({rel:.2f}x envelope) -> {pos}\n")
    print("Caveats: the envelope was established with FLRW-family reference")
    print("simulations; formal coverage of non-FLRW functional forms is an")
    print("extrapolation, and the pre-BBC endpoint's fit may sit differently.")


if __name__ == "__main__":
    main()
