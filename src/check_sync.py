#!/usr/bin/env python3
"""Cross-surface consistency gate: assert that the frozen numbers in
frozen_results.json appear in README.md and paper/main.tex. Run by CI on
every change to src/, README, or the manuscript, so no surface can drift
from the frozen results without failing the build."""
import json
import os
import sys

ROOT = os.path.join(os.path.dirname(__file__), "..")
M = json.load(open(os.path.join(ROOT, "frozen_results.json")))

# keys -> surfaces they must appear in (README, main.tex)
REQUIRED = {
    "README.md": ["delta_chi2_released", "delta_chi2_prebbc", "L_BBC",
                  "leverage_ratio_words", "L_wcdm", "L_eds",
                  "envelope_probe_mag", "envelope_eds_mag",
                  "modes_quadrature", "rho_C", "L_fixed", "L_adapt"],
    "paper/main.tex": ["delta_chi2_released", "delta_chi2_prebbc", "L_BBC",
                       "om_released", "eta_released", "om_prebbc",
                       "eta_prebbc", "chi2dof_released", "chi2dof_prebbc",
                       "rho_C", "rho_C_alt", "b_norm", "d_norm",
                       "L_fixed", "L_adapt", "L_wcdm", "L_eds",
                       "modes_quadrature", "envelope_probe_mag",
                       "envelope_eds_mag", "camilleri_scale_mag",
                       "tripp_alpha", "biascor_slope",
                       "N_default", "N_official"],
}

fail = 0
for surface, keys in REQUIRED.items():
    text = open(os.path.join(ROOT, surface)).read()
    for k in keys:
        v = M[k]
        # match with or without leading sign/LaTeX minus
        variants = [v, v.replace("-", "\u2212"), v.replace("-", ""),
                    v.replace("+", "")]
        if not any(x in text for x in variants):
            print(f"MISSING in {surface}: {k} = {v}")
            fail += 1
if fail:
    print(f"\n{fail} frozen value(s) missing -- surfaces out of sync.")
    sys.exit(1)
print("All frozen values present in README.md and paper/main.tex.")
