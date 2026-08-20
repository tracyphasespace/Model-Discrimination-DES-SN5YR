#!/usr/bin/env python3
"""Second-probe leverage tests and diagnostic sensitivities (review items).

1. rho_C evaluated at BOTH best-fit definitions of the model-separation
   vector d (released-vector fits and pre-BBC-vector fits).
2. Einstein-de Sitter probe (Om=1, no free shape parameter) vs free-Om LCDM:
   leverage L_BBC^EdS.
3. Flat wCDM probe with w=-0.5 fixed, Om as the single shape parameter:
   leverage L_BBC^wCDM.
4. Tripp regression coefficients with standard errors (sign-check follow-up).

The broader conclusion -- that BBC sensitivity is a property of
wide-separation comparisons rather than of the non-expanding probe -- is
licensed only if these probes reproduce substantial leverage.
"""
import numpy as np
from bracket import load, make_chi2, mu_modelB, mu_lcdm


def mu_wcdm(z, om, w=-0.5, n_grid=400):
    zs = np.linspace(0, z.max(), n_grid)
    zm = 0.5 * (zs[1:] + zs[:-1])
    E = np.sqrt(om * (1 + zm) ** 3 + (1 - om) * (1 + zm) ** (3 * (1 + w)))
    integ = np.concatenate([[0], np.cumsum(np.diff(zs) / E)])
    return 5 * np.log10((1 + z) * np.interp(z, zs, integ))


def main():
    hd, C = load()
    chi2 = make_chi2(C)
    Ci = np.linalg.inv(C)
    u = np.ones(len(hd))
    uCu = u @ Ci @ u
    z = hd.zHD.to_numpy()
    MU = hd.MU.to_numpy()
    b = hd.biasCor_mu.to_numpy()
    pre = MU + b

    def proj(x):
        return x - (u @ Ci @ x) / uCu * u

    def fitA(mu_obs):
        oms = np.arange(0.02, 1.2001, 0.01)
        cl = [chi2(mu_obs - mu_lcdm(z, o)) for o in oms]
        i = int(np.argmin(cl))
        return oms[i], cl[i]

    # --- 1. rho_C at both d definitions ---
    bt = proj(b)
    nb = np.sqrt(bt @ Ci @ bt)
    for tag, (om, eta) in [("released fits (0.35, 0.30)", (0.35, 0.30)),
                           ("pre-BBC fits  (0.50, 0.05)", (0.50, 0.05))]:
        d = proj(mu_modelB(z, eta) - mu_lcdm(z, om))
        nd = np.sqrt(d @ Ci @ d)
        rho = (bt @ Ci @ d) / (nb * nd)
        print(f"rho_C [{tag}]: {rho:+.3f}  (|b|={nb:.2f}, |d|={nd:.2f})")

    # --- 2 & 3. second probes vs free-Om LCDM ---
    for name, mu_probe_fn, has_param in [
            ("EdS (Om=1, no free shape)", lambda mo: chi2(mo - mu_lcdm(z, 1.0)), False),
            ("wCDM (w=-0.5, Om fitted)",
             lambda mo: min(chi2(mo - mu_wcdm(z, o)) for o in np.arange(0.02, 1.0001, 0.01)),
             True)]:
        ds = {}
        for vec_tag, mo in [("released", MU), ("preBBC", pre)]:
            omA, c2A = fitA(mo)
            c2P = mu_probe_fn(mo)
            ds[vec_tag] = c2P - c2A
            print(f"{name} | {vec_tag:9s}: Dchi2(probe - LCDM) = {c2P - c2A:+8.1f} "
                  f"(LCDM Om={omA:.2f})")
        L = ds["released"] - ds["preBBC"]
        print(f"{name} | L_BBC = {L:+.1f}\n")

    # --- 4. Tripp regression with uncertainties ---
    y = (MU + b - hd.mB.to_numpy())
    X = np.column_stack([hd.x1.to_numpy(), -hd.c.to_numpy(), np.ones(len(hd))])
    coef, res, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ coef
    dof = len(y) - 3
    s2 = (r @ r) / dof
    cov = s2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    print(f"Tripp regression on MU+biasCor: alpha = {coef[0]:.5f} +/- {se[0]:.5f}, "
          f"beta = {-coef[1]:.5f} +/- {se[1]:.5f}, rms = {r.std():.4f} mag")
    print("DES-reported global values: alpha = 0.16087, beta = 3.11780 (gamma 0.03754")
    print("not modeled here; its omission is part of the residual rms).")


if __name__ == "__main__":
    main()
