#!/usr/bin/env python3
"""
extract_raw_photometry.py -- audit tool: rebuild des_sn5yr_raw_photometry.csv
from DES's public 12 GB data release.

Usage:
    python3 extract_raw_photometry.py \
        --snana-dir /path/to/DES-SN5YR/SNANA/DES_RAW \
        --out data/des_sn5yr_raw_photometry.csv

What this does
--------------
DES-SN5YR releases raw photometry in SNANA ASCII or FITS format, one file
per SN, bundled in a 12 GB tarball. This script walks that directory,
parses each SN's light-curve header + observation table, and writes a
single flat CSV with the schema used by the replication experiments:

    snid, ra, dec, z, mjd, band, flux_nu_jy, flux_nu_jy_err,
    wavelength_eff_nm, snr, survey, source_dataset

The output CSV should be byte-identical (modulo row order) to the cached
copy hosted on Zenodo. Compare with:

    sort data/des_sn5yr_raw_photometry.csv | sha256sum
    sort /path/to/zenodo/copy.csv           | sha256sum

Dependencies: numpy, pandas, optionally astropy (for FITS). ASCII format
is handled with only the standard library + pandas.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# DES griz effective wavelengths (nm), published values.
BAND_NM = {"g": 472.0, "r": 642.0, "i": 784.0, "z": 926.0, "Y": 1005.0}


def _fluxcal_to_jy(fluxcal: np.ndarray) -> np.ndarray:
    """SNANA FLUXCAL is normalized so 27.5 mag = 1.0 count.

    m_AB = 27.5 - 2.5 log10(FLUXCAL)  =>  flux_nu_jy = 3631 * 10^(-0.4 m_AB)
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        m_AB = 27.5 - 2.5 * np.log10(np.where(fluxcal > 0, fluxcal, np.nan))
        return 3631.0 * 10.0 ** (-0.4 * m_AB)


def parse_snana_ascii(path: Path) -> pd.DataFrame | None:
    """Parse a single SNANA .DAT ASCII file. Returns a DataFrame or None."""
    hdr = {}
    obs_rows = []
    with open(path, "r") as f:
        in_obs = False
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("SNID:"):
                hdr["snid"] = line.split(":", 1)[1].strip()
            elif line.startswith("RA:"):
                hdr["ra"] = float(line.split()[1])
            elif line.startswith(("DEC:", "DECL:")):
                hdr["dec"] = float(line.split()[1])
            elif line.startswith("REDSHIFT_FINAL:") or line.startswith("REDSHIFT_HELIO:"):
                hdr["z"] = float(line.split()[1])
            elif line.startswith("SURVEY:"):
                hdr["survey"] = line.split()[1]
            elif line.startswith("VARLIST:"):
                varlist = line.split()[1:]
                in_obs = True
            elif line.startswith("OBS:") and in_obs:
                parts = line.split()[1:]
                obs_rows.append(parts)
            elif line.startswith("END_PHOTOMETRY") or line.startswith("END:"):
                in_obs = False
    if not obs_rows or "snid" not in hdr:
        return None

    df = pd.DataFrame(obs_rows, columns=varlist)
    # Coerce numeric columns
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="ignore")

    # Standardize column names
    rename = {"MJD": "mjd", "BAND": "band", "FLT": "band",
              "FLUXCAL": "fluxcal", "FLUXCALERR": "fluxcalerr"}
    df = df.rename(columns=rename)
    df["snid"] = hdr["snid"]
    df["ra"] = hdr.get("ra", np.nan)
    df["dec"] = hdr.get("dec", np.nan)
    df["z"] = hdr.get("z", np.nan)
    df["survey"] = hdr.get("survey", "DES")
    df["source_dataset"] = "DES-SN5YR"

    df["band"] = df["band"].astype(str).str[0].str.lower()
    df["flux_nu_jy"] = _fluxcal_to_jy(df["fluxcal"].astype(float).values)
    df["flux_nu_jy_err"] = (df["fluxcalerr"].astype(float)
                            * df["flux_nu_jy"] / df["fluxcal"].replace(0, np.nan))
    df["wavelength_eff_nm"] = df["band"].map(BAND_NM)
    with np.errstate(divide="ignore", invalid="ignore"):
        df["snr"] = df["fluxcal"] / df["fluxcalerr"].replace(0, np.nan)
    return df[["snid", "ra", "dec", "z", "mjd", "band",
               "flux_nu_jy", "flux_nu_jy_err", "wavelength_eff_nm",
               "snr", "survey", "source_dataset"]]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snana-dir", type=Path, required=True,
                    help="Directory containing DES SNANA per-SN .DAT files.")
    ap.add_argument("--out", type=Path, required=True,
                    help="Output CSV path.")
    ap.add_argument("--pattern", default="*.DAT",
                    help="Glob for per-SN files (default: *.DAT).")
    args = ap.parse_args()

    files = sorted(args.snana_dir.rglob(args.pattern))
    if not files:
        sys.exit(f"No files matching {args.pattern} under {args.snana_dir}")
    print(f"Found {len(files)} SNANA files under {args.snana_dir}")

    dfs = []
    fail = 0
    for i, f in enumerate(files):
        if i % 500 == 0:
            print(f"  [{i:5d}/{len(files)}] parsing... (fails={fail})")
        try:
            df = parse_snana_ascii(f)
        except Exception as e:
            fail += 1
            if fail < 10:
                print(f"    parse failure in {f.name}: {e}")
            continue
        if df is not None and len(df):
            dfs.append(df)

    print(f"Concatenating {len(dfs)} per-SN tables...")
    out = pd.concat(dfs, ignore_index=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"Wrote {len(out):,} rows -> {args.out}")
    print(f"Unique SNe: {out['snid'].nunique()}")


if __name__ == "__main__":
    main()
