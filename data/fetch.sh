#!/usr/bin/env bash
# Fetch the public DES-SN5YR release files needed for the bracket analysis.
# Pinned to release tag v1.2 — do NOT use main (Dovekie recalibration differs).
set -euo pipefail
B="https://raw.githubusercontent.com/des-science/DES-SN5YR/v1.2/4_DISTANCES_COVMAT"
for f in DES-SN5YR_HD+MetaData.csv DES-SN5YR_HD.csv STAT+SYS.txt.gz README.md; do
  [ -f "$f" ] || curl -fL --retry 3 -o "$f" "$B/$f"
done
echo "done. NOTE: STAT+SYS.txt.gz contains SYSTEMATICS ONLY despite the name;"
echo "the statistical term is in MUERR_FINAL (see README.md line ~81)."
# Raw griz light curves for all 8,293 candidates (optional; used by the
# extended replication): rebuild from the DES public release with
# ../tools/extract_raw_photometry.py, or fetch the cached CSV (Zenodo DOI TBD).
# Verify frozen v1.2 hashes (raw-photometry hash covers the in-repo copy)
( cd "$(dirname "$0")/.." && sha256sum -c data/SHA256SUMS --ignore-missing ) \
  && echo "hash verification OK"
