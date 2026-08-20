#!/usr/bin/env bash
# Fetch the twelve single-systematic covariance blocks used by
# src/bbc_systematic_modes.py (DES-SN5YR release, tag v1.2).
set -euo pipefail
mkdir -p SingleSYS_CovMatrix && cd SingleSYS_CovMatrix
B="https://raw.githubusercontent.com/des-science/DES-SN5YR/v1.2/4_DISTANCES_COVMAT/SingleSYS_CovMatrix"
for f in BS20 P21SYS1 P21SYS2 P21SYS3 P21_HOSTCOLOR W22_AGE HOSTEFFshift \
         SVAHOSTLIB CClikelihood SNNtraining SNIRF SCONE; do
  [ -f "$f.txt.gz" ] || curl -fL --retry 3 -o "$f.txt.gz" "$B/$f.txt.gz"
done
sha256sum *.txt.gz
echo "done."
