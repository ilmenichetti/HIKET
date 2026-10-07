#!/usr/bin/env bash
# =============================================================================
# 04_download_emep.sh — EMEP MSC-W rv5.3 trend runs with CAMS-REG emissions
# (Simpson, met.no; Zenodo 10.5281/zenodo.12580842, CC BY 4.0), 0.1 x 0.1 deg,
# ecosystem deposition. Whole archive (4 GB) kept in the GIS store.
# Run from repo root:  bash GIS/04_download_emep.sh
# =============================================================================
set -u
G="${HIKET_GIS_DIR:-/Volumes/NextGenC_SS/HIKET_GIS}"; dst="$G/raw/emep_ndep"; mkdir -p "$dst"
URL="https://zenodo.org/api/records/12580842/files/EMEP_Files.zip/content"
# A dropped connection (e.g. curl error 56) is resumed from the last byte, up to 50 times.
for i in $(seq 1 50); do
  curl -sSL --retry 5 --retry-all-errors --retry-delay 10 -C - -o "$dst/EMEP_Files.zip" "$URL" && break
  echo "attempt $i interrupted at $(stat -f%z "$dst/EMEP_Files.zip") bytes; resuming"; sleep 15
done
unzip -l "$dst/EMEP_Files.zip" > "$dst/EMEP_Files_listing.txt"
echo "DONE $(date)"
