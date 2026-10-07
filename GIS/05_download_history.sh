#!/usr/bin/env bash
# =============================================================================
# 05_download_history.sh — early-20th-century human-impact maps for Finland
# Aakala, Kulha & Kuuluvainen 2023, Landscape Ecology 38:2417-2431; figshare
# 10.6084/m9.figshare.23257562, CC BY 4.0, all EPSG:3067:
#   slash_and_burn.gpkg     share of area under slash-and-burn ("kaskimaita"),
#                           commonness 1860 / 1913 (Heikinheimo 1915 classes)
#   population_1925.gpkg    population points (field Asukkaita2 = persons/point)
#   parishes_1902.gpkg, railroads_1925.gpkg, state_forests_1925.gpkg
# Run from repo root:  bash GIS/05_download_history.sh
# =============================================================================
set -u
G="${HIKET_GIS_DIR:-/Volumes/NextGenC_SS/HIKET_GIS}"; dst="$G/raw/aakala2023_history"; mkdir -p "$dst"
curl -s -m 60 https://api.figshare.com/v2/articles/23257562 > "$dst/figshare_metadata.json"
python3 - "$dst" <<'PY'
import json, sys, subprocess
d = sys.argv[1]; r = json.load(open(d + "/figshare_metadata.json"))
for f in r["files"]:
    subprocess.run(["curl", "-sSL", "--retry", "5", "-o", f"{d}/{f['name']}", f["download_url"]], check=False)
    print("got", f["name"])
PY
