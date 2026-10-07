#!/usr/bin/env bash
# =============================================================================
# 01_download_national.sh — national GIS layers that do not depend on plot
# locations, downloaded WHOLE into the HIKET GIS store on the external drive.
#
#   Metsäkeskus forest use declarations (metsänkäyttöilmoitukset), by region,
#     GeoPackage, CC BY 4.0 — harvest notices with cutting type, completion
#     year and forest-damage qualifier; usable from 2004
#   Metsäkeskus KEMERA / Metka subsidised and completed works, by region
#   Luke Topographic Wetness Index, 16 m, 2016 (single GeoTIFF + overview)
#
# Store: $HIKET_GIS_DIR (default /Volumes/NextGenC_SS/HIKET_GIS). Re-runnable:
# curl -C - resumes partial files, complete files are skipped. Metsäkeskus files
# come 6 at a time (their storage gives ~80 kB/s per connection).
# Run from repo root:  bash GIS/01_download_national.sh
# =============================================================================
set -u
G="${HIKET_GIS_DIR:-/Volumes/NextGenC_SS/HIKET_GIS}"
[ -d "$G" ] || { echo "GIS store not found: $G"; exit 1; }

fetch() {  # url dest
  local url="$1" dest="$2"
  local remote; remote=$(curl -sIL -m 60 "$url" | grep -i '^content-length' | tail -1 | tr -dc '0-9')
  if [ -f "$dest" ] && [ -n "$remote" ] && [ "$(stat -f%z "$dest")" = "$remote" ]; then
    echo "ok   $(basename "$dest")"; return; fi
  echo "get  $(basename "$dest")"
  local i; for i in $(seq 1 50); do              # resume dropped connections
    curl -sSL --retry 5 --retry-all-errors --retry-delay 10 -C - -o "$dest" "$url" && return
    sleep 15
  done; echo "FAIL $url"
}

# --- Luke TWI first: fast server, single file -----------------------------------
TW="https://www.nic.funet.fi/index/geodata/luke/twi"
for f in TWI_16m_Finland_NA_lakes_int.tif TWI_16m_Finland_NA_lakes_int.tif.ovr \
         TWI_16m_Finland_NA_lakes_int.tif.aux.xml Readme_Topographical_Wetness_Index_2016.txt \
         TWI_metadata_description.docx supplementtwimetadatadescription.pdf; do
  fetch "$TW/$f" "$G/raw/twi_16m/$f"
done

# --- Metsäkeskus: ~80 kB/s per connection, so 6 files at a time ----------------
MK="https://avoin.metsakeskus.fi/aineistot"
export -f fetch
for prod in Metsankayttoilmoitukset:metsakeskus_declarations Kemera:metsakeskus_kemera; do
  src=${prod%%:*}; dst="$G/raw/${prod##*:}"; mkdir -p "$dst"
  for f in $(curl -s -m 60 "$MK/$src/Maakunta/" | grep -oE 'href="[^"]+\.zip"' | sed 's/href="//;s/"$//'); do
    name=$(python3 -c 'import urllib.parse,sys;print(urllib.parse.unquote(sys.argv[1]))' "$(basename "$f")")
    enc=$(python3 -c 'import urllib.parse,sys;print(urllib.parse.quote(sys.argv[1]))' "$name")
    printf '%s\t%s\n' "$MK/$src/Maakunta/$enc" "$dst/$name"
  done
done | xargs -P 6 -L 1 bash -c 'fetch "$0" "$1"'

echo "DONE $(date)"
