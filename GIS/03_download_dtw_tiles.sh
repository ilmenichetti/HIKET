#!/usr/bin/env bash
# =============================================================================
# 03_download_dtw_tiles.sh — Luke depth-to-water (DTW) 2 m, 2023, all five
# stream-initiation thresholds (0.5/1/2/4/10 ha = very moist ... dry), ONLY the
# tiles that cover a plot (list from 02_plots_and_dtw_tiles.R). ~9 MB per tile.
# Source: https://www.nic.funet.fi/index/geodata/luke/dtw/2023/ (Salmivaara, Luke)
# Re-runnable: a tile is skipped only if its size matches the server's (a
# truncated GeoTIFF can still open). 4 parallel downloads.
# Run from repo root:  bash GIS/03_download_dtw_tiles.sh
# =============================================================================
set -u
G="${HIKET_GIS_DIR:-/Volumes/NextGenC_SS/HIKET_GIS}"
D="https://www.nic.funet.fi/index/geodata/luke/dtw/2023"
L="$G/plots/dtw_tiles_needed.txt"; [ -f "$L" ] || { echo "run 02 first"; exit 1; }
get_tile() {   # $1 = tile name; uses $D $v $dst
  local f="$dst/$1" u="$D/DTW_INT_CMv2_$v/$1" r
  if [ -s "$f" ]; then                        # complete = same size as on the server
    r=$(curl -sI -m 60 "$u" | grep -i "^content-length" | tr -dc "0-9")
    [ -n "$r" ] && [ "$(stat -f%z "$f")" = "$r" ] && return 0
    echo "refetch $v $1 (local $(stat -f%z "$f"), remote $r)"
  fi
  curl -sS --retry 5 --retry-all-errors --retry-delay 10 -o "$f" "$u" || echo "FAIL $v $1"
}
export -f get_tile
for v in 050 1 2 4 10; do
  dst="$G/raw/dtw_2m/DTW_INT_CMv2_$v"; mkdir -p "$dst"
  export D v dst
  # keep the xargs command short: macOS xargs -I caps the command at 255 bytes
  xargs -P 4 -I{} bash -c 'get_tile "$1"' _ {} < "$L"
  echo "variant $v: $(ls "$dst"/*.tif 2>/dev/null | wc -l) tiles"
done
echo "DONE $(date)"
