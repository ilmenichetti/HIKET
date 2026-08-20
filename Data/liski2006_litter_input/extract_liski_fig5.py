#!/usr/bin/env python3
"""
Extract the litter-input series from Fig. 5 of Liski et al. (2006), Ann. For. Sci. 63:687-697,
"Carbon accumulation in Finland's forests 1922-2004".

NOT a digitisation. The figure is VECTOR in the publisher PDF, so the curve vertices are read
exactly out of the drawing; the only approximation is that a curve is polyline-simplified by the
original plotting software. Axis calibration is likewise exact: it is taken from the tick-mark
coordinates in the same vector path that draws the frame.

Usage:  python3 extract_liski_fig5.py <paper.pdf> <outdir>
Requires: pdftocairo (poppler).
"""
import re, sys, subprocess, csv, os
from collections import defaultdict

PAGE = 5                       # journal p. 691 = PDF page 5
YFLIP = 799.0                  # page is 629 x 799 pt; pdftocairo emits matrix(1,0,0,-1,0,799)
# Axis calibration, read from the frame/tick path on that page (see README).
YT = [247.74, 227.07, 206.06, 185.44, 164.42, 143.80, 122.79, 102.16]   # ticks
YV = [-5, 0, 5, 10, 15, 20, 25, 30]                                     # Tg C /yr
XT = [347.51, 375.27, 403.02, 430.78, 458.54, 485.90, 513.66, 541.42, 569.18]
XV = [1922, 1932, 1942, 1952, 1962, 1972, 1982, 1992, 2002]
# Series identified by stroke grey level + width (they are distinct, so no eyeballing).
SERIES = {(65.5, 1.189): "tree_litter",
          (22.4, 1.189): "harvest_residues",
          (13.7, 1.586): "ground_vegetation",   # dotted -> many short dashes
          (13.7, 0.793): "natural_mortality",   # dashed, near-constant
          (31.0, 0.396): "luc_transfer"}        # transfers between forest and other land uses

def main(pdf, outdir):
    svg = os.path.join(outdir, "_liski_p5.svg")
    subprocess.run(["pdftocairo", "-svg", "-f", str(PAGE), "-l", str(PAGE), pdf, svg], check=True)
    s = open(svg).read()
    ky = (YV[-1] - YV[0]) / (YT[-1] - YT[0]); kx = (XV[-1] - XV[0]) / (XT[-1] - XT[0])
    tg = lambda y: YV[1] + (y - YT[1]) * ky
    yr = lambda x: XV[0] + (x - XT[0]) * kx
    buckets = defaultdict(list)
    for m in re.finditer(r"<path[^>]*/?>", s):
        p = m.group(0); d = re.search(r'd="([^"]*)"', p)
        if not d: continue
        P = [(float(a), YFLIP - float(b)) for a, b in
             re.findall(r"[ML]\s*(-?\d+\.?\d*)[ ,]+(-?\d+\.?\d*)", d.group(1))]
        g = re.search(r'stroke="rgb\(([\d.]+)%', p); w = re.search(r'stroke-width="([\d.]+)"', p)
        if not P or not g: continue
        xs = [q[0] for q in P]; ys = [q[1] for q in P]
        if min(xs) < 345 or max(xs) > 576 or min(ys) < 95 or max(ys) > 255: continue   # Fig-5 panel
        key = (round(float(g.group(1)), 1), float(w.group(1)) if w else None)
        if key in SERIES: buckets[SERIES[key]].extend((yr(x), tg(y)) for x, y in P)
    # collapse to one value per year (dashed series contribute several dashes per year)
    out = {}
    for nm, P in buckets.items():
        byyear = defaultdict(list)
        for t, v in P: byyear[round(t)].append(v)
        out[nm] = {y: sum(v) / len(v) for y, v in byyear.items()}
    years = list(range(1922, 2005))
    def interp(d):                                   # linear fill across gaps
        ks = sorted(d)
        res = []
        for y in years:
            if y in d: res.append(d[y]); continue
            lo = [k for k in ks if k < y]; hi = [k for k in ks if k > y]
            if not lo or not hi: res.append(d[ks[0] if not lo else ks[-1]]); continue
            a, b = lo[-1], hi[0]
            res.append(d[a] + (d[b] - d[a]) * (y - a) / (b - a))
        return res
    cols = {nm: interp(d) for nm, d in out.items()}
    # TREE BASIS = tree litter + harvest residues + natural mortality. This is the
    # RECOMMENDED driver (decision 2026-08-20): its composition is identical to the
    # post-1985 Tupek product (tree litter incl. residues and mortality, understorey
    # EXCLUDED), so sigma_input applies the same understorey correction on both sides
    # of 1985. Including ground vegetation would make sigma_input mean different things
    # before and after the join -- a discontinuity no diagnostic would reveal.
    # It also removes the mixed-basis problem: ground vegetation is the only
    # upland-restricted term in Fig. 5, the rest cover all forest land.
    tree = [sum(cols[k][i] for k in ("tree_litter", "harvest_residues", "natural_mortality")
                if k in cols) for i in range(len(years))]
    allv = [sum(cols[k][i] for k in ("tree_litter", "harvest_residues",
                                     "ground_vegetation", "natural_mortality")
                if k in cols) for i in range(len(years))]
    path = os.path.join(outdir, "liski2006_fig5_input_to_soil.csv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["year"] + sorted(cols)
                   + ["total_input_tree_basis_TgC_yr", "total_input_all_TgC_yr"])
        for i, y in enumerate(years):
            w.writerow([y] + [round(cols[k][i], 4) for k in sorted(cols)]
                       + [round(tree[i], 4), round(allv[i], 4)])
    print("wrote", path)
    os.remove(svg)
    return years, cols, tree

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
