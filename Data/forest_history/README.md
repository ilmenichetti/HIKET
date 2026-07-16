# Finnish forest history — NFI growing stock

`nfi_growing_stock.csv` — total growing stock volume (million m³) on productive +
poorly productive forest land, whole Finland, National Forest Inventories NFI1–NFI13.

## Source
Korhonen K.T., Räty M., Haakana H., Heikkinen J., Hotanen J.-P., Kuronen M.,
Pitkänen J. (2024). *Forests of Finland 2019–2023 and their development 1921–2023.*
Silva Fennica 58(5) art. 24045. https://doi.org/10.14214/sf.24045 (open access).
PDF archived in `literature/Korhonen_etal_2024_ForestsOfFinland_SilvaFennica_24045.pdf`.

## Provenance of each point (see the CSV `provenance` column)
The paper reports the series **only as Figure 10a** (a plot) plus exact endpoints in
the text/tables — it does not tabulate the per-inventory series. Therefore:

- **Exact, from the paper:** NFI1 = 1400 (text: "1.4 G m³", recalculated with modern
  volume functions), NFI11 = 2356, NFI12 = 2475 (= NFI13 − the reported 77 M m³
  NFI12→NFI13 gain), NFI13 = 2552 (Table 55 / text "2.6 G m³", SE 13).
- **Digitized from Fig 10a (Total curve), ±~50 M m³:** NFI2–NFI10. The load-bearing
  feature is the **shape** — a depleted, near-stationary base (~1400–1500) into the
  ~1970s, then a sustained rise; the paper states "most of the increase … has taken
  place after the end of the 1960s." Read from the paper's own recalculated curve so
  the series is internally consistent with the recalculated NFI1 anchor.

## Consumers (single source of truth)
- `manuscript/figures/build_F10b_forest_history.R` — the F10b history figure.
- The transient pre-run in the six model wrappers (C3): the 1917→1985 litter-input
  *shape* is derived from this series (endpoints `J_1917`/`J_1985` stay calibrated;
  only the shape between them follows growing stock). Only NFI1–NFI7 (+interp to 1985)
  are used by the pre-run window.

The CSV is checked into git (small, citable); the source PDF is gitignored
(`literature/*.pdf`, copyright) but catalogued in `literature/README.md`.
