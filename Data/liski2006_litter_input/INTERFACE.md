# Interface — how to swap in a better input series

The pre-run shape is built by `Model_functions_real_data_transient/preinit_input_shape.R` and
consumed identically by all six `run_*_transient_calibration.R`:

    PREINIT_SHAPE <- growing_stock_preinit_shape()   # 1917->1985, 68 steps, [0,1]

**Contract: a length-68 numeric vector on [0,1], 0 at 1917 and 1 at 1985.** Anything satisfying that
drops in without touching the calibration scripts.

## To plug in a replacement series

1. Put it in this folder as a CSV with a `year` column and a total-input column in any consistent
   flux unit (the shape is normalised, so absolute units cancel).
2. Divide by a **per-hectare** basis before normalising — either the series is already per hectare,
   or use `liski2006_upland_soil_area.csv` (13.90 → 15.10 M ha over 1917–1985, growth placed in
   1965–1980). ⚠ Skipping this imports the forest-area expansion, which is the defect this whole
   exercise exists to remove.
3. Keep the composition matching the post-1985 driver — tree litter incl. residues and mortality,
   understorey excluded. A pre/post mismatch silently changes what `sigma_input` means at the join.
3. Interpolate to 1917–1985 annual and normalise to [0,1].

## Decisions in force (2026-08-20)

1. **Use `total_input_tree_basis_TgC_yr`**, not `total_input_all_TgC_yr` — its composition matches
   the post-1985 Tupek product, so `sigma_input` means the same thing on both sides of 1985.
   Ground vegetation is excluded deliberately; it is absorbed by the sigmas, as it is after 1985.
2. **Area timing 1965–1980** (Korhonen-sourced), with "no correction" carried as the reported bound.
   ⚠ Constant area is not the neutral choice — it is a claim Liski's own numbers falsify.

## Status

- ✅ data extracted, validated, documented (`README.md`)
- ✅ area basis recorded (`liski2006_upland_soil_area.csv`)
- ⬜ **NOT WIRED** — `preinit_input_shape.R` still derives the shape from growing stock alone, and
  nothing reads these files. Wiring means adding a second shape function plus a switch
  (`HIKET_PREINIT_LINEAR=1` already exists as the precedent for such a switch).
- ⬜ the area *timing* (1965–1980) is an assumption from Korhonen et al. 2024, not data.

## If the real series arrives from A. Lehtonen

Replace `liski2006_fig5_input_to_soil.csv` with it, keeping the same column contract, and note the
provenance change in `README.md`. The extraction here then becomes an independent cross-check on his
numbers rather than the source — worth keeping either way, since it validated to <5% on four tests.
