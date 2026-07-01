# doublechecks/

Scratch working folder for sanity-checking the HIKET pipeline. Nothing here
feeds the production pipeline — these are stand-alone diagnostic scripts plus
the figures they produce (`figures/`).

**Convention:** every script is self-contained and is run from the project
root, e.g.

```bash
Rscript doublechecks/check_litter_inputs.R
```

| Script | What it checks |
|--------|----------------|
| `check_litter_inputs.R` | The litter forcing each model integrates is one and the same object, and matches the raw Zenodo source CSV (no aggregation/units artefact). Plots the cross-plot mean litter trajectory from the source directly vs. from each model's input bundle. |
