# HIKET — release notes

**What this file is for.** Everything that has to be said, done or checked when the
code and data are **packaged for publication** — repository release, Zenodo deposit,
data availability statement, and anything a reader of the paper would need in order to
run or understand the code. Collected here as it comes up, so that the packaging step
is not a memory exercise at the end.

⚠ **Not a changelog.** Scientific decisions, corrections and run history live in
`CLAUDE.md`. Put an item here only if it is something the *release* must carry.

---

## Naming

- **`MRT` in the code means the mean transit time (`MTT`) of the manuscript.** The
  manuscript follows the terminology of Sierra et al. (2017, *Glob. Change Biol.*
  23:1763–1773), who show that "residence time" has been used for at least three
  different quantities and discourage the term. The quantity we compute is unchanged
  and is their Eqn (2), `MTT = 1'(-A)^-1 b` with `b` the normalised input vector —
  equivalently their Eqn (4), the steady-state stock divided by the input flux.
  The identifiers `intrinsic_mrt.R`, `intrinsic_mrt.rds`, the `mrt` columns and the
  `build_*_mrt_*.R` filenames were deliberately **not** renamed, because the figure
  staleness guards key on those names. ⬜ **On release: say this in the code README**,
  or rename the identifiers in one sweep and re-verify the guards.

## Software environment

- R 4.3.1; `BayesianTools` 0.1.8 (DEzs sampler); `parallel`, `compiler`, `coda`.
- Yasso Fortran (`yasso07.f90`, `yasso15.f90`) is compiled locally and the `.so` files
  are gitignored. ⬜ **On release: ship the build instructions**, including the rule
  that each `.f90` gets its OWN `R CMD SHLIB` call — a combined call silently produces
  one `.so` and breaks Yasso15 and Yasso20.

## Data availability

- Litter input: Tupek et al., Zenodo doi 10.5281/zenodo.19736499 (public).
- **SOC baseline: not a published dataset, and deliberately NOT cited as one.** An
  unpublished internal compilation is not a reference, and citing it would dress our own
  data up as a source (decided with Lorenzo, 2026-09-11). The Methods therefore describe
  what the data are — the layer-wise measurements of the three campaigns, put on one set
  of conventions for this study — and name no source. ⬜ **On release: the data
  availability statement and the acknowledgements carry the provenance instead.** Settle
  with the people who hold the measurements what can be deposited and how they wish to be
  credited. A `\gap` marker in the draft tracks this.
- Climate: `nfi_plot_weather_data_1961_2025.nc`. **Provisionally** FMI's 10 km gridded
  daily climatology (Aalto, Pirinen & Jylhä 2016, doi 10.1002/2015JD024651), extracted to
  plot coordinates by nearest neighbour. ⚠ **Not confirmed**, and the file we hold is a
  *modified* version whose modification is undescribed; the NetCDF records only
  `CRS: YKJ-KKJ` and the extraction method. ⬜ Confirm the product and describe the
  modification before submission, and settle whether the extracted file can be
  redistributed.

## Reproducing the sensitivity analysis on the input prior (added 2026-09-23)

- The σ_input prior centre is set by `HIKET_SIGMA_INPUT_CENTRE` in all six
  `Prior_specs/*_priors.R` (unset = **1.08**, the production value, Lehtonen & Heikkinen
  2015). Setting it to **1.27** (Liski et al. 2006) reproduces the sensitivity calibration
  reported in the paper; nothing else differs between the two. On the CSC cluster the
  variable must be passed as `SINGULARITYENV_HIKET_SIGMA_INPUT_CENTRE` to reach R inside the
  container. Before this switch the value was edited by hand, which the released code could
  not reproduce.
- The Yasso matrix exponential (`matrixexp` in the Fortran) is FMI's Yasso15 core code
  (Järvenpää) with two changes of ours: 20 instead of 10 Taylor terms, and a cap on the
  number of scaling steps. The method is standard scaling and squaring (Moler & Van Loan 2003).
