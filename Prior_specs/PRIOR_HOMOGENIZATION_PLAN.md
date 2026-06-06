# Prior specification & homogenization — design note

> **STATUS: IMPLEMENTED in `Prior_specs/*_priors.R` (2026-06-06).**
> All six prior files now carry the locked scheme: Yasso07 climate/woody widths
> on the published 1σ scale, all 12 fractions at logit SD 0.4 (decision #5:
> **explicit per-fraction listing**), SP1/TP2/TP3 climate widths on the Yasso07
> scale, TP2/TP3 fractions at logit 0.4. Verified: a prior-predictive draw shows
> the new Yasso07 `beta2` collapses the `xi` overflow rate from ~33 % to 0 %
> (p99 of the temperature factor 3e31 → 8.9). Locally de-risked at full N for all
> six models via `preflight_prior_pushforward.R` (§8); methods-doc table done.
> **Still pending:** one Puhti production run (all 6, full N, full iter) to
> confirm MCMC convergence behaviour — climate R-hat should improve while
> fraction R-hat stays poor (§4.2). Decision record (all five resolved) in §7.
>
> Width convention (decision #1): published "±" limits are read **conservatively
> as 1σ** (NOT divided by 1.96) → broader, less-informative priors. Centres
> (decision #2): **GUI centres kept**, published widths adopted. Woody widths
> (decision #3): from **Tuomi 2011 Table 4** (paper now on disk, see §7b).
> Fraction logit SD (decision #4): **0.4** (narrower preferred).

---

## 1. Why this note exists (the trigger)

The Yasso07 transient calibration of run `20260603_053323` failed to converge
(univariate R-hat 1.8–38.9, multivariate psrf 162, ESS 1–13, 25–46 % of MCMC
proposals returning `-Inf`). Root cause: the prior SD on `beta2` was set to
**0.05**, copied from the hand-set "weakly informative" block used by
SP1/TP2/TP3, *not* from any Yasso07 calibration.

`beta2` multiplies **T²** in the climate modifier

```
k_i(C) = alpha_i * exp(beta1*T + beta2*T^2) * (1 - exp(gamma*P_a))      (Tuomi 2009, Eq. 3)
```

In boreal climate the seasonal temperature quadrature points reach |T| ≈ 25 °C,
so T² ≈ 600. With `beta2` free to roam to ±0.1 (the 95 % span of an SD-0.05
prior), `beta2·T² ≈ ±60` and `exp(±60)` spans ~50 orders of magnitude. The
climate modifier `xi` therefore explodes or collapses, the transient pre-run
(1917→1985) compounds it through three integration stages, predicted SOC goes
to 0 or Inf on some plots, the all-or-nothing likelihood
(`calibration_engine_transient.R`: any one plot `-Inf` ⇒ whole proposal `-Inf`)
rejects ~half of proposals, DEzs starves, and the chains cannot mix.

Yasso15/20 do **not** have this problem because their `beta2` widths are
empirical posterior SDs (0.00014 / 0.00054) — ~100–360× tighter — so
`beta2·T²` stays bounded. Their R-hat is "bad but not crazy" (all parameters
in Warning, worst ≈ 3–5, multivariate psrf ≈ 10) for an unrelated reason:
genuine **structural non-identifiability** of ~19 free parameters against only
2 SOC observations per plot. That is a scientific finding, not a bug, and is
deliberately *not* fixed by tightening priors.

---

## 2. The homogenization concept (the important part)

The HIKET scientific question is whether divergent SOC projections across
models reflect **genuine structural differences** or **calibration artefacts**.
To answer it, the calibration setup must not itself inject asymmetry that
confounds structure with calibration choices. The Yasso07 `beta2=0.05` bug is
exactly such an artefact: one model's climate prior was 100× looser than the
others' for no defensible reason.

**Key principle: priors are homogeneous in *construction method, scale,
constraints, and degrees of freedom* — NOT in identical raw numbers.**

Identical raw numbers would be *wrong*, because the parameters live on
model-specific scales (different published MAPs, different identified
uncertainties). Forcing the same number is literally what broke Yasso07.
Homogeneity is therefore defined at the level of *how each prior is built*, so
that any residual difference in posterior behaviour or projection is
attributable to model **structure** and **data**, never to an arbitrary prior
asymmetry.

The three Yasso models (07, 15, 20) share the same structural skeleton — five
pools (A, W, E, N, H), the same compartmental ODE, the same climate-response
*functional form*, the same woody-size modifier, the same stick-break simplex
on inter-pool transfers — so homogeneity is meaningful and required.

### Three-tier scheme

| Tier | Parameters | Centre | Width | Homogeneous in… |
|---|---|---|---|---|
| **1 — climate & size** | `beta1/2`, `gamma`, `betaN1/2`, `gammaN`, `betaH1/2`, `gammaH`, `delta1/2`, `r` | each model's own published MAP | each model's own published posterior SD (published "±" limit read as 1σ; see §6) | **method + scale** (per-model numbers, same provenance rule, same order of magnitude) |
| **2 — transfer fractions** | the 12 `p_*` | each model's own published structure | **identical**: logit SD = **0.4** (was 1.0) | **method + scale + number** (width literally identical) |
| **3 — auxiliary uncertainty** | `sigma_init`, `sigma_input` | 1.0 | **identical**: log SD = 0.50 | **method + scale + number** |

### What is homogeneous (must be)
- The **rule** for every prior: Tier-1 centre & width from that model's own
  calibration; Tier-2 a common weak logit prior; Tier-3 a common log prior.
- The **scale** of Tier-1 widths: all on the empirical order of magnitude —
  no outliers like the old `beta2=0.05`.
- The **transforms / constraints**: stick-break simplex for fractions
  (budget `B = 1 − p_H`), `log` for strictly-positive params
  (`beta*1`, `delta2`, `r`, `sigma_*`), `unconstrained` for sign-free params
  (`beta*2`, `gamma*`, `delta1`). Bit-identical across the three models.
- The **degrees of freedom**: all 12 transfer fractions free in every model
  (the Yasso20 structural fix); same auxiliary parameters everywhere.

### What legitimately differs (and should — it is structure/data, not asymmetry)
- **Tier-1 centres and widths**, because each comes from a different published
  calibration on a different dataset.
- **The parameter set itself**: Yasso07 has a single AWE climate response and
  no separate N/H climate parameters; Yasso15/20 add pool-specific climate
  (`betaN*`, `betaH*`, `gammaN`, `gammaH`). That is a genuine structural
  difference, not a calibration choice.
- **Decomposition rate constants** (`alpha_*`), fixed per model.

### The tier split is driven by *consistency of the original calibrations*

The tier of a parameter follows its role in the **HIKET** calibration, not what
the original authors did — so the same class sits in the same tier across all
three models (`beta2` is Tier 1 everywhere; all 12 fractions are Tier 2
everywhere). But *why* the line falls between climate/size and fractions comes
directly from what each original study left free:

| Parameter class | Yasso07 (Tuomi 2009) | Yasso15 | Yasso20 (Viskari 2022) | Consistent? | Tier |
|---|---|---|---|---|---|
| Climate AWE (`beta1/2`, `gamma`) | free | free | free | yes | 1 |
| Climate N/H (`betaN*`, `betaH*`, `gammaN/H`) | not in model | free | free | yes (where present) | 1 |
| Size (`delta1/2`, `r`) | free (2011 woody) | free | free | yes | 1 |
| 12 transfer fractions | free (all 12) | free (all 12) | **only 3 free; 6 fixed at 0, 3 algebraic** | **NO** | 2 |
| `sigma_init/input` | HIKET-only | HIKET-only | HIKET-only | — | 3 |

- **Climate + size were freely calibrated in every original model** → their
  published posteriors have consistent coverage → the "own MAP + own posterior
  SD" recipe applies uniformly → **Tier 1**.
- **Fractions were calibrated inconsistently** — Viskari (2022) fixed/derived 9
  of Yasso20's 12. For those 9 there is *no usable published width* (a frozen
  parameter has zero posterior spread). Applying the Tier-1 recipe would freeze
  Yasso20's 9 fractions (SD≈0, re-imposing the warm-site structure) while
  Yasso15's same fractions stay free — a model-to-model asymmetry, exactly the
  calibration artefact HIKET exists to avoid. So fractions **deliberately do
  not use published posteriors at all** (not even Yasso07/15's, where they
  exist); they get a common weak prior instead → **Tier 2**.

In one line: *Tier 1 uses published posteriors because coverage is consistent
across models; Tier 2 overrides them because Yasso20 kept most fractions fixed,
and a common prior restores symmetry (equal degrees of freedom).*

### Why fractions get an *identical* prior but climate does not
For the climate/size parameters we *have* trustworthy per-model posteriors and
want each model to carry its own identified uncertainty. For the transfer
fractions we deliberately do **not** want any model's published structure to
dominate — six of Yasso20's fractions were structural zeros in Viskari (2022),
and the whole point of HIKET is to let Finnish data, not the original warm-site
constraints, decide the flow architecture. So all fractions get the *same*
weakly-informative prior, centred on published structure but wide enough to be
data-driven, identical across models for symmetry. The logit scale gives every
fraction the same *relative* exploration window regardless of whether its
centre is 6e-4 or 0.5 — so a previously-fixed near-zero fraction explores a
nearby region (~×2.7 of its centre at ±2σ with SD 0.5), not a different order
of magnitude (~×7 with the old SD 1.0). This is the "narrow a little, stay
honest" middle: fraction R-hat will remain poor (the non-identifiability
finding survives), but the chains stop wandering into pathological,
forward-model-exploding regions.

---

## 3. The climate response and what the widths must guarantee

```
xi = mean_over_seasons[ exp(beta1*T + beta2*T^2) ] * (1 - exp(gamma*P_a/1000))
```
(`*_wrapper_transient.R::compute_xi_*`; precipitation divided by 1000 in code.)

The prior on the T² coefficient (`beta2`, `betaN2`, `betaH2`) is the only one
that can detonate `xi`, because of the T² lever. The homogenization must keep
the 95 % prior span of `beta2·T²` at boreal |T|≈25 °C within ~±2, i.e.
`xi`'s temperature factor within roughly a 3–6× band — not 50 orders of
magnitude. The numbers in §4 satisfy this for all three models.

Illustration (AWE `beta2` at T² = 625; 95 % span = centre ± 1.96·SD):

| width source | `beta2` 95 % span | `beta2·T²` span | `exp(beta2·T²)` span |
|---|---|---|---|
| **current Yasso07 (SD 0.05)** | −0.0996 … +0.0964 | −62 … +60 | 1e-27 … 1e26  ← detonates |
| **proposed Yasso07 (SD 0.00065, 1σ)** | −0.0028 … −0.0003 | −1.78 … −0.19 | 0.17 … 0.83  ← bounded |
| Yasso15 (SD 0.00014) | −0.00049 … +0.00006 | −0.31 … +0.04 | 0.74 … 1.04 |
| Yasso20 (SD 0.00054) | −0.0031 … −0.0009 | −1.93 … −0.58 | 0.15 … 0.56 |

---

## 4. Full prior tables

Notation: `sigma_ppm` is the prior SD in **unconstrained (transformed)**
space. Physical 95 % range = inverse-transform of `centre_unconstrained ±
1.96·sigma_ppm`. For `log` params, `sigma_ppm` ≈ relative SD, and the physical
range is `centre · exp(±1.96·sigma_ppm)`. For `unconstrained` params,
physical = transformed (direct). For `stick_break` fractions, `sigma_ppm` is
the SD of `logit(p_i / remaining_stick)`.

### 4.1 Tier 1 — climate & size

**Yasso07** — centres = GUI / `YASSO07_FREE_DEFAULTS` (decision #2). Climate
widths from Tuomi 2009 Table 3; woody widths from Tuomi 2011 Table 4. **"±"
read as 1σ** (decision #1), so the width equals the published limit directly
(climate) or that limit ÷ centre for `log` params (relative SD, delta method).

| Param | transform | centre | **current SD** | **proposed SD** | derivation (1σ) |
|---|---|---|---|---|---|
| `beta1` | log | 0.09873 | 0.20 | **≈ 0.26** | T3: 7.6 ±2.0 ×10⁻²; log SD = 0.020 / 0.076 (rel. to paper MAP) |
| `beta2` | unconstrained | −0.001572 | **0.05** | **≈ 0.00065** | T3: −8.9 ±6.5 ×10⁻⁴; direct |
| `gamma` | unconstrained | −1.27169 | 0.30 | **≈ 0.20** | T3: −1.27 ±0.20; direct |
| `delta1` | unconstrained | −1.70841 | 0.15 | **≈ 0.16** | T4: −1.71 ±0.16; direct |
| `delta2` | log | 0.85856 | 0.10 | **≈ 0.12** | T4: 0.86 ±0.10; log SD = 0.10 / 0.86 |
| `r` | log | 0.30680 | 0.015 | **≈ 0.042** | T4: −0.306 ±0.013; log SD = 0.013 / 0.306 |

> Decisions applied here: **#1** widths read as 1σ (conservative/broader);
> **#2** GUI centres kept; **#3** woody widths from Tuomi 2011 Table 4 (paper
> on disk — see §7b). Provenance line for methods: *"centres from the Yasso07
> GUI (Tuomi et al. 2011, EMS); climate-prior widths from the published
> posterior limits (Tuomi et al. 2009 Table 3); woody-size-prior widths from
> Tuomi et al. 2011 (Ecol. Modelling) Table 4; limits read as 1σ."*
>
> Centre note (FYI, no action): the Table-3 *leaf-only* MAP (`beta1`=0.076,
> `beta2`=−0.00089) differs from the GUI centres; the GUI reflects the later
> leaf+wood combined fit. Woody Table-4 MAPs (`delta1/delta2/r`) **match** the
> GUI centres exactly, confirming the GUI woody params come from the 2011 paper.
> The current hand-set `r`=0.015 was *too tight* vs the published 0.042.

**Yasso15** (centres & widths = `Yasso15.dat` posterior; **NO CHANGE** — these
are already empirical and correctly scaled):

| Param | transform | centre | SD (current = proposed) |
|---|---|---|---|
| `beta1` | log | 0.09062 | 0.04593 |
| `beta2` | unconstrained | −0.000215 | 0.00014 |
| `gamma` | unconstrained | −1.80897 | 0.07022 |
| `betaN1` | log | 0.04878 | 0.10502 |
| `betaN2` | unconstrained | −0.0000792 | 0.00008 |
| `gammaN` | unconstrained | −1.17294 | 0.15543 |
| `betaH1` | log | 0.03518 | 0.13477 |
| `betaH2` | unconstrained | −0.000208 | 0.00016 |
| `gammaH` | unconstrained | −12.54094 | 1.50000 (capped; near-unidentifiable at Finnish precip) |
| `delta1` | unconstrained | −0.43883 | 0.32810 |
| `delta2` | log | 1.26838 | 0.28643 |
| `r` | log | 0.25687 | 0.05504 |

**Yasso20** (centres & widths = `Yasso20.dat` posterior; **NO CHANGE**):

| Param | transform | centre | SD (current = proposed) |
|---|---|---|---|
| `beta1` | log | 0.158 | 0.10355 |
| `beta2` | unconstrained | −0.002 | 0.00054 |
| `gamma` | unconstrained | −1.44 | 0.14583 |
| `betaN1` | log | 0.17 | 0.06246 |
| `betaN2` | unconstrained | −0.005 | 0.00041 |
| `gammaN` | unconstrained | −2.00 | 0.03532 |
| `betaH1` | log | 0.067 | 0.09118 |
| `betaH2` | unconstrained | 0.000 | 0.00009 |
| `gammaH` | unconstrained | −6.90 | 1.39056 |
| `delta1` | unconstrained | −2.55 | 0.43530 |
| `delta2` | log | 1.24 | 0.21098 |
| `r` | log | 0.25 | 0.05446 |

### 4.2 Tier 2 — transfer fractions (all three models)

Stick-break simplex, 4 columns (A/W/E/N), budget `B = 1 − p_H`. Centres are
each model's own published 12 fractions (Yasso07 from Tuomi 2009 / GUI;
Yasso15 & Yasso20 share the FMI Ryassofortran defaults). Width is **[TODO]**
changed from the implicit `rep(1.0)` default to a common **0.4** (decision #4 —
narrower preferred).

| | centre | logit SD (current) | logit SD (proposed) |
|---|---|---|---|
| all 12 `p_*` | per-model published | **1.0** (the `rep(1.0)` default) | **0.4** |

Relative exploration window (why narrower / 0.4):

| logit SD | mid-stick (centre 0.50), ±1σ | small fraction (centre 0.01), +2σ reach |
|---|---|---|
| 1.0 (current) | 0.27 – 0.73 | ~0.069 (≈ 7× centre) |
| 0.5 | 0.38 – 0.62 | ~0.027 (≈ 2.7× centre) |
| **0.4 (chosen)** | 0.40 – 0.60 | ~0.022 (≈ 2.2× centre) |

This is corroborated by Yasso07's own posterior: Tuomi 2009 Table 3 lists most
fractions as `0.00 (+0.07, −0.00)` — the data put them near 0 with a 95 %
*upper* limit ≈ 0.07, i.e. a small nearby region, not 0–0.5.

> **Guardrail — do not manufacture convergence.** The fractions are genuinely
> non-identifiable; 0.4 is meant to stop pathological roaming, *not* to force
> R-hat down. After the test run, **verify the fraction R-hat is still poor**
> (Warning-level). If it suddenly looks healthy, the prior is too tight and is
> imposing structure — back off toward 0.5. (User steer: "narrower better", but
> 0.3 is the practical floor.)
>
> Implementation note: because the fractions currently rely on the
> `sigma_ppm <- rep(1.0, N_FREE)` default in each run script, setting them to
> 0.4 means **either** adding the 12 `p_*` names to each `*_SIGMA_PPM` vector
> **or** changing the default fill from `1.0` to `0.4` and overriding only
> Tier-1/Tier-3. **Mechanism = decision #5, LEFT OPEN for next session**
> (recommendation: explicit per-fraction listing for traceability).

### 4.3 Tier 3 — auxiliary uncertainty (all models, no change)

| Param | transform | centre | SD |
|---|---|---|---|
| `sigma_init` | log | 1.00 | 0.50 |
| `sigma_input` | log | 1.00 | 0.50 |

`sigma_init` = proportional SD on the initial carbon state (enters the
transient pre-run physically); `sigma_input` = multiplier on litter inputs.
Both identical across all six models.

---

## 5. The simple models (SP1, TP2, TP3)

These share the *identical* AWE climate form `exp(beta1*T + beta2*T^2)` (see
`sp1_wrapper_transient.R`, `tp2_wrapper_transient.R`) and the *same*
`beta2=0.05` hand-set width — so they carry the **identical latent landmine**
and will explode the same way under transient init. They have no published
calibration of their own, and their climate centres were already borrowed from
the Yasso07 MAP.

**Plan [TODO]:** put their climate widths on the same Yasso07 empirical scale
(now the 1σ numbers) — `beta1 ≈ 0.26` (log), `beta2 ≈ 0.00065`, `gamma ≈ 0.20`
— matching their already-borrowed Yasso07 centres. Their rate/structural
parameters (`alpha_*`, `p_H`, `p_S`) keep their model-specific
weakly-informative widths, narrowed only if they prove to be harmful-wide.
SP1 has no transfer fractions; TP2/TP3 have `p_H`/`p_S` (logit) that should
follow the Tier-2 logit-0.4 rule. **Verify each model's `beta1` transform type
before applying** (assumed `log`, as in the Yasso models).

---

## 6. Conversion conventions (so numbers are reproducible)

1. **Published "±" limits read as 1σ** (decision #1, conservative/broader). The
   Tuomi 2009 Table 3 / 2011 Table 4 captions say "95 % confidence limits", but
   we deliberately treat the half-width as 1σ (no ÷1.96) to keep the climate
   priors broad and avoid over-claiming certainty on identifiable parameters.
   The fix holds regardless: even this broader `beta2`=0.00065 is ~77× tighter
   than the broken 0.05 and keeps `xi` bounded (§3).
2. **`log`-transform parameters** (`beta1`, `delta2`, `r`): physical SD → log
   SD via the delta method, `sigma_log ≈ sigma_phys / centre`. Width is
   therefore relative. For `beta1` the centre used in the ratio is the *paper
   MAP* (0.076), preserving the published relative uncertainty even though the
   prior is centred on the GUI value.
3. **`unconstrained` parameters** (`beta2`, `gamma`, `delta1`): transformed =
   physical, so `sigma_ppm` is the physical SD (= the 1σ limit) directly.

---

## 7. Decisions (resolved 2026-06-05)

1. **± convention** — ✅ **read as 1σ** (conservative / broader). No ÷1.96.
2. **Yasso07 centres** — ✅ **keep GUI centres**, adopt published widths.
3. **Yasso07 woody widths** (`delta1/delta2/r`) — ✅ **from Tuomi 2011 Table 4**
   (paper located on disk, §7b): `delta1`=0.16, `delta2`=0.12 (log),
   `r`=0.042 (log). No borrowing needed.
4. **Tier-2 fraction width** — ✅ **0.4** (narrower preferred; floor 0.3,
   guardrail in §4.2).
5. **Fraction width mechanism** — ✅ **explicit per-fraction listing** (decided
   2026-06-06). The 12 `p_*` are written out at logit 0.4 in each Yasso
   `*_SIGMA_PPM`; TP2/TP3 keep their already-explicit `p_H`/`p_S` entries.

---

## 7b. Local source assets (what's already on disk)

- **Tuomi 2009 paper (Table 3) — downloaded:**
  `Model_functions_real_data_transient/Decomposition_functions/Yasso_original/papers/1-s2.0-S030438000900386X-main.pdf`
  (Table 3 on p. 3367 = PDF page 6). Elsevier blocks automated re-fetch, so use
  this local copy; the extracted climate numbers are already in §4.1.
- **Tuomi 2011 woody paper (Table 4) — on disk:**
  `Model_functions_real_data_transient/Decomposition_functions/Yasso_original/papers/ecomodo2011.pdf`
  (Table 4 on p. 713 = PDF page 5). Source of the Yasso07 `delta1/delta2/r`
  widths (decision #3). doi:10.1016/j.ecolmodel.2010.10.025.
- **Yasso15/20 posterior samples:**
  `Model_functions_real_data/Decomposition_functions/Yasso_original/Yasso15.dat`,
  `Yasso20.dat`, `Yasso20_sample_parameters.rda`. `Priors_model_matching.R`
  already reads SDs from these (no Yasso07 equivalent exists).
- **Prior centres source (Yasso07):** `y07par_gui.csv` (single MAP row only —
  no uncertainty; that is why Table 3 is needed for widths).

## 8. Implementation checklist (when approved)

- [x] **Decision #5 settled** — explicit per-fraction listing (2026-06-06).
- [x] **`sigma_ppm` transformed-space convention honoured** — `log` params
      (`beta1`, `delta2`, `r`) given relative SDs (delta method, §6);
      `unconstrained` params (`beta2`, `gamma`, `delta1`) given the physical 1σ
      limit directly; fractions the stick-break-logit SD 0.4.
- [x] `Prior_specs/Yasso07_priors.R` — `beta1`→**0.26**, `beta2`→**0.00065**,
      `gamma`→**0.20**, `delta1`→**0.16**, `delta2`→**0.12**, `r`→**0.042**;
      12 `p_*` at logit **0.4**; header comment fixed (Tuomi 2009 T3 climate,
      Tuomi 2011 T4 woody, 1σ).
- [x] `Prior_specs/Yasso15_priors.R` — 12 `p_*` at logit **0.4** (Tier-1 unchanged).
- [x] `Prior_specs/Yasso20_priors.R` — 12 `p_*` at logit **0.4** (Tier-1 unchanged).
- [x] `Prior_specs/SP1_priors.R` — `beta1`→**0.26**, `beta2`→**0.00065**,
      `gamma`→**0.20**. `beta1` transform confirmed `log`.
- [x] `Prior_specs/TP2_priors.R` — climate widths to Yasso07 scale; `p_H` logit **0.4**.
- [x] `Prior_specs/TP3_priors.R` — climate widths to Yasso07 scale; `p_S`/`p_H` logit **0.4**.
- [x] All six files source cleanly; sigma names align 1:1 with `*_FREE_DEFAULTS`.
- [x] Prior-predictive verification: new `beta2` collapses `xi` overflow rate
      ~33 %→0 % (temperature-factor p99 3e31→8.9). The §3 detonation is defused.
- [x] **Local de-risk done (2026-06-06)** — `preflight_prior_pushforward.R`
      sources each model's real calibration script up to the MCMC launch and
      pushes prior draws through the genuine `ll_fn` at full N (no Puhti). The
      `beta2` detonation is gone in all six models; a `beta2` sweep shows the
      failure was prior-driven, not engine-driven (ll ≤ −10⁹ for positive
      excursions; −Inf only past +0.12). NEW-prior forward-blowup rates:
      SP1/TP2/Yasso15/Yasso20 0 %, Yasso07 ~2–3 %, TP3 ~10–12 % (genuine
      under-constraint: `alpha` ~5 %, climate ~3 %, fractions & `sigma` 0 % —
      left untouched per §5).
- [x] **Methods-doc table done** — `Calibration_real_data_transient/documentation/`
      `HIKET_calibration.Rmd` §"Priors" now carries the parameter-class ×
      prior-criteria table plus updated numeric tables (the file does exist; it
      lives in that `documentation/` subdir, not the repo root).
- [ ] **PENDING (Puhti):** one production run, all 6 models in parallel, full N,
      full iter — the only check local testing cannot give is the MCMC
      convergence behaviour: climate R-hat should *improve* while fraction R-hat
      *stays poor* (§4.2 guardrail — we are not forcing convergence).
- [x] Removed the stale `sigma_init` CI comment in
      `run_Yasso15_transient_calibration.R` (`[0.04, 0.27]`→`[0.37, 2.72]`).
