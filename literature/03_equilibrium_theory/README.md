# 03 — The theory of equilibrium initialisation

**Serves:** `HIKET_storyline_v3.tex` §"What steady state assumes" (`:201`) — the theory paragraph
requested in Lorenzo's annotation (p5): *what an equilibrium state is for a linear donor-controlled
system, what it assumes about the constancy of inputs and climate, what properties follow — C = J/k,
the loss of all timing information, and the fact that the assumption is about the **history**, not
about the soil.*

**This is not a literature search so much as a citation list for things we already know.** The
paragraph should be short, exact, and carry four or five canonical references. Below: the maths as it
should appear, the reference for each step, and one derivation that connects the theory directly to
our own central result.

---

## A. The statement, with its references

Every model in this study is a **linear, donor-controlled compartmental system**:

$$\frac{d\mathbf{C}}{dt} = \mathbf{b}\,J(t) \;-\; \mathbf{A}\,\xi(t)\,\mathbf{C}$$

with `C` the pool vector, `J(t)` the litter input, `b` its partitioning, `A` the (constant) transfer
and decay matrix and `ξ(t)` the climate modifier. Donor-controlled means each flux depends only on
the donor pool, which is what makes the system linear and gives it a unique, globally attracting
fixed point whenever `A` is invertible and the forcing is constant.

| Step | Statement | Reference |
|---|---|---|
| The fixed point | With `J(t) ≡ J̄` and `ξ(t) ≡ ξ̄`, the unique steady state is `C* = (A ξ̄)⁻¹ b J̄`. For one pool this is the familiar **C\* = J/k**. | Olson 1963, *Ecology* 44:322–331 — the origin of the litter `k` and of `X = L/k` at steady state |
| The identity behind it | `C* = J̄ × τ`, where `τ` is the mean residence time of the system. **Equilibrium determines only the product**, never its factors. | Ågren & Bosatta 1998; Sierra et al. 2017 on the vocabulary (age / turnover / transit / residence) |
| What is being assumed | Not a property of the soil but of its **history**: that `J` and `ξ` have been effectively constant over a window long compared with the slowest pool's residence time, and that no disturbance has occurred inside it. | Wutzler & Reichstein 2007 (folder 02) states the consequence when this fails |
| What is discarded | The stationary state is a *function of the parameters alone*. All information about trajectory — when input changed, how fast the soil has been responding — is discarded, and with it the age structure of the carbon. | Sierra, Müller, Metzler, Manzoni & Trumbore 2017; Metzler, Müller & Sierra 2018 (PNAS) for the age/transit-time formalism |
| The decomposition | Transient storage = **storage capacity** (the moving equilibrium of current forcing) + **storage potential** (the distance still to travel). An equilibrium start sets the second term to **zero by construction**. | Luo et al. 2017, *Biogeosciences* 14:145–161 |
| Why it cannot show a levelling-off | Zero potential means zero approach; the model can only track the moving capacity. A rise that decelerates is a statement about the potential term. | follows from the above — **our own argument, present it as such** |

---

## B. ⚠ The derivation that matters most for this paper

The equilibrium condition gives one equation in two unknowns:

$$C^* = J\,\tau \qquad\Longrightarrow\qquad \log C^* = \log J + \log \tau$$

so a calibration that is informed **only by the stock level** can identify `log J + log τ` and nothing
about `log J − log τ`. That is not an inconvenience of our sampler — **it is the ridge, derived from
the model definition**, and it is the same object as the measured `MRT × σ_input` trade-off
(CLAUDE.md; `F14_mrt_ridge`, `S13_mrt_ridge_benchmark`).

It also gives the storyline its corrected claim, exactly (and this replaces the wrong version that
was removed on 2026-09-03):

> The ridge is a property of the **equilibrium condition**, not of a particular initialisation. A
> transient start does not remove it; it **narrows** it, because the *approach* to the attractor
> depends on `τ` alone — the rate constant sets how fast the gap closes, while `J` sets only where it
> closes to. Observing the trajectory therefore adds information about `τ` that the level cannot
> contain. Measured here, the narrowing is partial: MRT and `σ_input` still correlate −0.64 to −0.73.

That is a clean, defensible, *theoretically grounded* version of the paper's central methodological
claim, and it needs no citation beyond the compartmental-systems literature.

---

## C. ⬜ A small computation worth doing (cheap, and it fills a real gap)

The slowest eigenvalue of `Aξ̄` sets the e-folding time of the approach to equilibrium:
`t_f = −ln(1−f)/λ_min` to close a fraction `f` of the gap. **Compute λ_min per model at the posterior
median, at the dataset-mean climate** — the same reference state `doublechecks/intrinsic_mrt.R`
already builds, so this is a few lines on top of it.

Why it is worth it:

- It converts "soils lag" from a literature claim into **our own number**, per model, and it is the
  quantity folder 04 §C is currently missing.
- It should reconcile the tension in that folder: bulk MRT is 22–27 yr, but λ_min is set by the
  **humus** pool and will be far longer — which is why a century of regrowth has not been enough.
- It is a prediction the reader can check against the fitted `σ_init`, and it is exactly the kind of
  thing a reviewer will ask for when the paper claims an equilibrium start is inadequate.

⚠ Report it as an **e-folding time of the slowest mode**, not as "the time to equilibrium" — a linear
system approaches asymptotically and never arrives, which is itself part of the argument.

---

## D. Draft text — §"What steady state assumes"

> All six models here are linear donor-controlled compartmental systems, in which each carbon flux
> depends only on the pool it leaves. Under constant litter input and constant climate such a system
> has a single, globally attracting fixed point, obtained by setting the derivative to zero: in the
> one-pool case the classical `C* = J/k` (Olson, 1963), and in general a stock equal to the input
> multiplied by the mean residence time of the system. Initialising a simulation at that fixed point
> is therefore not an assumption about the soil. It is an assumption about the soil's **history** —
> that input and climate have been effectively constant over a period long relative to the residence
> time of the slowest pool, and that no disturbance has intervened within it.
>
> Two properties follow, and both matter here. The first is that the equilibrium stock constrains only
> the *product* of input and residence time; the two are exchangeable at fixed stock, so a level
> alone cannot separate a fast soil receiving much litter from a slow soil receiving little. The
> second is that the stationary state is a function of the parameters alone, so all information about
> the trajectory is discarded: an equilibrium start has, by construction, no distance left to travel
> (Luo et al., 2017). It can track a moving target as climate and litter change, but it cannot
> express a soil that is still catching up — and therefore cannot produce an accumulation that
> decelerates, which is the dynamic the inventory most needs to project.

---

## E. What is still missing

1. ⬜ **Check Luo et al. 2017's exact vocabulary** ("carbon storage capacity" / "carbon storage
   potential") before adopting it — it is a good frame, but only if used as they define it.
   OA: `https://www.biogeosciences.net/14/145/2017/`.
2. ⬜ **Decide whether the derivation in §B goes in the Introduction or the Methods.** It is the
   theoretical statement of the paper's own result, so there is a case for stating it early and
   referring back — but it must not pre-empt the Results.
3. ⬜ The λ_min computation (§C).
4. ⬜ Confirm Ågren & Bosatta (1998, *Theoretical Ecosystem Ecology*, CUP) is the right general
   citation, or replace with Manzoni & Porporato 2009 (*SBB* 41:1355–1379).

## D. Search 2026-09-21 — who has put a NUMBER on the time to equilibrium, and for which models

Prompted by Lorenzo's Introduction note ("SOC is slow to react … how much time?"). Decided: no own
computation (§C stays optional); the sentence rests on the pool-specific timescale with the sources below.

| source | models | what it gives | verdict |
|---|---|---|---|
| **Wutzler & Reichstein 2007**, *Biogeosciences* 4:125–136 (open; PDF in folder 02) | Yasso (5 pools) | Defines **t95 = −ln(0.05)/k** for the slowest pool (their Eq. 6). Standard Yasso `k_hum2 = 1.2e-3 /yr` ⇒ t95 ≈ **2 500 yr**; with the relaxed rates needed to fit German chronosequences *"the times to reach equilibrium could span millennia"*. Also: calibrating to current stocks under an equilibrium assumption *overestimates the decay rate of the slowest pool*. | ✅ the citation for "centuries"; already in the bib |
| **⭐ Peltoniemi, Thürig, Ogle, Palosuo, Schrumpf, Wutzler, … Liski, Smith & Mäkipää 2007**, *Silva Fennica* 41(3):575–602, doi 10.14214/sf.290 (open; PDF in folder 02) | **7 models** reviewed: Yasso, RothC, Century, ROMUL, SOILN, Forest-DNDC, (+ a statistical class) | §2.2.5 *Model Initialization*: *"The assumption of a soil being in a steady state equilibrium with respect to current inputs is likely to be violated in most applications"*; and — our thesis, in 2007 — *"The soil C pool of a calibration site(s) may be far from equilibrium, but still the model parameters are calibrated so that the modelled equilibrium state matches the measured C pool. Therefore, a further correction of the parameters and the pools is necessary to avoid the underestimation of century-term soil C accumulation."* Table 2: Century *"typically initialized with spin-up (several 1000 years)"*; Yasso: *"initialization dominated the uncertainty of soil C balance but the use of longer spin-up period clearly reduced its effect (de Wit 2006, Peltoniemi 2006)"*. | ⭐⭐ **NOT YET CITED ANYWHERE.** Multi-model, Finnish, Mikko first author with Liski, Palosuo, Mäkipää and Wutzler. Cite in the Introduction (§"Steady state as a starting condition") as the group's own earlier statement of the problem — lineage, not correction. |
| **Xia, Luo, Wang, Weng & Hararuk 2012**, *GMD* 5:1259–1271 (open) | CABLE (CASA-CNP, 9 pools) | Traditional spin-up runs *"thousands of simulation years"*; the **passive SOM pool sets the spin-up time**; a semi-analytical steady state cuts it 92–96%. | ✅ the ESM-side statement that the slowest pool sets the timescale |
| **Exbrayat, Pitman & Abramowitz 2014**, *GMD* 7:2683–2692 (open) | CMIP5 ESMs | The 6-fold range in present-day soil C *"already exists at the beginning of the historical simulations"* and is set during spin-up; it persists to 2100 almost unchanged. | ✅ initialisation as a source of persistent between-model spread (the multi-model analogue of our F4) |
| **Wieder, Hartman, Sulman, Wang, Koven & Bonan 2018**, *GCB* 24:1563–1579 | CASA-CNP, MIMICS, CORPSE on a common testbed | Same initial stocks (~1 400 Pg) from common forcing, different turnover–temperature relations ⇒ divergent 20th-century trajectories (±20 Pg). | ✅ multi-model; spin-up equalises the LEVEL, not the RESPONSE — same pattern as our six |
| **Sierra, Hoyt, He & Trumbore 2018**, *GBC* 32:1574–1588 | "a wide range of SOC models" | Age vs transit-time distributions; ages centuries–millennia, transit times years–decades, for the same model. The formal reason a stock equilibrates on the AGE timescale while the output flux turns over on the TRANSIT timescale. | ✅ companion to Sierra 2017 (already cited); use if a reviewer asks why MTT 20–30 yr and "centuries" are not in conflict |
| Foereid, Bellamy, Holden & Kirk 2012, *Eur. J. Soil Sci.* 63:32–41 | RothC | Initialisation method changes predicted SOC change for England & Wales. | single model; optional |
| Nemo et al. 2017, *Environ. Model. Assess.* 22:215–229; Dimassi et al. 2018, *Geoderma* 311:25–36 | RothC; Century | Initialisation-scenario sensitivity on long-term experiments (cropland). | single model, cropland; optional |

**Resolution of the folder-04 tension** (60-yr post-fire recovery vs "centuries"): Sierra 2018 gives the formal
statement — the fast pools carry the transit time (years–decades), the slow pool carries the age (centuries),
and the stock follows the age. Nothing to compute.
