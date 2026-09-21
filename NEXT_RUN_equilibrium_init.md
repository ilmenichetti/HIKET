# The equilibrium-init counterfactual — run design (drafted 2026-09-04)

**Why.** The Introduction's load-bearing claim — *an equilibrium start cannot produce a rise, so it
cannot show a rise levelling off* — is currently a **theoretical** result from `C* = Jτ`. No
experiment in `doublechecks/` tests it: there are ablations for campaign weighting, pre-run shape,
priors and forward arms, but **no equilibrium-initialised arm**. A reviewer will ask for it.

**Why a forward-only version was rejected** (Lorenzo, 2026-09-04). Swapping the init under posteriors
that were *calibrated* with a transient init is a rigged comparison: the equilibrium arm fails
because it was never allowed to re-fit, not because equilibrium initialisation is inadequate. The
test only means something if the equilibrium arm gets the same freedom — `σ_input`, `σ_init`'s
replacement, fractions, climate — to do its best.

## Design

**One factor: the initial state. Nothing else moves.**

| | transient arm (exists) | equilibrium arm (new) |
|---|---|---|
| initial state | 1917 pre-run + 68-yr ramp, `σ_init` free | steady state at each draw's parameters and the **1985** litter flux |
| `σ_init` | free | **removed** — it has no referent |
| `σ_input` | free, arm-B prior | free, **identical prior** |
| everything else | arm B | byte-identical |
| target, likelihood, σ | corrected target, correlated likelihood | identical |

⚠ **`σ_init` cannot simply be fixed to 1** — under equilibrium init there is no 1917 state to scale,
so the parameter must be *dropped* and the free set becomes one shorter. That is a genuine asymmetry
in degrees of freedom and it **favours the transient arm**; say so, and check whether the fit
difference survives an information criterion that charges for it.

⚠ Which flux defines the equilibrium is itself a choice — 1985 alone, or a 1960–1990 mean as the NID
does. **Use the NID's own convention** so the arm is the convention as practised, not a straw man.

## What it should show, pre-registered

1. **The level: both arms reach it.** `σ_input` scales `C*` directly, so the equilibrium arm can buy
   the 1985 level. **If it cannot, something is wrong with the arm, not with the convention.**
2. **The trajectory: only the transient arm bends.** The equilibrium arm should be near-flat, tracking
   only the moving capacity as climate and litter change.
3. **⭐ The number to watch: the equilibrium arm's 1985–2024 rate.** Liski et al. 2006 started at
   equilibrium in 1922 and got ~**+0.024** tC ha⁻¹ yr⁻¹ against our observed **+0.259**. If our
   equilibrium arm lands near theirs, the Introduction's antecedent and the Results close a loop and
   the paper has its single most persuasive figure.
4. **`σ_input` should go UP in the equilibrium arm**, because the only way to raise a stock with no
   disequilibrium left is more input. That would also make the ridge argument concrete.

**Falsifier:** if the equilibrium arm reproduces the observed rate within its interval, the paper's
central claim is wrong and must be restated. Pre-register that.

## Cost and sequencing

Six calibration jobs, same footprint as any production run. ⚠ **Do not bundle it** with the
`σ_input` anchor correction (`literature/UNDERSTOREY_ANCHOR.md`, 1.08 → ~1.21) — that is a second
factor and the project's one-factor discipline exists for exactly this reason. Two runs, in this
order:

1. **Anchor correction** (σ_input centre 1.08 → ~1.21, arm-A-equivalent), because it changes what the
   input finding *is* and therefore what §missingflux says.
2. **Equilibrium-init counterfactual**, on whichever anchor is then current.

⚠ The pending **Fortran guard recompile** (0.9999 → 1.0) must land first; the point of that run was
the `.so` rebuild — one `R CMD SHLIB` per `.f90` on Roihu.

## New figure

**F16 — the counterfactual.** Two model ensembles over 1985–2024 against the three campaign means:
transient arm rising and bending, equilibrium arm flat, observed markers with CI. One panel. It is
the paper's thesis in a single picture and it is currently missing.
