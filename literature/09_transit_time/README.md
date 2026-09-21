# 09 — Transit time, age, turnover: getting the terminology right

**Serves:** the Materials and methods subsection *"Mean transit time as a common measure of a
model"*, and the rename of that quantity throughout the manuscript (2026-09-11).

**Primary source, read in full:** Sierra, C. A., Müller, M., Metzler, H., Manzoni, S., Trumbore,
S. E. (2017), *The muddle of ages, turnover, transit, and residence times in the carbon cycle*,
**Global Change Biology 23:1763–1773**, doi 10.1111/gcb.13556. Local copy
`Sierra_etal_2017_muddle_of_ages_GCB.pdf` (gitignored); text dump `sierra_extract.txt`.

---

## Why this changed the paper

**Our formula was already theirs.** `doublechecks/intrinsic_mrt.R` computes
`1' (-A)^-1 b` with `b` the normalised input vector. That is their **Eqn (2)**:

> MTT = (1,…,1) · A⁻¹ · I/ΣI

and their **Eqn (4)** is our numerical route:

> "It is important to note that the expression for mean transit time above (Eqn 2) is equivalent to
> MTT = Σx_ss / ΣI … That is, the ratio of the total stocks at steady state to the total inputs is
> equivalent to the mean transit time. This ratio is the turnover time as defined previously, and it
> is **only equivalent to mean transit time in the autonomous case at steady state**."

⭐ **Eqn (4) also fixed a separate defect** — the Methods said "the sum of the stocks is the mean
residence time, in years", which is a mass, not a time. The quantity is a stock **divided by** a
flux; ours reads as a bare stock only because ΣI = 1.

## The two quotes that carry the rename

> "All these different uses of the term **residence time** make it difficult to unambiguously apply
> it in carbon cycle research. We therefore do not consider this term any further in this manuscript
> and **discourage its use in further research** unless clear definitions are presented that differ
> from those adopted for system age and transit time."

> "In **autonomous multipool systems at steady state, the turnover time is equivalent to the mean
> transit time**."

The second is exactly our construction — the model frozen at the reference climate, at steady state,
multipool — so all three names describe the same number for us, and only for us.

## The distinctions, so they are not muddled again

| Concept | What it is |
|---|---|
| **System age** | age of the carbon **currently in** the system |
| **Transit time** | age of the carbon **leaving** the system — "the time it takes for a particle to transit a system". **This is ours.** |
| **Turnover time** | stock ÷ flux. They *adopt* this definition, but **discourage** its other use as the inverse of a first-order rate |
| **Residence time** | used in the literature for all three. **Discouraged.** |

⚠ From the Fig. 1 caption: *"the concepts of system age and transit time do not rely on assumptions
about model structure, steady state, or whether the system is autonomous."* Because transit time is
defined for time-varying systems too, the Methods must say — and now do — that we report the transit
time of the **frozen** system at the reference condition, not of the real soil.

## Secondary

Rasmussen, M., Hastings, A., Smith, M. J. *et al.* (2016), *Transit times and mean ages for
nonautonomous and autonomous compartmental systems*, **J. Math. Biol. 73:1379–1398**,
doi 10.1007/s00285-016-0990-8 — the source of the closed forms Sierra et al. quote. Cited in the
Methods alongside them. ⬜ Not read; cited only for the formula Sierra attributes to it.

---

⚠ **The code still says `MRT`** — deliberate, because the figure staleness guards key on those
filenames. See `RELEASE_NOTES.md`.
