# Scramjet Inlet Shock-Train — Reduced-Order Engineering Study

## Project purpose

A rigorous, portfolio-quality, reduced-order engineering analysis of
inlet compression and shock-train behavior for a **generic** hypersonic
scramjet inlet. This is a generic reduced-order study, not a
reconstruction of any real engine, and not experimentally validated
unless actual validation data are introduced in a future milestone.

## Final project objective

Progress from verified compressible-flow/shock fundamentals into inlet
compression, shock-train modeling, operating sweeps, sensitivity
analysis, and final engineering synthesis, culminating in:

**Pressure recovery + shock-train position vs. Mach.**

## Milestone 1 scope (this repository state)

Milestone 1 establishes the compressible-flow and thermodynamic
foundation: a calorically-perfect-gas model, isentropic stagnation
relations, and the classical single-normal-shock relations (including
the exact stagnation-pressure recovery). **No oblique shocks, inlet
geometry, multiple shocks, or shock-train modeling yet** — see
[Limitations](#limitations) below and [DESIGN.md](DESIGN.md) for the
full scope boundary and source audit.

## Representative operating point

An illustrative generic hypersonic upstream state (not tied to any real
vehicle, altitude, or flight condition):

```
M1    = 6.0
T1    = 220 K
p1    = 2500 Pa
gamma = 1.4
R     = 287.05 J/(kg K)
```

## Headline M1 = 6 results

| Quantity | Value |
|---|---|
| Downstream Mach `M2` | 0.4042 |
| **Static pressure ratio** `p2/p1` | **41.83** |
| Density ratio `rho2/rho1` | 5.268 |
| Temperature ratio `T2/T1` | 7.941 |
| **Stagnation-pressure recovery** `p0,2/p0,1` | **0.02965** |
| Stagnation-pressure loss | 97.03 % |
| `T0` conservation residual | ~0 (machine precision) |

### Static pressure rise vs. stagnation-pressure recovery

These are **not the same quantity** and are never conflated in this
project:

- **`p2/p1` (static pressure rise, > 1):** the flow physically
  compresses across the shock. At `M1=6`, the static pressure rises by
  a factor of ~42.
- **`p0,2/p0,1` (stagnation-pressure recovery, < 1):** the fraction of
  total pressure retained after the irreversible entropy rise through
  the shock. At `M1=6`, under 3% of the upstream total pressure
  survives.

`p2/p1` is never called "pressure recovery" anywhere in this project.
See [DESIGN.md](DESIGN.md) §4 for the full discussion.

## Equations / model summary

Calorically perfect gas (`p = rho R T`, constant `gamma`):

```
cp = gamma R / (gamma - 1)
a  = sqrt(gamma R T)
V  = M a
```

Isentropic stagnation relations (any `M >= 0`):

```
T0/T   = 1 + (gamma-1)/2 * M^2
p0/p   = [1 + (gamma-1)/2 * M^2]^(gamma/(gamma-1))
rho0/rho = [1 + (gamma-1)/2 * M^2]^(1/(gamma-1))
```

Normal-shock relations (`M1 > 1` only):

```
M2^2      = [1 + (gamma-1)/2 M1^2] / [gamma M1^2 - (gamma-1)/2]
p2/p1     = 1 + [2 gamma/(gamma+1)] (M1^2 - 1)
rho2/rho1 = [(gamma+1) M1^2] / [(gamma-1) M1^2 + 2]
T2/T1     = (p2/p1) / (rho2/rho1)
p0,2/p0,1 = {[(gamma+1)M1^2]/[(gamma-1)M1^2+2]}^(gamma/(gamma-1))
            * {(gamma+1)/[2 gamma M1^2 - (gamma-1)]}^(1/(gamma-1))
```

Full derivation notes and source audit: [DESIGN.md](DESIGN.md).

## Repository layout

```
src/scramjet_inlet/
    gas_dynamics.py    # perfect-gas thermodynamics + isentropic relations
    normal_shock.py    # classical normal-shock relations
    oblique_shock.py   # theta-beta-M relation, single-ramp oblique shock
tests/                 # independent verification (97 tests)
scripts/
    manual_check.py           # full M1=6 normal-shock calculation, printed
    normal_shock_study.py     # M1 in [1.05, 8] sweep + sample table
    generate_m1_figures.py    # regenerates Milestone 1 figures/*.png
    oblique_shock_study.py    # theta_max / weak-strong / M1=6 sweeps
    generate_m2_figures.py    # regenerates Milestone 2 figures/*.png
figures/               # generated PNGs (see below)
DESIGN.md              # conventions, equations, source audit, verification
README.md              # this file
```

## Install / run instructions

```bash
pip install -e ".[dev]"
```

Run the verification suite and linter (must both pass cleanly):

```bash
pytest -W error -q
ruff check .
```

Run the scripts:

```bash
python scripts/manual_check.py
python scripts/normal_shock_study.py
python scripts/generate_m1_figures.py
```

`generate_m1_figures.py` regenerates `figures/*.png` deterministically
(no randomness) — running it twice produces byte-identical files.

## Figures

- **`figures/fig1_shock_state_ratios_vs_M1.png`** — `M2`, `p2/p1`,
  `rho2/rho1`, `T2/T1` vs. `M1`, as a 2x2 multi-panel layout (the
  quantities have incompatible scales, so they are not forced onto one
  shared axis).
- **`figures/fig2_stagnation_pressure_recovery_vs_M1.png`** —
  `p0,2/p0,1` vs. `M1`, with the representative `M1=6` point marked.
  Makes the engineering consequence explicit: strong normal shocks
  cause severe total-pressure loss in hypersonic flow.
- **`figures/fig3_representative_M1_6_state_comparison.png`** — engineering
  summary panel for the `M1=6` operating point: normalized
  upstream/downstream property comparison (log scale, since the changes
  span three orders of magnitude), Mach/velocity comparison, and a text
  summary including the `T0` conservation check.

All figures are generated deterministically and each carries an explicit
"generic perfect-gas reduced-order study — not experimentally validated"
caption.

## Assumptions

- Calorically perfect ideal gas (`gamma`, `R` constant; no real-gas,
  dissociation, or vibrational-excitation effects).
- Inviscid, adiabatic, steady 1-D flow.
- Single normal shock; no boundary layers, geometry, or heat addition.

## Limitations

Milestone 1 does **not** yet model:

- Oblique shocks or Prandtl-Meyer expansions
- Inlet ramp or cowl geometry
- Multi-shock inlet compression or shock trains
- Isolator flow or boundary-layer/shock interaction
- Inlet unstart
- Mass-flow capture/spillage
- Combustion or heat addition (Rayleigh flow)
- Friction (Fanno flow)
- Viscous CFD/RANS
- Calibration against any real engine or flight vehicle

All results are for a generic, calorically-perfect-gas idealization and
are not experimentally validated. See [DESIGN.md](DESIGN.md) for the
full scope boundary, source audit, and verification strategy.

---

# Milestone 2 — Single-Ramp Oblique-Shock Compression

## Objective

Implement and independently verify the classical oblique-shock relations
and the theta-beta-M relation for a **single** attached shock generated
by one 2-D compression ramp: quantify compression and stagnation-
pressure recovery vs. ramp angle and Mach number, and establish the
weak/strong solution structure and detachment limit.

**Milestone 2 is still NOT an inlet model.** No multi-ramp geometry,
cowl, reflected shocks, or shock train is implemented — see
[Limitations](#milestone-2-limitations) below.

## theta-beta-M relation

```
tan(theta) = 2*cot(beta) * (M1^2*sin^2(beta) - 1)
                            / (M1^2*(gamma + cos(2*beta)) + 2)
```

where `M1` = upstream Mach, `theta` = ramp/deflection angle, `beta` =
shock angle (from the upstream flow direction), and the Mach angle
`mu = asin(1/M1)`. For an attached shock: `mu < beta < pi/2`, and there
are generally **two** roots (weak/strong branches) for a given
`(M1, theta)` below the detachment limit `theta_max(M1, gamma)`. Full
derivation and source audit: [DESIGN.md](DESIGN.md) §M2-1 – M2-3.

## Weak / strong branches and detachment

- **Weak branch** (`mu < beta < beta*`): the physically expected
  operating branch for external hypersonic inlet compression surfaces.
- **Strong branch** (`beta* < beta < pi/2`): the alternative mathematical
  root, not the typical external-inlet operating branch (physically
  relevant behind blunt bodies / high back-pressure conditions).
- The two branches **merge** at `theta = theta_max(M1, gamma)` (the
  maximum-deflection / detachment limit); for `theta > theta_max`, no
  attached-shock solution exists (`oblique_shock`/`shock_angle` raise
  `ValueError` rather than fabricating a result).

Both branches are located via robust *bracketed* root-finding
(`scipy.optimize.brentq`) inside intervals determined by a bracketed
scalar minimization for `theta_max` (`scipy.optimize.minimize_scalar`,
`method="bounded"`) — never an unconstrained Newton iteration. See
[DESIGN.md](DESIGN.md) §M2-4 – M2-5.

## Representative M1=6, theta=10° result

Same upstream state as Milestone 1 (`M1=6.0, T1=220 K, p1=2500 Pa,
gamma=1.4, R=287.05 J/(kg K)`), with an illustrative 10° compression-ramp
angle (well inside `theta_max(M1=6) ≈ 42.4°`):

| Quantity | Weak oblique | Strong oblique | M1 normal shock (ref.) |
|---|---|---|---|
| `beta` | 17.59° | 87.61° | 90° (by definition) |
| `M2` | 4.648 | 0.414 | 0.404 |
| `p2/p1` (static rise) | 3.668 | 41.760 | 41.833 |
| **`p0,2/p0,1` (recovery)** | **0.8069** | 0.02976 | 0.02965 |
| Stagnation-pressure loss | **19.31 %** | 97.02 % | 97.04 % |

At the *same* 10° deflection, the **weak oblique shock retains ~81% of
upstream total pressure — roughly 27x more than either the strong
oblique root or an equivalent normal shock (~3%)**. This is the physical
motivation for distributed oblique compression in hypersonic inlets.
This single-ramp comparison is illustrative, not a complete inlet
result. Full numeric detail: [DESIGN.md](DESIGN.md) §M2-10.

## Run commands

```bash
python scripts/oblique_shock_study.py
python scripts/generate_m2_figures.py
```

## Figures

- **`figures/fig4_theta_beta_M_diagram.png`** — theta-beta-M diagram for
  `M1 in {2,3,4,6,8}`: weak/strong branches, Mach-angle limit, and
  maximum-deflection (detachment) points.
- **`figures/fig5_beta_vs_theta_M1_6.png`** — `beta_weak(theta)` and
  `beta_strong(theta)` at `M1=6`, showing the dual-root structure, the
  10° representative case, and branch coalescence at `theta_max`.
- **`figures/fig6_compression_and_recovery_vs_theta_M1_6.png`** — `p2/p1`
  and `p0,2/p0,1` vs. `theta` at `M1=6` (weak branch, two aligned
  panels), with the 10° case and the `M1` normal-shock reference marked.
- **`figures/fig7_weak_vs_strong_vs_normal_M1_6_theta10.png`** — compact
  engineering comparison of weak oblique / strong oblique / normal shock
  at `M1=6, theta=10°`: shock angle, `M2`, `p2/p1`, and stagnation-
  pressure loss.

All figures are generated deterministically and carry the same "generic
… not experimentally validated" caveat as the Milestone 1 figures.

## Verification

`tests/test_oblique_shock.py` adds 34 independent-verification tests
(hand-derived reference formulas, not re-calls of production code),
covering the theta-beta-M residual, weak/strong bracket ordering and
merging at `theta_max`, detachment rejection, tangential-velocity
conservation, direct-vs-reconstructed stagnation-pressure recovery,
weak/strong limiting behavior, monotonicity, and no-NaN/Inf sweeps. 97
tests pass in total (63 Milestone 1 + 34 Milestone 2). Full strategy:
[DESIGN.md](DESIGN.md) §M2-11 – M2-12.

## Milestone 2 Limitations

Does **not** implement: multiple ramps / sequential oblique shocks /
multi-shock compression; cowl-generated or reflected shocks; shock-shock
interaction; shock trains; isolator flow; Fanno or Rayleigh flow;
boundary-layer interaction or shock-induced separation; inlet unstart;
mass-flow capture/spillage; combustor coupling; viscous CFD/RANS; or
calibration against any real engine or vehicle. Geometry is limited to
ONE ideal 2-D compression ramp producing ONE attached oblique shock. See
[DESIGN.md](DESIGN.md) §M2-14 for the full boundary.
