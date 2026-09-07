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
    gas_dynamics.py   # perfect-gas thermodynamics + isentropic relations
    normal_shock.py   # classical normal-shock relations
tests/                # independent verification (63 tests)
scripts/
    manual_check.py           # full M1=6 calculation, printed
    normal_shock_study.py     # M1 in [1.05, 8] sweep + sample table
    generate_m1_figures.py    # regenerates figures/*.png deterministically
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
