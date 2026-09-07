# DESIGN.md — Milestone 1: Compressible-Flow and Single Normal-Shock Foundation

This document records the conventions, gas model, equations, source audit,
and verification strategy underlying Milestone 1 of the scramjet inlet
shock-train project. It intentionally overclaims nothing: this is a
generic, reduced-order, calorically-perfect-gas study, not a validated
reconstruction of any real vehicle or engine.

## 1. Conventions

- SI units throughout: Pa, K, kg/m^3, m/s, J/(kg K).
- Subscript `1` = upstream (pre-shock) state; subscript `2` = downstream
  (post-shock) state.
- Subscript/naming `0` (written `p0`, `T0`, `rho0` in code and `p_0`,
  `T_0` in figures) = stagnation (total) property.
- `gamma` = ratio of specific heats (cp/cv); `R` = specific gas constant.
- **STATIC PRESSURE RISE** (`p2/p1`) and **STAGNATION-PRESSURE RECOVERY**
  (`p0,2/p0,1`) are always named explicitly and never conflated. `p2/p1`
  is never called "pressure recovery" anywhere in this project's code,
  tests, scripts, or documentation.

## 2. Gas model

Calorically perfect ideal gas: `p = rho * R * T`, constant `gamma`.
Defaults: `gamma = 1.4`, `R = 287.05 J/(kg K)` (representative of air).

Implemented in `src/scramjet_inlet/gas_dynamics.py`:

- `cp = gamma * R / (gamma - 1)`
- `a = sqrt(gamma * R * T)` (speed of sound; requires `T > 0`)
- `V = M * a` (requires `M >= 0`, `T > 0`)
- `rho = p / (R * T)` (ideal-gas density; requires `p > 0`, `T > 0`)
- Isentropic stagnation relations (valid for any `M >= 0`, since they
  follow from the adiabatic-energy definition of `T0` combined with the
  isentropic relation between `T0/T` and `p0/p`, `rho0/rho`):
  - `T0/T = 1 + (gamma-1)/2 * M^2`
  - `p0/p = [1 + (gamma-1)/2 * M^2]^(gamma/(gamma-1))`
  - `rho0/rho = [1 + (gamma-1)/2 * M^2]^(1/(gamma-1))`

No input is silently clipped. Invalid `gamma <= 1`, `R <= 0`, `T <= 0`,
`p <= 0`, or `M < 0` raises `ValueError` with a specific message.
`M = 0` is explicitly allowed (a valid stagnation state) for
`velocity_from_mach`, `T0_over_T`, `p0_over_p`, `rho0_over_rho`.

## 3. Normal-shock equations

Implemented in `src/scramjet_inlet/normal_shock.py`, valid for `M1 > 1`
only (a normal shock cannot exist for `M1 <= 1` in steady flow of a
calorically perfect gas — attempting one for `M1 <= 1` raises
`ValueError`):

```
M2^2      = [1 + (gamma-1)/2 * M1^2] / [gamma*M1^2 - (gamma-1)/2]

p2/p1     = 1 + [2*gamma/(gamma+1)] * (M1^2 - 1)

rho2/rho1 = [(gamma+1)*M1^2] / [(gamma-1)*M1^2 + 2]

T2/T1     = (p2/p1) / (rho2/rho1)

p0,2/p0,1 = { [(gamma+1)*M1^2] / [(gamma-1)*M1^2 + 2] }^(gamma/(gamma-1))
            * { (gamma+1) / [2*gamma*M1^2 - (gamma-1)] }^(1/(gamma-1))
```

### Derivation basis for `p0,2/p0,1`

The exact stagnation-pressure ratio across a normal shock follows
directly from the definition `p0 = p * (p0/p)_isentropic(M)` applied on
each side of the shock, combined with the shock's own static-pressure
ratio:

```
p0,2/p0,1 = (p2/p1) * (p0/p)_isentropic(M2) / (p0/p)_isentropic(M1)
```

Substituting the closed-form `p2/p1` and the isentropic `p0/p(M)`
relation and simplifying (using the shock's `M2^2` and `rho2/rho1`
relations to eliminate `M2`) yields the compact closed form shown above,
written entirely in terms of `M1` and `gamma`. This algebraic route is
exactly the "reconstruction" strategy independently exercised in
`tests/test_normal_shock.py::test_reconstructed_recovery_matches_direct`
and demonstrated numerically in `scripts/manual_check.py`: given
upstream `p1, M1`, compute `p0,1` via the isentropic relation; compute
downstream `p2` via the shock static-pressure ratio; compute `p0,2` from
`p2` and `M2` via the isentropic relation; and confirm
`p0,2/p0,1` (reconstructed) agrees with the closed-form `p0,2/p0,1`
(direct) to numerical precision (tested at `rel=1e-9`).

Stagnation temperature is conserved across an adiabatic shock
(`T0,2 = T0,1`) because a normal shock does no work and adds no heat —
the energy equation `h + V^2/2 = h0 = cp*T0` is unchanged across the
discontinuity for an adiabatic, no-work process, even though the process
is irreversible (entropy rises, so `p0,2 < p0,1`, while `T0` is
unaffected). This is verified in
`test_stagnation_temperature_conserved`.

## 4. Static vs. stagnation distinction

This is the single most important physical distinction in this project
and is treated as a first-class design constraint, not an afterthought:

- **Static pressure rise `p2/p1 > 1`**: the flow physically compresses
  across the shock. This quantity increases without bound as `M1`
  increases (see Fig. 1).
- **Stagnation-pressure recovery `p0,2/p0,1 < 1`**: the fraction of
  "useful" total pressure retained after the irreversible entropy rise
  through the shock. This quantity monotonically *decreases* toward 0 as
  `M1` increases (see Fig. 2) — the shock becomes an increasingly lossy,
  irreversible compression process even though it is compressing the gas
  more strongly in a static sense.

Confusing these two ratios is a common engineering error; this project's
code, tests, docstrings, script output, and figures consistently label
and distinguish them.

## 5. Source audit

Sources actually inspected for this milestone (fetched and read via web
tools during this session):

1. **NASA Glenn Research Center, "Normal Shock Wave"**
   (<https://www.grc.nasa.gov/www/k-12/airplane/normal.html>) — provides
   the calorically-perfect-gas normal-shock relations for `M2`, `p2/p1`,
   `rho2/rho1`, `T2/T1`, and `p0,2/p0,1` in terms of upstream Mach
   number and `gamma`. (Note: this particular NASA page labels its own
   upstream/downstream stations with subscripts `0`/`1` in a
   non-standard way for a *different* purpose — it is not using `0` for
   "stagnation" on that page. We cross-checked algebraically that,
   translated into this project's `1`=upstream/`2`=downstream,
   `0`=stagnation convention, the relations are identical in structure
   to the standard forms below.)
2. **NASA Glenn Research Center, "Isentropic Flow Relations"**
   (<https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html>) —
   provides the isentropic stagnation-to-static `T/T0`, `p/p0`,
   `rho/rho0` relations for a calorically perfect gas.

These were cross-checked against the well-established closed forms
presented in standard compressible-flow references (Anderson, *Modern
Compressible Flow with Historical Perspective*; the normal-shock and
isentropic-flow chapters), which this project's author has independent
working knowledge of; the exact numeric agreement of the two independent
verification paths in `tests/test_normal_shock.py` (production formula
vs. hand-written reference formula vs. reconstructed-from-isentropic-
relations pathway) provides strong internal confirmation that the
implemented equations are the standard perfect-gas forms. No
discrepancies were found between the NASA source and the standard
textbook forms once station-labeling conventions were reconciled.

No PDFs were downloaded; only public HTML reference pages were fetched
via automated web-content retrieval and are not stored in this
repository (per the project's no-downloaded-references policy).

## 6. Representative operating point

Illustrative, generic hypersonic upstream state (not tied to any real
vehicle, altitude, or flight condition):

```
M1    = 6.0
T1    = 220 K
p1    = 2500 Pa
gamma = 1.4
R     = 287.05 J/(kg K)
```

Headline results (`scripts/manual_check.py`):

| Quantity | Value |
|---|---|
| `M2` | 0.4042 |
| `p2/p1` (static rise) | 41.83 |
| `rho2/rho1` | 5.268 |
| `T2/T1` | 7.941 |
| `p0,2/p0,1` (stagnation recovery) | 0.02965 |
| Stagnation-pressure loss | 97.03 % |
| `T0` conservation residual | 0.0 (exact, to machine precision) |

Physical read: the flow compresses by a factor of ~42 in static
pressure, but retains under 3% of its upstream total pressure — a stark
illustration of why single strong normal shocks are a poor hypersonic
inlet-compression strategy, motivating later milestones' move to
oblique-shock and multi-shock compression.

## 7. Independent-verification strategy

`tests/test_gas_dynamics.py` and `tests/test_normal_shock.py` write
reference formulas directly in the test bodies (not by re-calling the
production functions) wherever a hand/closed-form check is meaningful,
per the project's requirement that verification be genuinely
independent. Specific highlights:

- Every closed-form ratio (`M2`, `p2/p1`, `rho2/rho1`, `T2/T1`,
  `p0,2/p0,1`, `T0/T`, `p0/p`, `rho0/rho`) is checked against a
  hand-written formula evaluated separately in the test.
- The `p0,2/p0,1` stagnation-pressure recovery is additionally checked
  via full reconstruction: `p1, M1 -> p0,1` (isentropic) and
  `p2 (shock), M2 (shock) -> p0,2` (isentropic), then compared against
  the closed-form production value
  (`test_reconstructed_recovery_matches_direct`).
- Physical sanity constraints: `M2 < 1`, `p2 > p1`, `T2 > T1`,
  `rho2 > rho1`, `0 < p0,2/p0,1 < 1`, `T0,2 ≈ T0,1`, weak-shock limit
  `M1 -> 1+` recovers the identity state, monotonic worsening of
  recovery and monotonic increase of static compression with `M1`.
- Strong-shock density-ratio limit `(gamma+1)/(gamma-1) = 6` (for
  `gamma=1.4`) is checked as an asymptotic trend, not equality, at
  `M1=8` and `M1=20`.
- Input validation: negative/zero `T`, `p`, `R`; `gamma <= 1`; `M1 <= 1`
  for the shock solver; negative `M`.
- Scalar/vector (NumPy broadcasting) consistency.
- No NaN/Inf across representative sweeps `M1 in [1.05, 8]` and
  `M in [0, 8]`.

63 tests currently pass (`pytest -W error -q`).

## 8. Numerical / physical sanity audit

Performed over `M1 in [1.05, 8]` via `scripts/normal_shock_study.py` and
`tests/test_normal_shock.py::test_no_nan_inf_over_shock_sweep`:

- `M2` remains strictly subsonic (`< 1`) for all `M1 > 1` tested.
- `p2/p1`, `rho2/rho1`, `T2/T1` all remain `> 1` and increase
  monotonically with `M1`.
- `p0,2/p0,1` remains strictly in `(0, 1)`, approaches 1 as `M1 -> 1+`,
  and falls to ~0.85 % at `M1=8` — severe but finite and physically
  interpretable (no zero-crossing, no sign flip, no blow-up).
- `T0` is conserved to machine precision (residual `~1e-16` relative) at
  every sampled point.
- `rho2/rho1` trends toward, but stays below, the strong-shock limiting
  value `(gamma+1)/(gamma-1) = 6` (reaches ~5.57 at `M1=8`, ~5.83 at
  `M1=20`), consistent with the known asymptotic behavior of a perfect
  gas normal shock.
- No discontinuities, NaNs, or infinities appear anywhere on the tested
  sweep.

## 9. Design decisions

- **No SciPy root-finding used.** All Milestone 1 relations are exact
  closed forms in `M1`/`M`/`gamma`; no iterative solve is needed or
  used. SciPy remains a declared dependency for future milestones
  (e.g., oblique-shock beta-theta-M implicit relations) but is not
  imported anywhere in Milestone 1 code.
- **`NormalShockResult` dataclass** groups the five shock ratios
  (`M1, M2, p2_p1, rho2_rho1, T2_T1, p02_p01`) as a single return value
  for clarity and to keep the static/stagnation distinction explicit in
  the attribute names themselves (`p2_p1` vs. `p02_p01`).
- **Explicit `ValueError`s, never silent clipping**, for all physically
  invalid inputs (`gamma <= 1`, `R <= 0`, `T <= 0`, `p <= 0`, `M < 0`,
  `M1 <= 1` for the shock solver) — required so downstream milestones
  (inlet sweeps, sensitivity studies) fail loudly rather than silently
  producing nonphysical output.
- **Figure 1 uses a 2x2 multi-panel layout** rather than one shared axis,
  because `M2` (O(1), decreasing), `p2/p1` (O(10-100), increasing),
  `rho2/rho1` (O(1-6), increasing, bounded), and `T2/T1` (O(1-15),
  increasing, unbounded) have incompatible scales; forcing them onto one
  axis would be visually misleading.
- **Figure 3's normalized-ratio panel uses a log y-axis**, because the
  property changes at `M1=6` span roughly 0.03x (`p0`) to 42x (`p`) —
  more than three orders of magnitude — and a linear axis would render
  the `p0` bar invisible next to the `p` bar.

## 10. Explicit limitations (Milestone 1)

Milestone 1 does **not** model, and no code, test, script, or figure in
this milestone should be read as implying otherwise:

- Oblique shocks or Prandtl-Meyer expansions
- Inlet ramp or cowl geometry
- Multi-shock inlet compression or shock trains
- Isolator flow or boundary-layer/shock interaction
- Inlet unstart
- Combustion or heat addition (Rayleigh flow)
- Friction (Fanno flow)
- Viscous CFD/RANS
- Calibration against any real engine or flight vehicle

All results are for a generic, calorically-perfect-gas, inviscid,
single-normal-shock idealization and are not experimentally validated.
