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

---

# DESIGN.md — Milestone 2: Single-Ramp Oblique-Shock Compression

Milestone 2 extends the Milestone 1 foundation with the classical
oblique-shock relations and the theta-beta-M relation for a SINGLE
attached shock generated by one 2-D compression ramp. **Milestone 2 is
still NOT an inlet model** — no multi-ramp geometry, cowl, reflected
shocks, or shock train is implemented (see §17 for the full boundary).

## M2-1. Sources actually inspected

1. **NASA Glenn Research Center, "Oblique Shock Wave"**
   (<https://www.grc.nasa.gov/www/k-12/airplane/oblique.html>) — states
   the theta-beta-M relation in the algebraically equivalent form
   `cot(theta) = tan(beta) * [ (gamma+1)*M1^2 / (2*M1^2*sin^2(beta) - 1) - 1 ]`,
   describes the tangential-velocity-conserved / normal-velocity-shocked
   decomposition, and gives the downstream-Mach relation via the normal
   component. Cross-checked algebraically: taking the reciprocal and
   applying standard trig identities converts this to the classical
   `tan(theta) = 2*cot(beta) * (M1^2*sin^2(beta)-1) / (M1^2*(gamma+cos(2*beta))+2)`
   form implemented in this project; the two are the same relation.
2. **NASA Glenn Research Center, "Mach Angle"**
   (<https://www.grc.nasa.gov/www/k-12/airplane/machang.html>) — confirms
   `mu = asin(1/M)`.

Cross-checked against the closed-form theta-beta-M relation and normal-
shock-based oblique-shock construction presented in standard
compressible-flow references (Anderson, *Modern Compressible Flow*,
oblique-shock chapter), which this project's author has independent
working knowledge of. The representative M1=6, theta=10 deg numbers this
milestone reports (beta_weak ~= 17.6 deg, M2 ~= 4.65, p2/p1 ~= 3.67)
match the well-known values for this classical textbook case. No
discrepancies were found. No PDFs were downloaded; only public HTML
reference pages were fetched and are not stored in this repository.

## M2-2. Angle and sign conventions

- `M_1` = upstream Mach number (`M_1 > 1`).
- `theta` = flow-deflection / compression-ramp angle, `0 < theta < theta_max(M_1, gamma)`.
- `beta` = shock angle, measured from the upstream flow direction, `mu < beta < pi/2`.
- `mu` = Mach angle = `asin(1/M_1)`.
- All angles are handled in **radians internally** throughout
  `oblique_shock.py`. Degrees are used ONLY at the convenience/reporting
  layer (scripts and figures), always via explicit `np.radians` /
  `np.degrees` conversions, never silently mixed with radians.
- `beta = pi/2` recovers the normal-shock limit (`theta -> 0` on the
  strong branch); `beta = mu` recovers the zero-strength Mach-wave limit
  (`theta -> 0` on the weak branch).

## M2-3. theta-beta-M relation

Implemented in `theta_from_beta(M1, beta, gamma)`:

```
tan(theta) = 2*cot(beta) * (M1^2*sin^2(beta) - 1)
                            / (M1^2*(gamma + cos(2*beta)) + 2)
```

For `beta` in the open admissible domain `(mu, pi/2)`, the numerator
`M1^2*sin^2(beta) - 1` is strictly positive (since `sin(beta) > sin(mu) =
1/M1` there) and the denominator `M1^2*(gamma+cos(2*beta)) + 2` is
strictly positive (since `cos(2*beta) >= -1` and `gamma > 1`, so
`gamma+cos(2*beta) > 0` for `gamma=1.4`), so `tan(theta) > 0` and
`arctan` is evaluated without singularity anywhere in the domain.
`theta(beta)` is 0 at both endpoints (`beta=mu` and `beta=pi/2`) and has
a single interior maximum `theta_max` at some `beta* in (mu, pi/2)` — the
classical theta-beta-M "bubble" shape (see Fig. 4).

## M2-4. Numerical branch-solving method (`shock_angle`)

No unconstrained Newton iteration is used anywhere in this module.
Branch selection is via explicit interval bracketing:

1. `detachment_limit(M1, gamma)` locates `beta*` (the shock angle at
   `theta_max`) using `scipy.optimize.minimize_scalar(..., method=
   "bounded")` over the bracket `(mu+eps, pi/2-eps)` — a robust bracketed
   scalar minimizer appropriate for the smooth, unimodal `theta(beta)`
   objective on this domain (`eps = 1e-9`, `xatol = 1e-14`).
2. Given `beta*`, the weak-branch root for a target `theta < theta_max`
   is found via `scipy.optimize.brentq` (a robust bracketed root finder)
   on `(mu+eps, beta*)`, where `theta(beta) - theta` is guaranteed to
   change sign (negative at `mu`, positive at `beta*`, since `theta(beta)`
   is monotonically increasing there).
3. The strong-branch root is found via `brentq` on `(beta*, pi/2-eps)`,
   where `theta(beta)-theta` again changes sign (positive at `beta*`,
   negative at `pi/2`, since `theta(beta)` is monotonically decreasing
   there).
4. `theta > theta_max` (checked against `detachment_limit`, with a
   `1e-10` numerical tolerance) raises `ValueError` — an explicitly
   DETACHED condition; no beta is fabricated. `theta <= 0` raises
   `ValueError` (this routine is compression-shock-only). `M1 <= 1`
   raises `ValueError`.
5. `theta` within `1e-9` rad of `theta_max` returns `beta*` directly
   (branch coalescence), since both brackets degenerate to zero width at
   that point.
6. If `f(lo)*f(hi) > 0` on either bracket (which should not occur given
   the unimodal structure and the `theta<=theta_max` guard), the solver
   raises `RuntimeError` rather than returning a possibly-wrong root —
   solver failures are reported honestly, never silently swallowed.

This satisfies the project's requirement that the solver "never silently
return the wrong branch": each branch's root is bracketed by construction
(the sign of `theta(beta)-theta` at each bracket endpoint is known
analytically from the monotonicity of `theta(beta)` on that sub-interval),
so `brentq` cannot cross into the other branch's interval.

## M2-5. theta_max(M1, gamma) method

See M2-4 item 1. `detachment_limit` returns both `theta_max` and
`beta_at_theta_max` from the same bounded-minimization call (no
duplicated work). Verified independently in
`tests/test_oblique_shock.py` via: (a) direct hand-formula
theta-beta-M evaluation at the returned `beta_at_theta_max` reproduces
`theta_max`; (b) a fine grid search over `beta` cross-checks the
representative-case `beta_weak` result at M1=6, theta=10 deg; (c)
smoothness/monotonic-increase of `theta_max(M1)` is checked over
`M1 in [1.5, 8]`.

## M2-6. Oblique-shock state derivation through Mn1

Given a solved `beta` (either branch):

```
Mn1 = M1 * sin(beta)
```

The **normal-shock relations from Milestone 1 are reused directly**
(`normal_shock(Mn1, gamma)`, no duplicated formulas in production code)
to obtain `Mn2`, `p2/p1`, `rho2/rho1`, `T2/T1`, `p0,2/p0,1` — because the
oblique shock, viewed in the frame that decomposes the upstream velocity
into components normal and tangential to the shock, behaves EXACTLY like
a normal shock in the normal direction; the tangential velocity
component is conserved across the (ideal, inviscid) shock and carries no
entropy change. This is why the oblique shock's static and stagnation
RATIOS are identical in form to the normal-shock ratios evaluated at
`Mn1` — this project's oblique-shock static/stagnation ratios are
therefore trivially consistent with (indeed, literally computed from)
the already-verified Milestone 1 normal-shock relations.

## M2-7. M2 reconstruction

```
M2 = Mn2 / sin(beta - theta)
```

Verified independently in
`test_tangential_velocity_conservation`: velocities are reconstructed
from scratch (via `speed_of_sound` and `M*a`, not via any
oblique-shock-internal helper) at each station, and
`V1*cos(beta) == V2*cos(beta-theta)` is checked directly — i.e. the
tangential velocity component is confirmed conserved across the shock,
which is the physical justification for the `M2` reconstruction formula
above.

## M2-8. Total-pressure-recovery derivation and verification

`p0,2/p0,1` for the oblique shock is computed as the Milestone-1
normal-shock stagnation-pressure-recovery formula evaluated at `Mn1`
(§M2-6). This is independently reconstructed in
`test_reconstructed_recovery_matches_direct` via the same two-step
pathway used in Milestone 1: `p1, M1 -> p0,1` (isentropic relation) and
`p2 (oblique shock), M2 (full downstream Mach) -> p0,2` (isentropic
relation), and the ratio is checked against the direct value (agreement
to `rel=1e-6`, looser than Milestone 1's `1e-9` only because this path
additionally compounds through the `M2` reconstruction of §M2-7).

As with Milestone 1, the same fundamental distinction is maintained
throughout:

```
STATIC compression:      p2/p1 > 1
TOTAL-pressure recovery: p0,2/p0,1 < 1
```

`p2/p1` is never called "pressure recovery" anywhere in this project.

## M2-9. Weak/strong branch interpretation

For a given `M1` and attached `theta < theta_max`, TWO mathematically
valid shock angles solve the theta-beta-M relation:

- **Weak branch** (`mu < beta_weak < beta*`): the physically realized
  solution for external supersonic/hypersonic compression surfaces (a
  ramp on the outside of a vehicle, an inlet ramp with a downstream
  region open to atmospheric-like back pressure). Downstream flow
  usually remains supersonic (`M2 > 1`) except very close to
  `theta_max`. Retains substantially more total pressure than the strong
  branch or an equivalent normal shock.
- **Strong branch** (`beta* < beta_strong < pi/2`): the alternative
  mathematical root. Downstream flow is always subsonic. Physically
  realized only when a downstream boundary condition (e.g. a blunt body
  or sufficiently high back pressure) forces it; NOT the typical
  operating branch for an external hypersonic inlet ramp. This project
  reports the strong branch purely as the alternative mathematical
  solution for comparison (§M2-10), not as the expected inlet-design
  branch.

Both branches satisfy the identical theta-beta-M relation for the same
`(M1, theta)`; they differ only in which root of the (generally
two-valued) equation is selected. They merge into a single root at
`theta = theta_max`.

## M2-10. Representative numerical case: M1=6, theta=10 deg

Upstream state: `M1=6.0, T1=220 K, p1=2500 Pa, gamma=1.4, R=287.05
J/(kg K)` (same as Milestone 1). `theta = 10 deg` chosen as a
representative, illustrative external-compression-ramp angle (well
inside `theta_max(M1=6) ~= 42.4 deg`, i.e. a physically reasonable
single-ramp deflection).

| Quantity | Weak branch | Strong branch | M1 normal shock (ref.) |
|---|---|---|---|
| `beta` | 17.59 deg | 87.61 deg | 90 deg (by definition) |
| `Mn1` | 1.813 | 5.995 | 6.000 |
| `Mn2` | 0.614 | 0.404 | 0.404 |
| `M2` | 4.648 | 0.414 | 0.404 |
| `p2/p1` | 3.668 | 41.760 | 41.833 |
| `rho2/rho1` | 2.380 | 5.267 | 5.268 |
| `T2/T1` | 1.541 | 7.928 | 7.941 |
| `p0,2/p0,1` | **0.8069** | 0.02976 | 0.02965 |
| Stagnation-pressure loss | **19.31 %** | 97.02 % | 97.04 % |
| `T0` conservation residual | ~1e-15 | ~1e-15 | 0 (exact) |

Engineering read: at the *same* 10-degree deflection, the weak oblique
shock retains **~81% of upstream total pressure** versus **~3%** for
either the strong-branch root or an equivalent normal shock — a factor
of ~27x more retained total pressure. This is the physical motivation
for distributed oblique compression in hypersonic inlet design.
**This single-ramp result is illustrative, not a complete inlet
result** — it does not yet represent an inlet's overall pressure
recovery (that requires multiple ramps / a full compression system,
which is out of scope for Milestone 2).

## M2-11. Independent-verification strategy

`tests/test_oblique_shock.py` contains 34 tests, all writing reference
formulas/reconstructions directly in the test bodies rather than
re-calling the production functions under test. Highlights beyond those
already covered in M2-3 through M2-8:

- Weak/strong bracket ordering (`mu < beta_weak < beta* < beta_strong <
  pi/2`) and both roots satisfying the *same* input `theta`.
- Branch coalescence near `theta_max` (both branches converge to `beta*`
  within `1e-3` rad as `theta -> theta_max`).
- Rejection of `theta > theta_max` (detached), `theta <= 0`, `M1 <= 1`,
  invalid `gamma`, invalid `branch` string.
- `Mn1 > 1` (supersonic normal component) on both branches.
- Normal-shock ratios at `Mn1` matching an independent hand formula
  (not `normal_shock()` called a second time — a separately written
  closed-form reference in the test file).
- Weak-branch `theta -> 0+` limit: `beta -> mu`, `p2/p1 -> 1`,
  `p0,2/p0,1 -> 1`, `M2 -> M1`.
- Strong-branch `theta -> 0+` limit: `beta -> pi/2` and the ratios
  approach the Milestone-1 normal-shock solution at the same `M1`.
- Monotonicity on the weak branch: `p2/p1` increases and `p0,2/p0,1`
  decreases with `theta`.
- Continuity/finiteness of `beta_weak` over `M1 in [1.5, 8]` at fixed
  `theta` (skipping combinations where `theta` would be detached at that
  `M1`).
- Representative-case cross-check against an independent fine grid
  search over `beta` (200,000 points), not the production Brent solve.
- No NaN/Inf over an admissible `(M1, theta)` grid up to 80% of
  `theta_max`.
- Smoothness of `theta_max(M1)` across `M1 in [1.5, 8]`.

97 tests pass in total (`pytest -W error -q`): 63 from Milestone 1
(unchanged) + 34 new for Milestone 2.

## M2-12. Physical sanity audit

Explicitly checked (see `scripts/oblique_shock_study.py` output and the
project's automated sanity-check pass) across `M1 in {2,3,4,5,6,8}`:

- Weak `beta` lies only slightly above the Mach angle `mu` at small
  `theta` (e.g. at `M1=6`, small `theta`: `beta_weak ~= 11.0 deg` vs.
  `mu = 9.59 deg`).
- Strong `beta` lies much closer to 90 deg at small `theta` (e.g. at
  `M1=6`: `beta_strong ~= 89.5 deg`).
- Weak and strong roots merge at `theta_max` for every tested `M1`.
- `p2/p1` increases and `p0,2/p0,1` decreases monotonically with `theta`
  on the weak branch for every tested `M1`.
- The strong branch has substantially larger total-pressure loss than
  the weak branch at the same `(M1, theta)`, for every tested case.
- The normal shock produces the largest loss of the three compared
  representative M1=6, theta=10 deg cases (though only marginally more
  than the strong oblique root, as expected since the strong root
  approaches the normal-shock limit).
- `M2` remains physically meaningful (supersonic on the weak branch
  except very near `theta_max`; subsonic on the strong branch and for
  the normal shock) across all tested cases.
- `T0` is conserved to within ~1e-15 relative at every sampled point.
- No `beta` below the Mach angle and no attached solution above
  `theta_max` was ever returned (`shock_angle` raises `ValueError`
  instead).
- No discontinuities or branch-swapping observed in any sweep, including
  at high Mach (`M1=8`).

## M2-13. Design decisions

- **SciPy is now used** (`scipy.optimize.brentq`,
  `scipy.optimize.minimize_scalar`) — the first genuine use of SciPy in
  this project, justified by the theta-beta-M relation's dual-root
  branch structure, which an exact closed-form inverse does not cleanly
  provide for arbitrary `(M1, theta)`. Both solvers are bracketed
  methods, per the project's numerical-robustness requirement; no
  unconstrained Newton iteration is used anywhere.
- **`shock_angle`/`oblique_shock` are scalar-only by design.** Branch
  bracketing requires locating `theta_max` per-`(M1, gamma)` pair, which
  is not cleanly vectorizable through SciPy's scalar bracketed solvers
  without looping in Python anyway; callers needing a sweep (scripts,
  tests) loop explicitly over scalar calls. The underlying
  `theta_from_beta` building block IS vectorizable (plain NumPy
  arithmetic) and is exercised at array `beta` in the figure-generation
  scripts.
- **`ObliqueShockResult` carries optional dimensional fields** (`p1, T1,
  p2, rho2, T2, V1, V2, T0_1, T0_2`), populated only when `p1`/`T1` are
  supplied to `oblique_shock`, keeping the non-dimensional core API
  usable without requiring a full thermodynamic state.
- **Figure 4 (theta-beta-M diagram) caps the theta-axis at 50 deg**
  and the beta-axis at 90 deg for readability across `M1 in
  {2,3,4,6,8}`; `M1=5` was intentionally omitted from the plotted curves
  (though included in the numerical sweep) to keep the 5-curve legend
  legible per the project's "not an unreadable spaghetti plot"
  requirement.

## M2-14. Explicit limitations (Milestone 2)

Milestone 2 is **still NOT an inlet model**. It does not implement, and
no code/test/script/figure in this milestone should be read as implying:

- Multiple ramps, sequential oblique shocks, or multi-shock compression
- Cowl-generated shocks, reflected shocks, or shock-shock interaction
- Shock trains or isolator flow
- Fanno flow (friction) or Rayleigh flow (heat addition)
- Boundary-layer interaction or shock-induced separation
- Inlet unstart
- Mass-flow capture/spillage
- Combustor coupling or heat addition
- Viscous CFD/RANS
- Calibration against any real engine or flight vehicle

Geometry is limited to ONE ideal 2-D compression ramp producing ONE
attached oblique shock. All results remain generic, calorically-perfect-
gas, inviscid idealizations and are not experimentally validated.
