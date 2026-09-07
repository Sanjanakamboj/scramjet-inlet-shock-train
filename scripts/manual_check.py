#!/usr/bin/env python3
"""Print the complete representative M1 = 6 normal-shock calculation.

Representative generic hypersonic upstream state (illustrative only --
not tied to any specific vehicle, altitude, or flight condition):

    M1 = 6.0
    T1 = 220 K
    p1 = 2500 Pa
    gamma = 1.4
    R = 287.05 J/(kg K)

This is a generic perfect-gas reduced-order study -- not experimentally
validated.
"""

from __future__ import annotations

import numpy as np

from scramjet_inlet.gas_dynamics import (
    T0_over_T,
    cp_from_gamma_R,
    density_ideal_gas,
    p0_over_p,
    speed_of_sound,
    velocity_from_mach,
)
from scramjet_inlet.normal_shock import normal_shock

GAMMA = 1.4
R = 287.05
M1 = 6.0
T1 = 220.0
P1 = 2500.0


def main() -> None:
    print("=" * 72)
    print("MILESTONE 1 -- REPRESENTATIVE NORMAL-SHOCK CALCULATION (M1 = 6)")
    print("Generic perfect-gas reduced-order study -- NOT experimentally validated")
    print("=" * 72)

    cp = cp_from_gamma_R(GAMMA, R)
    a1 = speed_of_sound(T1, GAMMA, R)
    V1 = velocity_from_mach(M1, T1, GAMMA, R)
    rho1 = density_ideal_gas(P1, T1, R)
    T0_1 = T1 * T0_over_T(M1, GAMMA)
    p0_1 = P1 * p0_over_p(M1, GAMMA)

    print("\n--- Upstream state (station 1) ---")
    print(f"  gamma                 = {GAMMA:.4f}")
    print(f"  R                     = {R:.4f} J/(kg K)")
    print(f"  cp                    = {cp:.4f} J/(kg K)")
    print(f"  M1                    = {M1:.4f}")
    print(f"  T1                    = {T1:.4f} K")
    print(f"  p1                    = {P1:.4f} Pa")
    print(f"  a1 (speed of sound)   = {a1:.4f} m/s")
    print(f"  V1 (velocity)         = {V1:.4f} m/s")
    print(f"  rho1 (density)        = {rho1:.6f} kg/m^3")
    print(f"  T0,1 (stagnation T)   = {T0_1:.4f} K")
    print(f"  p0,1 (stagnation p)   = {p0_1:.4f} Pa")

    shock = normal_shock(M1, GAMMA)
    M2 = float(shock.M2)
    p2_p1 = float(shock.p2_p1)
    rho2_rho1 = float(shock.rho2_rho1)
    T2_T1 = float(shock.T2_T1)
    p02_p01 = float(shock.p02_p01)

    T2 = T1 * T2_T1
    p2 = P1 * p2_p1
    rho2 = rho1 * rho2_rho1
    a2 = speed_of_sound(T2, GAMMA, R)
    V2 = velocity_from_mach(M2, T2, GAMMA, R)
    T0_2 = T2 * T0_over_T(M2, GAMMA)
    p0_2 = p2 * p0_over_p(M2, GAMMA)

    print("\n--- Normal-shock jump (station 1 -> station 2) ---")
    print(f"  M2                    = {M2:.6f}")
    print(f"  p2/p1 (static rise)   = {p2_p1:.6f}")
    print(f"  rho2/rho1             = {rho2_rho1:.6f}")
    print(f"  T2/T1                 = {T2_T1:.6f}")

    print("\n--- Downstream state (station 2) ---")
    print(f"  p2                    = {p2:.4f} Pa")
    print(f"  rho2                  = {rho2:.6f} kg/m^3")
    print(f"  T2                    = {T2:.4f} K")
    print(f"  a2 (speed of sound)   = {a2:.4f} m/s")
    print(f"  V2 (velocity)         = {V2:.4f} m/s")
    print(f"  T0,2 (stagnation T)   = {T0_2:.4f} K")
    print(f"  p0,2 (stagnation p)   = {p0_2:.4f} Pa")

    T0_resid = abs(T0_2 - T0_1) / T0_1
    loss_pct = (1.0 - p02_p01) * 100.0

    print("\n--- Stagnation-pressure recovery (NOT the same as p2/p1!) ---")
    print(f"  p0,2/p0,1 (direct)          = {p02_p01:.6f}")
    print(f"  p0,2/p0,1 (from p0_2/p0_1)  = {p0_2 / p0_1:.6f}")
    print(f"  Stagnation-pressure loss    = {loss_pct:.4f} %")
    print(f"  T0 conservation residual    = {T0_resid:.3e}  (should be ~0)")

    print("\n--- Summary ---")
    print(f"  STATIC pressure rise:      p2/p1     = {p2_p1:.3f}  (> 1, flow compresses)")
    print(f"  STAGNATION pressure recovery: p0,2/p0,1 = {p02_p01:.3f}  (< 1, entropy rise)")
    print("  These are NOT the same quantity -- see README.md / DESIGN.md.")

    assert np.isfinite([M2, p2_p1, rho2_rho1, T2_T1, p02_p01]).all()


if __name__ == "__main__":
    main()
