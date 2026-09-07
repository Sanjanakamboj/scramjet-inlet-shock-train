#!/usr/bin/env python3
"""Sweep single-ramp oblique-shock behavior across Mach number and ramp
angle, and compare the representative M1=6, theta=10 deg weak/strong
oblique-shock cases against the M1=6 normal shock from Milestone 1.

Generic perfect-gas reduced-order study -- not experimentally validated.
This is still NOT an inlet model: a single 2-D compression ramp only.
"""

from __future__ import annotations

import numpy as np

from scramjet_inlet.normal_shock import normal_shock
from scramjet_inlet.oblique_shock import detachment_limit, oblique_shock

GAMMA = 1.4
M1_LIST = [2.0, 3.0, 4.0, 5.0, 6.0, 8.0]
SAFE_FRACTIONS = [0.25, 0.5, 0.75, 0.9]

M1_REP = 6.0
T1_REP = 220.0
P1_REP = 2500.0
THETA_REP_DEG = 10.0


def section(title: str) -> None:
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)


def main() -> None:
    section("SINGLE-RAMP OBLIQUE-SHOCK SWEEP -- generic perfect-gas study, NOT validated")

    section(f"theta_max(M1) and weak-branch samples for M1 in {M1_LIST}")
    header = (
        f"{'M1':>5} | {'theta_max(deg)':>14} | {'theta(deg)':>10} | {'beta_weak(deg)':>14} | "
        f"{'p2/p1':>8} | {'p0,2/p0,1':>10} | {'M2':>7}"
    )
    print(header)
    print("-" * len(header))
    for M1 in M1_LIST:
        limit = detachment_limit(M1, GAMMA)
        theta_max_deg = np.degrees(limit.theta_max)
        for frac in SAFE_FRACTIONS:
            theta = frac * limit.theta_max
            r = oblique_shock(M1, theta, GAMMA, branch="weak")
            print(
                f"{M1:5.2f} | {theta_max_deg:14.4f} | {np.degrees(theta):10.4f} | "
                f"{np.degrees(r.beta):14.4f} | {r.p2_p1:8.4f} | {r.p02_p01:10.6f} | {r.M2:7.4f}"
            )
        print("-" * len(header))

    section(f"theta sweep at M1={M1_REP} (small angle up to ~95% of theta_max), weak branch")
    limit_rep = detachment_limit(M1_REP, GAMMA)
    theta_max_rep_deg = np.degrees(limit_rep.theta_max)
    beta_at_theta_max_deg = np.degrees(limit_rep.beta_at_theta_max)
    print(f"theta_max(M1={M1_REP}) = {theta_max_rep_deg:.4f} deg")
    print(f"beta_at_theta_max      = {beta_at_theta_max_deg:.4f} deg")

    theta_fracs = np.linspace(0.02, 0.95, 12)
    header2 = (
        f"{'theta(deg)':>10} | {'beta_weak(deg)':>14} | {'p2/p1':>8} | "
        f"{'p0,2/p0,1':>10} | {'M2':>7}"
    )
    print("\n" + header2)
    print("-" * len(header2))
    for frac in theta_fracs:
        theta = frac * limit_rep.theta_max
        r = oblique_shock(M1_REP, theta, GAMMA, branch="weak")
        print(
            f"{np.degrees(theta):10.4f} | {np.degrees(r.beta):14.4f} | {r.p2_p1:8.4f} | "
            f"{r.p02_p01:10.6f} | {r.M2:7.4f}"
        )

    section(
        f"Representative case: M1={M1_REP}, theta={THETA_REP_DEG} deg "
        f"(T1={T1_REP} K, p1={P1_REP} Pa)"
    )
    theta_rep = np.radians(THETA_REP_DEG)
    r_weak = oblique_shock(M1_REP, theta_rep, GAMMA, branch="weak", p1=P1_REP, T1=T1_REP)
    r_strong = oblique_shock(M1_REP, theta_rep, GAMMA, branch="strong", p1=P1_REP, T1=T1_REP)
    ns = normal_shock(M1_REP, GAMMA)

    def report(label: str, r) -> None:
        print(f"\n--- {label} ---")
        print(f"  beta                  = {np.degrees(r.beta):.4f} deg")
        print(f"  Mn1                   = {r.Mn1:.4f}")
        print(f"  Mn2                   = {r.Mn2:.4f}")
        print(f"  M2                    = {r.M2:.4f}")
        print(f"  p2/p1 (static rise)   = {r.p2_p1:.4f}")
        print(f"  rho2/rho1             = {r.rho2_rho1:.4f}")
        print(f"  T2/T1                 = {r.T2_T1:.4f}")
        if r.p2 is not None:
            print(f"  p2                    = {r.p2:.2f} Pa")
            print(f"  rho2                  = {r.rho2:.5f} kg/m^3")
            print(f"  T2                    = {r.T2:.2f} K")
            print(f"  V2                    = {r.V2:.2f} m/s")
        print(f"  p0,2/p0,1 (recovery)  = {r.p02_p01:.6f}")
        loss_pct = (1.0 - r.p02_p01) * 100.0
        print(f"  Stagnation-pressure loss = {loss_pct:.3f} %")
        if r.T0_1 is not None:
            resid = abs(r.T0_2 - r.T0_1) / r.T0_1
            print(f"  T0 conservation residual = {resid:.3e}")

    print(f"Mach angle mu           = {np.degrees(np.arcsin(1.0/M1_REP)):.4f} deg")
    print(f"theta_max(M1={M1_REP})   = {theta_max_rep_deg:.4f} deg")
    report("WEAK oblique shock (expected external-inlet operating branch)", r_weak)
    report(
        "STRONG oblique shock (alternative mathematical root; NOT typical external-inlet use)",
        r_strong,
    )

    print("\n--- M1 NORMAL SHOCK (Milestone 1 reference) ---")
    print(f"  M2                    = {float(ns.M2):.4f}")
    print(f"  p2/p1 (static rise)   = {float(ns.p2_p1):.4f}")
    print(f"  p0,2/p0,1 (recovery)  = {float(ns.p02_p01):.6f}")
    loss_ns = (1.0 - float(ns.p02_p01)) * 100.0
    print(f"  Stagnation-pressure loss = {loss_ns:.3f} %")

    section("WEAK OBLIQUE vs. STRONG OBLIQUE vs. NORMAL SHOCK -- comparison")
    print(
        f"{'Case':>28} | {'beta(deg)':>10} | {'M2':>7} | {'p2/p1':>8} | "
        f"{'p0,2/p0,1':>10} | {'loss %':>8}"
    )
    print("-" * 88)
    print(
        f"{'Weak oblique (theta=10deg)':>28} | {np.degrees(r_weak.beta):10.4f} | "
        f"{r_weak.M2:7.4f} | {r_weak.p2_p1:8.4f} | {r_weak.p02_p01:10.6f} | "
        f"{(1-r_weak.p02_p01)*100:8.3f}"
    )
    print(
        f"{'Strong oblique (theta=10deg)':>28} | {np.degrees(r_strong.beta):10.4f} | "
        f"{r_strong.M2:7.4f} | {r_strong.p2_p1:8.4f} | {r_strong.p02_p01:10.6f} | "
        f"{(1-r_strong.p02_p01)*100:8.3f}"
    )
    print(
        f"{'Normal shock':>28} | {'90 (n/a: normal)':>10} | "
        f"{float(ns.M2):7.4f} | {float(ns.p2_p1):8.4f} | {float(ns.p02_p01):10.6f} | "
        f"{loss_ns:8.3f}"
    )

    print(
        "\nEngineering takeaway: at the same ramp deflection, the weak oblique shock "
        "retains dramatically more total pressure than either the strong oblique root "
        "or a single terminal normal shock -- this is the physical motivation for "
        "distributed oblique compression in hypersonic inlets. This is NOT yet a "
        "complete inlet result (single ramp only; see README.md / DESIGN.md)."
    )


if __name__ == "__main__":
    main()
