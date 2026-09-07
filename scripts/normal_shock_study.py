#!/usr/bin/env python3
"""Sweep upstream Mach number M1 in [1.05, 8.0] through the normal-shock
relations and report representative values.

Reports both the static compression (p2/p1) and the stagnation-pressure
loss (1 - p0,2/p0,1) at each sample point -- these are distinct
quantities; see README.md / DESIGN.md.

Generic perfect-gas reduced-order study -- not experimentally validated.
"""

from __future__ import annotations

import numpy as np

from scramjet_inlet.normal_shock import normal_shock

GAMMA = 1.4
M1_MIN, M1_MAX = 1.05, 8.0
SAMPLE_POINTS = [1.5, 2.0, 3.0, 4.0, 6.0, 8.0]


def main() -> None:
    print("=" * 88)
    print(f"NORMAL-SHOCK SWEEP: M1 in [{M1_MIN}, {M1_MAX}], gamma = {GAMMA}")
    print("Generic perfect-gas reduced-order study -- NOT experimentally validated")
    print("=" * 88)

    M1_sweep = np.linspace(M1_MIN, M1_MAX, 400)
    sweep_result = normal_shock(M1_sweep, GAMMA)

    assert np.all(np.isfinite(sweep_result.M2))
    assert np.all(np.isfinite(sweep_result.p02_p01))
    assert np.all(sweep_result.M2 < 1.0)
    assert np.all(sweep_result.p02_p01 > 0.0) and np.all(sweep_result.p02_p01 < 1.0)

    print(f"\nSweep of {len(M1_sweep)} points completed; all M2 < 1, all p0,2/p0,1 in (0,1).")

    header = (
        f"{'M1':>6} | {'M2':>8} | {'p2/p1':>10} | {'rho2/rho1':>10} | "
        f"{'T2/T1':>8} | {'p0,2/p0,1':>10} | {'loss %':>8}"
    )
    print("\nRepresentative sample points:")
    print(header)
    print("-" * len(header))
    for M1 in SAMPLE_POINTS:
        r = normal_shock(M1, GAMMA)
        loss_pct = (1.0 - float(r.p02_p01)) * 100.0
        print(
            f"{M1:6.2f} | {float(r.M2):8.4f} | {float(r.p2_p1):10.4f} | "
            f"{float(r.rho2_rho1):10.4f} | {float(r.T2_T1):8.4f} | "
            f"{float(r.p02_p01):10.6f} | {loss_pct:8.3f}"
        )

    print("\nSTATIC COMPRESSION (p2/p1) increases monotonically with M1:")
    p2_p1_at_samples = [float(normal_shock(m, GAMMA).p2_p1) for m in SAMPLE_POINTS]
    for m, v in zip(SAMPLE_POINTS, p2_p1_at_samples, strict=True):
        print(f"  M1={m:4.2f}: p2/p1 = {v:10.4f}")

    print("\nSTAGNATION-PRESSURE LOSS (1 - p0,2/p0,1) worsens monotonically with M1:")
    loss_at_samples = [(1.0 - float(normal_shock(m, GAMMA).p02_p01)) * 100.0 for m in SAMPLE_POINTS]
    for m, v in zip(SAMPLE_POINTS, loss_at_samples, strict=True):
        print(f"  M1={m:4.2f}: loss = {v:7.3f} %")

    print(
        "\nEngineering takeaway: static pressure keeps rising with shock strength, "
        "but total-pressure recovery collapses -- a single strong normal shock is a "
        "very poor hypersonic compression strategy. This motivates oblique-shock / "
        "multi-shock inlet compression, which is out of scope for Milestone 1."
    )


if __name__ == "__main__":
    main()
