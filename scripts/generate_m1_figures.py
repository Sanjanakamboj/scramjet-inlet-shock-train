#!/usr/bin/env python3
"""Generate Milestone 1 figures (deterministic, no randomness):

  figures/fig1_shock_state_ratios_vs_M1.png
  figures/fig2_stagnation_pressure_recovery_vs_M1.png
  figures/fig3_representative_M1_6_state_comparison.png

All figures are generic perfect-gas reduced-order results and are not
experimentally validated.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scramjet_inlet.gas_dynamics import T0_over_T, p0_over_p, velocity_from_mach
from scramjet_inlet.normal_shock import normal_shock

GAMMA = 1.4
R = 287.05
FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"

CAPTION = "Generic perfect-gas reduced-order study — not experimentally validated"

M1_REP = 6.0
T1_REP = 220.0
P1_REP = 2500.0


def _stamp_caption(fig) -> None:
    fig.text(
        0.5,
        0.005,
        CAPTION,
        ha="center",
        va="bottom",
        fontsize=8,
        style="italic",
        color="0.35",
    )


def make_figure1() -> None:
    """Normal-shock state ratios vs. M1, multi-panel (incompatible scales)."""
    M1 = np.linspace(1.01, 8.0, 400)
    r = normal_shock(M1, GAMMA)

    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.5), sharex=True)
    ax_M2, ax_p, ax_rho, ax_T = axes.flat

    ax_M2.plot(M1, r.M2, color="#1f4e79", lw=2)
    ax_M2.axhline(1.0, color="gray", lw=0.8, ls=":")
    ax_M2.set_ylabel(r"$M_2$")
    ax_M2.set_title(r"Downstream Mach number $M_2$")

    ax_p.plot(M1, r.p2_p1, color="#b3261e", lw=2)
    ax_p.set_ylabel(r"$p_2 / p_1$")
    ax_p.set_title(r"Static pressure ratio $p_2/p_1$")

    ax_rho.plot(M1, r.rho2_rho1, color="#1b7a3d", lw=2)
    ax_rho.axhline(
        (GAMMA + 1.0) / (GAMMA - 1.0), color="gray", lw=0.8, ls=":", label=r"limit $=6$"
    )
    ax_rho.set_xlabel(r"Upstream Mach number $M_1$")
    ax_rho.set_ylabel(r"$\rho_2 / \rho_1$")
    ax_rho.set_title(r"Density ratio $\rho_2/\rho_1$")
    ax_rho.legend(loc="lower right", fontsize=8)

    ax_T.plot(M1, r.T2_T1, color="#7a4fa3", lw=2)
    ax_T.set_xlabel(r"Upstream Mach number $M_1$")
    ax_T.set_ylabel(r"$T_2 / T_1$")
    ax_T.set_title(r"Static temperature ratio $T_2/T_1$")

    for ax in axes.flat:
        ax.grid(alpha=0.3)
        ax.set_xlim(1.0, 8.0)

    fig.suptitle(
        r"Normal-shock state ratios vs. upstream Mach number $M_1$ ($\gamma=1.4$)",
        fontsize=13,
    )
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    _stamp_caption(fig)
    fig.savefig(FIGURES_DIR / "fig1_shock_state_ratios_vs_M1.png", dpi=200)
    plt.close(fig)


def make_figure2() -> None:
    """Stagnation-pressure recovery vs. M1, with M1=6 marked."""
    M1 = np.linspace(1.01, 8.0, 400)
    r = normal_shock(M1, GAMMA)
    rep = normal_shock(M1_REP, GAMMA)

    fig, ax = plt.subplots(figsize=(8.5, 6.0))
    ax.plot(M1, r.p02_p01, color="#1f4e79", lw=2.2, label=r"$p_{0,2}/p_{0,1}$ (normal shock)")
    ax.scatter(
        [M1_REP],
        [float(rep.p02_p01)],
        color="#b3261e",
        zorder=5,
        s=70,
        label=rf"$M_1={M1_REP:g}$: $p_{{0,2}}/p_{{0,1}}={float(rep.p02_p01):.4f}$",
    )
    ax.axvline(M1_REP, color="#b3261e", lw=0.8, ls="--", alpha=0.7)

    ax.set_xlabel(r"Upstream Mach number $M_1$")
    ax.set_ylabel(r"Stagnation-pressure recovery $p_{0,2}/p_{0,1}$")
    ax.set_title(
        "Normal-shock stagnation-pressure recovery vs. $M_1$\n"
        "Strong normal shocks cause severe total-pressure loss in hypersonic flow"
    )
    ax.set_xlim(1.0, 8.0)
    ax.set_ylim(0.0, 1.02)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", fontsize=9)
    ax.text(
        0.02,
        0.05,
        r"Note: $p_{0,2}/p_{0,1}$ is stagnation-pressure recovery, "
        r"NOT the static ratio $p_2/p_1$ (which instead rises $>1$).",
        transform=ax.transAxes,
        fontsize=8,
        color="0.3",
    )

    fig.tight_layout(rect=(0, 0.03, 1, 1))
    _stamp_caption(fig)
    fig.savefig(FIGURES_DIR / "fig2_stagnation_pressure_recovery_vs_M1.png", dpi=200)
    plt.close(fig)


def make_figure3() -> None:
    """Representative M1=6 state comparison: engineering summary panel."""
    r = normal_shock(M1_REP, GAMMA)
    M2 = float(r.M2)
    T2 = T1_REP * float(r.T2_T1)
    p2 = P1_REP * float(r.p2_p1)
    rho2 = P1_REP / (R * T1_REP) * float(r.rho2_rho1)

    T0_1 = T1_REP * float(T0_over_T(M1_REP, GAMMA))
    T0_2 = T2 * float(T0_over_T(M2, GAMMA))
    p0_1 = P1_REP * float(p0_over_p(M1_REP, GAMMA))
    p0_2 = p2 * float(p0_over_p(M2, GAMMA))

    V1 = float(velocity_from_mach(M1_REP, T1_REP, GAMMA, R))
    V2 = float(velocity_from_mach(M2, T2, GAMMA, R))

    fig = plt.figure(figsize=(11.5, 6.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[2.4, 0.7], hspace=0.42, wspace=0.30)
    ax_bar = fig.add_subplot(gs[0, 0])
    ax_M = fig.add_subplot(gs[0, 1])
    ax_text = fig.add_subplot(gs[1, :])
    ax_text.axis("off")

    labels = [
        r"$p\ (\mathrm{Pa})$",
        r"$T\ (\mathrm{K})$",
        r"$\rho\ (\mathrm{kg/m^3})$",
        r"$p_0\ (\mathrm{Pa})$",
        r"$T_0\ (\mathrm{K})$",
    ]
    upstream = [P1_REP, T1_REP, P1_REP / (R * T1_REP), p0_1, T0_1]
    downstream = [p2, T2, rho2, p0_2, T0_2]
    # Normalize each pair by the upstream value so they share one axis.
    # Use a log scale: the property changes span ~0.03x (p0) to ~42x (p),
    # more than three orders of magnitude, so a linear axis would render
    # the p0 bar invisible next to the p bar.
    norm_up = [1.0 for _ in upstream]
    norm_down = [d / u for d, u in zip(downstream, upstream, strict=True)]

    x = np.arange(len(labels))
    width = 0.35
    ax_bar.bar(x - width / 2, norm_up, width, label="Upstream (station 1)", color="#8fa8c4")
    ax_bar.bar(x + width / 2, norm_down, width, label="Downstream (station 2)", color="#b3261e")
    for xi, val in zip(x, norm_down, strict=True):
        va = "bottom" if val >= 1.0 else "top"
        offset = val * 1.15 if val >= 1.0 else val * 0.85
        ax_bar.text(xi + width / 2, offset, f"{val:.2f}×", ha="center", va=va, fontsize=8)
    ax_bar.axhline(1.0, color="gray", lw=0.8, ls=":")
    ax_bar.set_yscale("log")
    ax_bar.set_ylim(0.01, 100)
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels(labels)
    ax_bar.set_ylabel("Value / upstream value (log scale)")
    ax_bar.set_title("Property change across the shock (normalized)")
    ax_bar.legend(fontsize=8, loc="upper right")
    ax_bar.grid(axis="y", which="both", alpha=0.25)

    # Right panel: Mach number and velocity before/after
    cats = [r"$M$", r"$V\ (\mathrm{m/s},\ /10^3)$"]
    up_vals = [M1_REP, V1 / 1000.0]
    down_vals = [M2, V2 / 1000.0]
    xm = np.arange(len(cats))
    ax_M.bar(xm - width / 2, up_vals, width, label="Upstream", color="#8fa8c4")
    ax_M.bar(xm + width / 2, down_vals, width, label="Downstream", color="#b3261e")
    for xi, val in zip(xm - width / 2, up_vals, strict=True):
        ax_M.text(xi, val + 0.12, f"{val:.2f}", ha="center", fontsize=8)
    for xi, val in zip(xm + width / 2, down_vals, strict=True):
        ax_M.text(xi, val + 0.12, f"{val:.2f}", ha="center", fontsize=8)
    ax_M.set_xticks(xm)
    ax_M.set_xticklabels(cats)
    ax_M.set_ylim(0, max(up_vals) * 1.3)
    ax_M.set_title(rf"Mach number and velocity ($M_1={M1_REP:g}\rightarrow M_2={M2:.3f}$)")
    ax_M.legend(fontsize=8)
    ax_M.grid(axis="y", alpha=0.3)

    loss_pct = (1.0 - float(r.p02_p01)) * 100.0
    T0_resid = abs(T0_2 - T0_1) / T0_1
    summary = (
        rf"$p_2/p_1 = {float(r.p2_p1):.2f}$ (static pressure rise)"
        "        "
        rf"$p_{{0,2}}/p_{{0,1}} = {float(r.p02_p01):.4f}$ (stagnation-pressure recovery)"
        "        "
        rf"Stagnation-pressure loss $= {loss_pct:.1f}\%$"
        "\n"
        rf"$T_{{0,2}}/T_{{0,1}} = {T0_2/T0_1:.6f}$ (conserved; residual {T0_resid:.1e})"
        "        "
        "Static rise and stagnation recovery are distinct quantities -- see README.md"
    )
    ax_text.text(
        0.5,
        0.5,
        summary,
        transform=ax_text.transAxes,
        fontsize=9.5,
        ha="center",
        va="center",
        linespacing=1.8,
        bbox={"boxstyle": "round", "facecolor": "#f2f2f2", "edgecolor": "0.6"},
    )

    fig.suptitle(
        rf"Representative hypersonic normal shock: $M_1={M1_REP:g}$, "
        rf"$T_1={T1_REP:g}\,\mathrm{{K}}$, $p_1={P1_REP:g}\,\mathrm{{Pa}}$",
        fontsize=13,
    )
    fig.subplots_adjust(left=0.09, right=0.97, top=0.84, bottom=0.10)
    _stamp_caption(fig)
    fig.savefig(FIGURES_DIR / "fig3_representative_M1_6_state_comparison.png", dpi=200)
    plt.close(fig)


def main() -> None:
    FIGURES_DIR.mkdir(exist_ok=True)
    make_figure1()
    make_figure2()
    make_figure3()
    print(f"Figures written to: {FIGURES_DIR}")
    for f in sorted(FIGURES_DIR.glob("*.png")):
        print(f"  {f.name}")


if __name__ == "__main__":
    main()
