#!/usr/bin/env python3
"""Generate Milestone 2 figures (deterministic, no randomness):

  figures/fig4_theta_beta_M_diagram.png
  figures/fig5_beta_vs_theta_M1_6.png
  figures/fig6_compression_and_recovery_vs_theta_M1_6.png
  figures/fig7_weak_vs_strong_vs_normal_M1_6_theta10.png

All figures are generic perfect-gas reduced-order results for a SINGLE
compression ramp and are not experimentally validated. This is still NOT
an inlet model.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scramjet_inlet.normal_shock import normal_shock
from scramjet_inlet.oblique_shock import (
    detachment_limit,
    mach_angle,
    oblique_shock,
    theta_from_beta,
)

GAMMA = 1.4
FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"
CAPTION = (
    "Generic perfect-gas reduced-order study — not experimentally validated (single ramp only)"
)

M1_REP = 6.0
THETA_REP_DEG = 10.0
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


def _theta_beta_curve(M1: float, n: int = 600):
    """Return (beta_deg, theta_deg) full weak+strong curve for a given M1."""
    mu = mach_angle(M1)
    betas = np.linspace(mu + 1e-6, np.pi / 2.0 - 1e-6, n)
    thetas = np.array([theta_from_beta(M1, b, GAMMA) for b in betas])
    return np.degrees(betas), np.degrees(thetas)


def make_figure4() -> None:
    """theta-beta-M diagram for several Mach numbers."""
    M1_curves = [2.0, 3.0, 4.0, 6.0, 8.0]
    colors = ["#1f4e79", "#b3261e", "#1b7a3d", "#7a4fa3", "#c26b0a"]

    fig, ax = plt.subplots(figsize=(9.0, 7.0))
    for M1, color in zip(M1_curves, colors, strict=True):
        beta_deg, theta_deg = _theta_beta_curve(M1)
        ax.plot(theta_deg, beta_deg, color=color, lw=2.0, label=rf"$M_1={M1:g}$")

        limit = detachment_limit(M1, GAMMA)
        ax.plot(
            np.degrees(limit.theta_max),
            np.degrees(limit.beta_at_theta_max),
            marker="o",
            color=color,
            markersize=5,
            zorder=5,
        )
        mu = mach_angle(M1)
        ax.plot(0.0, np.degrees(mu), marker="|", color=color, markersize=10, zorder=5)

    ax.set_xlabel(r"Deflection angle $\theta$ (deg)")
    ax.set_ylabel(r"Shock angle $\beta$ (deg)")
    ax.set_title(
        r"$\theta$-$\beta$-$M$ diagram ($\gamma=1.4$): weak branch (lower), "
        "strong branch (upper)\n"
        "Dots mark maximum deflection (detachment); tick marks on the "
        r"$\beta$-axis mark the Mach angle $\mu$"
    )
    ax.set_xlim(0, 50)
    ax.set_ylim(0, 90)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", fontsize=9, title=r"$M_1$")

    fig.tight_layout(rect=(0, 0.03, 1, 1))
    _stamp_caption(fig)
    fig.savefig(FIGURES_DIR / "fig4_theta_beta_M_diagram.png", dpi=200)
    plt.close(fig)


def make_figure5() -> None:
    """beta_weak(theta) and beta_strong(theta) for M1=6."""
    limit = detachment_limit(M1_REP, GAMMA)
    theta_max_deg = np.degrees(limit.theta_max)

    beta_deg, theta_deg = _theta_beta_curve(M1_REP, n=2000)

    # Split the full curve into weak (theta increasing) and strong (theta
    # decreasing) branches at the theta_max index for clean plotting.
    idx_max = int(np.argmax(theta_deg))
    weak_theta, weak_beta = theta_deg[: idx_max + 1], beta_deg[: idx_max + 1]
    strong_theta, strong_beta = theta_deg[idx_max:], beta_deg[idx_max:]

    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    ax.plot(
        weak_theta, weak_beta, color="#1f4e79", lw=2.2,
        label=r"weak branch $\beta_{\rm weak}(\theta)$",
    )
    ax.plot(
        strong_theta, strong_beta, color="#b3261e", lw=2.2,
        label=r"strong branch $\beta_{\rm strong}(\theta)$",
    )

    r_weak = oblique_shock(M1_REP, np.radians(THETA_REP_DEG), GAMMA, branch="weak")
    r_strong = oblique_shock(M1_REP, np.radians(THETA_REP_DEG), GAMMA, branch="strong")
    ax.scatter(
        [THETA_REP_DEG, THETA_REP_DEG],
        [np.degrees(r_weak.beta), np.degrees(r_strong.beta)],
        color="black",
        zorder=6,
        s=45,
        label=rf"$\theta={THETA_REP_DEG:g}^\circ$ representative case",
    )
    ax.axvline(THETA_REP_DEG, color="gray", lw=0.8, ls="--", alpha=0.6)
    ax.plot(
        theta_max_deg,
        np.degrees(limit.beta_at_theta_max),
        marker="*",
        color="#c26b0a",
        markersize=14,
        zorder=6,
        label=rf"branch coalescence at $\theta_{{max}}={theta_max_deg:.1f}^\circ$",
    )

    ax.set_xlabel(r"Ramp deflection angle $\theta$ (deg)")
    ax.set_ylabel(r"Shock angle $\beta$ (deg)")
    ax.set_title(rf"Shock angle vs. ramp angle at $M_1={M1_REP:g}$: dual-root structure")
    ax.set_xlim(0, theta_max_deg * 1.05)
    ax.set_ylim(0, 90)
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right", fontsize=8.5)

    fig.tight_layout(rect=(0, 0.03, 1, 1))
    _stamp_caption(fig)
    fig.savefig(FIGURES_DIR / "fig5_beta_vs_theta_M1_6.png", dpi=200)
    plt.close(fig)


def make_figure6() -> None:
    """p2/p1 and p0,2/p0,1 vs theta for M1=6, weak branch, two aligned panels."""
    limit = detachment_limit(M1_REP, GAMMA)
    theta_fracs = np.linspace(0.01, 0.97, 200)
    thetas = theta_fracs * limit.theta_max
    p2_p1_vals = np.array([oblique_shock(M1_REP, t, GAMMA, branch="weak").p2_p1 for t in thetas])
    p02_p01_vals = np.array(
        [oblique_shock(M1_REP, t, GAMMA, branch="weak").p02_p01 for t in thetas]
    )
    theta_deg = np.degrees(thetas)

    theta_rep = np.radians(THETA_REP_DEG)
    r_rep = oblique_shock(M1_REP, theta_rep, GAMMA, branch="weak")
    ns = normal_shock(M1_REP, GAMMA)

    fig, (ax_p, ax_p0) = plt.subplots(2, 1, figsize=(8.5, 8.0), sharex=True)

    ax_p.plot(theta_deg, p2_p1_vals, color="#b3261e", lw=2.2)
    ax_p.scatter(
        [THETA_REP_DEG], [r_rep.p2_p1], color="black", zorder=5, s=50,
        label=rf"$\theta=10^\circ$: $p_2/p_1={r_rep.p2_p1:.2f}$",
    )
    ax_p.axhline(
        float(ns.p2_p1), color="gray", lw=1.0, ls=":",
        label=rf"$M_1$ normal shock: $p_2/p_1={float(ns.p2_p1):.1f}$",
    )
    ax_p.set_ylabel(r"Static pressure ratio $p_2/p_1$")
    ax_p.set_title(
        rf"Compression and stagnation-pressure recovery vs. ramp angle "
        rf"($M_1={M1_REP:g}$, weak branch)"
    )
    ax_p.grid(alpha=0.3)
    ax_p.legend(fontsize=8.5, loc="upper left")

    ax_p0.plot(theta_deg, p02_p01_vals, color="#1f4e79", lw=2.2)
    ax_p0.scatter(
        [THETA_REP_DEG], [r_rep.p02_p01], color="black", zorder=5, s=50,
        label=rf"$\theta=10^\circ$: $p_{{0,2}}/p_{{0,1}}={r_rep.p02_p01:.3f}$",
    )
    ax_p0.axhline(
        float(ns.p02_p01), color="gray", lw=1.0, ls=":",
        label=rf"$M_1$ normal shock: $p_{{0,2}}/p_{{0,1}}={float(ns.p02_p01):.3f}$",
    )
    ax_p0.set_xlabel(r"Ramp deflection angle $\theta$ (deg)")
    ax_p0.set_ylabel(r"Stagnation-pressure recovery $p_{0,2}/p_{0,1}$")
    ax_p0.set_ylim(0, 1.02)
    ax_p0.grid(alpha=0.3)
    ax_p0.legend(fontsize=8.5, loc="lower left")
    ax_p0.text(
        0.5,
        0.93,
        "Increasing ramp deflection increases compression but costs total pressure.",
        transform=ax_p0.transAxes,
        ha="center",
        va="top",
        fontsize=9,
        style="italic",
        color="0.3",
        bbox={
            "boxstyle": "round,pad=0.3",
            "facecolor": "white",
            "edgecolor": "0.75",
            "alpha": 0.85,
        },
    )

    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    _stamp_caption(fig)
    fig.savefig(FIGURES_DIR / "fig6_compression_and_recovery_vs_theta_M1_6.png", dpi=200)
    plt.close(fig)


def make_figure7() -> None:
    """Weak oblique vs. strong oblique vs. normal shock comparison, M1=6, theta=10 deg."""
    theta_rep = np.radians(THETA_REP_DEG)
    r_weak = oblique_shock(M1_REP, theta_rep, GAMMA, branch="weak", p1=P1_REP, T1=T1_REP)
    r_strong = oblique_shock(M1_REP, theta_rep, GAMMA, branch="strong", p1=P1_REP, T1=T1_REP)
    ns = normal_shock(M1_REP, GAMMA)

    cases = ["Weak oblique", "Strong oblique", "Normal shock"]
    colors = ["#1f4e79", "#c26b0a", "#b3261e"]

    beta_deg = [np.degrees(r_weak.beta), np.degrees(r_strong.beta), None]  # None: not meaningful
    M2_vals = [r_weak.M2, r_strong.M2, float(ns.M2)]
    p2_p1_vals = [r_weak.p2_p1, r_strong.p2_p1, float(ns.p2_p1)]
    p02_p01_vals = [r_weak.p02_p01, r_strong.p02_p01, float(ns.p02_p01)]
    loss_vals = [(1.0 - v) * 100.0 for v in p02_p01_vals]

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.0))
    ax_beta, ax_M2, ax_p, ax_loss = axes.flat
    x = np.arange(3)

    # beta panel: normal shock has no meaningful "shock angle" in the same
    # sense (it is trivially 90 deg / perpendicular); label it explicitly
    # rather than inventing a comparable number.
    beta_plot_vals = [beta_deg[0], beta_deg[1], 90.0]
    bars = ax_beta.bar(x, beta_plot_vals, color=colors)
    ax_beta.set_xticks(x)
    ax_beta.set_xticklabels(cases, fontsize=8.5)
    ax_beta.set_ylabel(r"Shock angle $\beta$ (deg)")
    ax_beta.set_title(r"Shock angle $\beta$")
    ax_beta.set_ylim(0, 100)
    for xi, val, case in zip(x, beta_plot_vals, cases, strict=True):
        label = f"{val:.1f}°" if case != "Normal shock" else "90° (by definition)"
        ax_beta.text(xi, val + 2, label, ha="center", fontsize=8)
    ax_beta.grid(axis="y", alpha=0.3)

    ax_M2.bar(x, M2_vals, color=colors)
    ax_M2.set_xticks(x)
    ax_M2.set_xticklabels(cases, fontsize=8.5)
    ax_M2.set_ylabel(r"Downstream Mach $M_2$")
    ax_M2.set_title(r"Downstream Mach number $M_2$")
    for xi, val in zip(x, M2_vals, strict=True):
        ax_M2.text(xi, val + 0.08, f"{val:.2f}", ha="center", fontsize=8)
    ax_M2.grid(axis="y", alpha=0.3)

    ax_p.bar(x, p2_p1_vals, color=colors)
    ax_p.set_xticks(x)
    ax_p.set_xticklabels(cases, fontsize=8.5)
    ax_p.set_ylabel(r"Static pressure ratio $p_2/p_1$")
    ax_p.set_title(r"Static pressure rise $p_2/p_1$")
    for xi, val in zip(x, p2_p1_vals, strict=True):
        ax_p.text(xi, val + 1.0, f"{val:.2f}", ha="center", fontsize=8)
    ax_p.grid(axis="y", alpha=0.3)

    ax_loss.bar(x, loss_vals, color=colors)
    ax_loss.set_xticks(x)
    ax_loss.set_xticklabels(cases, fontsize=8.5)
    ax_loss.set_ylabel("Stagnation-pressure loss (%)")
    ax_loss.set_title(r"Stagnation-pressure loss $1-p_{0,2}/p_{0,1}$")
    ax_loss.set_ylim(0, 105)
    for xi, val in zip(x, loss_vals, strict=True):
        ax_loss.text(xi, val + 2, f"{val:.1f}%", ha="center", fontsize=8)
    ax_loss.grid(axis="y", alpha=0.3)
    del bars

    fig.suptitle(
        rf"Weak oblique vs. strong oblique vs. normal shock: $M_1={M1_REP:g}$, "
        rf"$\theta={THETA_REP_DEG:g}^\circ$",
        fontsize=13,
    )
    fig.text(
        0.5, 0.085,
        "Weak-shock compression is normally preferred for external hypersonic inlet "
        "compression:\n"
        "it preserves far more total pressure, at the cost of substantially less "
        "compression per shock.",
        ha="center", va="bottom", fontsize=8.5, style="italic", color="0.3", linespacing=1.6,
    )
    fig.tight_layout(rect=(0, 0.13, 1, 0.93))
    _stamp_caption(fig)
    fig.savefig(FIGURES_DIR / "fig7_weak_vs_strong_vs_normal_M1_6_theta10.png", dpi=200)
    plt.close(fig)


def main() -> None:
    FIGURES_DIR.mkdir(exist_ok=True)
    make_figure4()
    make_figure5()
    make_figure6()
    make_figure7()
    print(f"Milestone 2 figures written to: {FIGURES_DIR}")
    for name in (
        "fig4_theta_beta_M_diagram.png",
        "fig5_beta_vs_theta_M1_6.png",
        "fig6_compression_and_recovery_vs_theta_M1_6.png",
        "fig7_weak_vs_strong_vs_normal_M1_6_theta10.png",
    ):
        print(f"  {name}")


if __name__ == "__main__":
    main()
