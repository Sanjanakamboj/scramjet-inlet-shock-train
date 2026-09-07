"""scramjet_inlet: reduced-order compressible-flow tools for a generic
hypersonic scramjet inlet study.

Milestone 1 scope: calorically-perfect-gas thermodynamics, isentropic
stagnation relations, and single normal-shock relations (including exact
stagnation-pressure recovery). No oblique shocks, inlet geometry, or
shock-train modeling yet -- see README.md and DESIGN.md.
"""

from scramjet_inlet.gas_dynamics import (
    DEFAULT_GAMMA,
    DEFAULT_R,
    T0_over_T,
    cp_from_gamma_R,
    density_ideal_gas,
    p0_over_p,
    rho0_over_rho,
    speed_of_sound,
    velocity_from_mach,
)
from scramjet_inlet.normal_shock import NormalShockResult, normal_shock
from scramjet_inlet.oblique_shock import (
    DetachmentLimit,
    ObliqueShockResult,
    detachment_limit,
    mach_angle,
    oblique_shock,
    shock_angle,
    theta_from_beta,
)

__all__ = [
    "DEFAULT_GAMMA",
    "DEFAULT_R",
    "cp_from_gamma_R",
    "speed_of_sound",
    "velocity_from_mach",
    "density_ideal_gas",
    "T0_over_T",
    "p0_over_p",
    "rho0_over_rho",
    "NormalShockResult",
    "normal_shock",
    "DetachmentLimit",
    "ObliqueShockResult",
    "mach_angle",
    "theta_from_beta",
    "detachment_limit",
    "shock_angle",
    "oblique_shock",
]
