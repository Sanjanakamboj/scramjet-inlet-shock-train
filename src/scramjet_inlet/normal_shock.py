"""Classical calorically-perfect-gas normal-shock relations.

Conventions: subscript 1 denotes the upstream (supersonic) state, subscript
2 the downstream (post-shock) state. All relations require M1 > 1 (a
normal shock cannot exist for M1 <= 1 in steady calorically-perfect-gas
flow). See DESIGN.md for the full derivation/source audit.

Fundamental distinction maintained throughout this module:

    STATIC PRESSURE RISE:        p2/p1 > 1   (flow compresses)
    STAGNATION-PRESSURE RECOVERY: p02/p01 < 1  (irreversible entropy rise)

p2/p1 is never referred to as "pressure recovery" anywhere in this project.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from scramjet_inlet.gas_dynamics import DEFAULT_GAMMA

ArrayLike = npt.ArrayLike


@dataclass
class NormalShockResult:
    """Normal-shock solution ratios for given upstream M1, gamma.

    Attributes
    ----------
    M1 : upstream Mach number (input, echoed back)
    M2 : downstream Mach number (< 1 for M1 > 1)
    p2_p1 : static pressure ratio p2/p1 (> 1)
    rho2_rho1 : density ratio rho2/rho1 (> 1)
    T2_T1 : static temperature ratio T2/T1 (> 1)
    p02_p01 : stagnation-pressure recovery p0,2/p0,1 (in (0, 1))
    """

    M1: np.ndarray
    M2: np.ndarray
    p2_p1: np.ndarray
    rho2_rho1: np.ndarray
    T2_T1: np.ndarray
    p02_p01: np.ndarray


def _validate_M1(M1: ArrayLike) -> np.ndarray:
    M1_arr = np.asarray(M1, dtype=float)
    if np.any(M1_arr <= 1.0):
        raise ValueError(
            "Normal-shock relations require M1 > 1 (supersonic upstream flow); "
            f"got M1={M1}"
        )
    return M1_arr


def normal_shock(M1: ArrayLike, gamma: ArrayLike = DEFAULT_GAMMA) -> NormalShockResult:
    """Solve the calorically-perfect-gas normal-shock jump conditions.

    Requires M1 > 1. Implements the exact closed-form relations:

        M2^2 = [1 + (gamma-1)/2 * M1^2] / [gamma * M1^2 - (gamma-1)/2]

        p2/p1 = 1 + [2*gamma/(gamma+1)] * (M1^2 - 1)

        rho2/rho1 = [(gamma+1) * M1^2] / [(gamma-1) * M1^2 + 2]

        T2/T1 = (p2/p1) / (rho2/rho1)

        p0,2/p0,1 = [ ((gamma+1) M1^2) / ((gamma-1) M1^2 + 2) ]^(gamma/(gamma-1))
                    * [ (gamma+1) / (2*gamma*M1^2 - (gamma-1)) ]^(1/(gamma-1))

    Source: Anderson, Modern Compressible Flow (normal-shock chapter);
    cross-checked against NASA Glenn Research Center's normal-shock
    reference page. See DESIGN.md for the full audit and derivation notes.

    Returns a NormalShockResult. Accepts scalar or array-like M1, gamma
    with standard NumPy broadcasting.
    """
    M1_arr = _validate_M1(M1)
    gamma_arr = np.asarray(gamma, dtype=float)
    if np.any(gamma_arr <= 1.0):
        raise ValueError(f"gamma must be > 1; got {gamma}")

    g = gamma_arr
    M1sq = M1_arr**2

    M2sq = (1.0 + 0.5 * (g - 1.0) * M1sq) / (g * M1sq - 0.5 * (g - 1.0))
    M2 = np.sqrt(M2sq)

    p2_p1 = 1.0 + (2.0 * g / (g + 1.0)) * (M1sq - 1.0)

    rho2_rho1 = ((g + 1.0) * M1sq) / ((g - 1.0) * M1sq + 2.0)

    T2_T1 = p2_p1 / rho2_rho1

    p02_p01 = (
        (((g + 1.0) * M1sq) / ((g - 1.0) * M1sq + 2.0)) ** (g / (g - 1.0))
    ) * (
        ((g + 1.0) / (2.0 * g * M1sq - (g - 1.0))) ** (1.0 / (g - 1.0))
    )

    return NormalShockResult(
        M1=M1_arr,
        M2=M2,
        p2_p1=p2_p1,
        rho2_rho1=rho2_rho1,
        T2_T1=T2_T1,
        p02_p01=p02_p01,
    )
