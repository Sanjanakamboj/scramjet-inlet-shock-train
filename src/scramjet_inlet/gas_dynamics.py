"""Calorically-perfect-gas thermodynamics and isentropic flow relations.

All functions use SI units (Pa, K, kg/m^3, m/s, J/(kg K)) and a calorically
perfect ideal gas: p = rho * R * T, with constant gamma = cp/cv.

Inputs may be Python scalars or NumPy arrays; return types follow NumPy's
usual broadcasting behavior. Invalid inputs raise ValueError with a
physically meaningful message -- values are never silently clipped.

Default gas properties (representative of air at moderate temperature):
    DEFAULT_GAMMA = 1.4
    DEFAULT_R     = 287.05 J/(kg K)

Domain notes (see DESIGN.md for derivations/sources):
    - cp_from_gamma_R, speed_of_sound, density_ideal_gas: require T > 0
      (and p > 0 for density). No Mach-number restriction.
    - velocity_from_mach: accepts M >= 0 (M = 0 is a valid, physical
      stagnation state).
    - T0_over_T, p0_over_p, rho0_over_rho: accept M >= 0. These are general
      isentropic relations and are NOT restricted to supersonic flow.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

DEFAULT_GAMMA: float = 1.4
DEFAULT_R: float = 287.05  # J/(kg K), dry air


ArrayLike = npt.ArrayLike


def _validate_gamma(gamma: ArrayLike) -> np.ndarray:
    gamma_arr = np.asarray(gamma, dtype=float)
    if np.any(gamma_arr <= 1.0):
        raise ValueError(
            f"gamma must be > 1 for a calorically perfect gas; got {gamma}"
        )
    return gamma_arr


def _validate_R(R: ArrayLike) -> np.ndarray:
    R_arr = np.asarray(R, dtype=float)
    if np.any(R_arr <= 0.0):
        raise ValueError(f"Specific gas constant R must be > 0 J/(kg K); got {R}")
    return R_arr


def _validate_positive_T(T: ArrayLike) -> np.ndarray:
    T_arr = np.asarray(T, dtype=float)
    if np.any(T_arr <= 0.0):
        raise ValueError(f"Temperature T must be > 0 K; got {T}")
    return T_arr


def _validate_positive_p(p: ArrayLike) -> np.ndarray:
    p_arr = np.asarray(p, dtype=float)
    if np.any(p_arr <= 0.0):
        raise ValueError(f"Pressure p must be > 0 Pa; got {p}")
    return p_arr


def _validate_nonnegative_mach(M: ArrayLike) -> np.ndarray:
    M_arr = np.asarray(M, dtype=float)
    if np.any(M_arr < 0.0):
        raise ValueError(f"Mach number M must be >= 0; got {M}")
    return M_arr


def cp_from_gamma_R(gamma: ArrayLike = DEFAULT_GAMMA, R: ArrayLike = DEFAULT_R):
    """Specific heat at constant pressure: cp = gamma * R / (gamma - 1).

    Valid for any gamma > 1, R > 0. No temperature or Mach restriction.
    Returns J/(kg K).
    """
    gamma_arr = _validate_gamma(gamma)
    R_arr = _validate_R(R)
    return gamma_arr * R_arr / (gamma_arr - 1.0)


def speed_of_sound(T: ArrayLike, gamma: ArrayLike = DEFAULT_GAMMA, R: ArrayLike = DEFAULT_R):
    """Speed of sound for a calorically perfect gas: a = sqrt(gamma * R * T).

    Requires T > 0 K. Returns m/s.
    """
    T_arr = _validate_positive_T(T)
    gamma_arr = _validate_gamma(gamma)
    R_arr = _validate_R(R)
    return np.sqrt(gamma_arr * R_arr * T_arr)


def velocity_from_mach(
    M: ArrayLike, T: ArrayLike, gamma: ArrayLike = DEFAULT_GAMMA, R: ArrayLike = DEFAULT_R
):
    """Flow velocity from Mach number and static temperature: V = M * a.

    Requires M >= 0 (M = 0 is a valid stagnation state) and T > 0 K.
    Returns m/s.
    """
    M_arr = _validate_nonnegative_mach(M)
    a = speed_of_sound(T, gamma, R)
    return M_arr * a


def density_ideal_gas(p: ArrayLike, T: ArrayLike, R: ArrayLike = DEFAULT_R):
    """Ideal-gas density: rho = p / (R * T).

    Requires p > 0 Pa and T > 0 K. Returns kg/m^3.
    """
    p_arr = _validate_positive_p(p)
    T_arr = _validate_positive_T(T)
    R_arr = _validate_R(R)
    return p_arr / (R_arr * T_arr)


def T0_over_T(M: ArrayLike, gamma: ArrayLike = DEFAULT_GAMMA):
    """Isentropic stagnation-to-static temperature ratio.

    T0/T = 1 + (gamma - 1)/2 * M^2

    Valid for any M >= 0 (including M = 0, where T0/T = 1). This relation
    follows from the definition of stagnation temperature via the energy
    equation for a calorically perfect gas and holds for any locally
    adiabatic flow (isentropic or not), but is used here in the isentropic
    context consistent with p0/p and rho0/rho below.
    """
    M_arr = _validate_nonnegative_mach(M)
    gamma_arr = _validate_gamma(gamma)
    return 1.0 + 0.5 * (gamma_arr - 1.0) * M_arr**2


def p0_over_p(M: ArrayLike, gamma: ArrayLike = DEFAULT_GAMMA):
    """Isentropic stagnation-to-static pressure ratio.

    p0/p = [1 + (gamma - 1)/2 * M^2]^(gamma/(gamma - 1))

    Valid for any M >= 0. Requires isentropic (reversible, adiabatic) flow.
    """
    gamma_arr = _validate_gamma(gamma)
    ratio_T = T0_over_T(M, gamma_arr)
    return ratio_T ** (gamma_arr / (gamma_arr - 1.0))


def rho0_over_rho(M: ArrayLike, gamma: ArrayLike = DEFAULT_GAMMA):
    """Isentropic stagnation-to-static density ratio.

    rho0/rho = [1 + (gamma - 1)/2 * M^2]^(1/(gamma - 1))

    Valid for any M >= 0. Requires isentropic (reversible, adiabatic) flow.
    """
    gamma_arr = _validate_gamma(gamma)
    ratio_T = T0_over_T(M, gamma_arr)
    return ratio_T ** (1.0 / (gamma_arr - 1.0))
