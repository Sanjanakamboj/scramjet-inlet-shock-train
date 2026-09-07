"""Independent verification of gas_dynamics.py.

Reference values/formulas are written directly in these tests (not by
calling the production functions a second time), per the project's
independent-verification requirement.
"""

import numpy as np
import pytest

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

GAMMA = 1.4
R = 287.05


# 1. cp from gamma and R
def test_cp_hand_formula():
    # cp = gamma R / (gamma - 1); hand value for air at gamma=1.4, R=287.05
    expected = 1.4 * 287.05 / 0.4  # = 1004.675
    assert cp_from_gamma_R(1.4, 287.05) == pytest.approx(expected, rel=1e-12)
    assert cp_from_gamma_R(1.4, 287.05) == pytest.approx(1004.675, rel=1e-9)


# 2. Speed of sound
def test_speed_of_sound_hand_value():
    # a = sqrt(gamma R T); at T=220 K: sqrt(1.4*287.05*220)
    T = 220.0
    expected = np.sqrt(1.4 * 287.05 * 220.0)
    assert speed_of_sound(T, GAMMA, R) == pytest.approx(expected, rel=1e-12)
    assert speed_of_sound(288.15, GAMMA, R) == pytest.approx(340.29, rel=1e-3)


def test_speed_of_sound_rejects_nonpositive_T():
    with pytest.raises(ValueError):
        speed_of_sound(0.0)
    with pytest.raises(ValueError):
        speed_of_sound(-10.0)


# 3. Mach-to-velocity conversion
def test_velocity_from_mach_hand_value():
    M, T = 6.0, 220.0
    a_hand = np.sqrt(1.4 * 287.05 * 220.0)
    expected = M * a_hand
    assert velocity_from_mach(M, T, GAMMA, R) == pytest.approx(expected, rel=1e-12)


def test_velocity_from_mach_zero_allowed():
    assert velocity_from_mach(0.0, 300.0) == pytest.approx(0.0, abs=1e-12)


def test_velocity_from_mach_rejects_negative():
    with pytest.raises(ValueError):
        velocity_from_mach(-1.0, 300.0)


# 4. Ideal-gas density
def test_density_ideal_gas_hand_value():
    p, T = 2500.0, 220.0
    expected = p / (287.05 * T)
    assert density_ideal_gas(p, T, R) == pytest.approx(expected, rel=1e-12)


def test_density_rejects_nonpositive():
    with pytest.raises(ValueError):
        density_ideal_gas(0.0, 300.0)
    with pytest.raises(ValueError):
        density_ideal_gas(101325.0, 0.0)


# 5. Isentropic T0/T
def test_T0_over_T_hand_formula():
    M = 6.0
    expected = 1.0 + 0.5 * (GAMMA - 1.0) * M**2  # = 1 + 0.2*36 = 8.2
    assert T0_over_T(M, GAMMA) == pytest.approx(expected, rel=1e-12)
    assert T0_over_T(M, GAMMA) == pytest.approx(8.2, rel=1e-12)


def test_T0_over_T_at_M0_is_one():
    assert T0_over_T(0.0, GAMMA) == pytest.approx(1.0, abs=1e-12)


# 6. Isentropic p0/p
def test_p0_over_p_hand_formula():
    M = 6.0
    T_ratio = 1.0 + 0.5 * (GAMMA - 1.0) * M**2  # 8.2
    expected = T_ratio ** (GAMMA / (GAMMA - 1.0))  # 8.2^3.5
    assert p0_over_p(M, GAMMA) == pytest.approx(expected, rel=1e-12)


# 7. Isentropic rho0/rho
def test_rho0_over_rho_hand_formula():
    M = 6.0
    T_ratio = 1.0 + 0.5 * (GAMMA - 1.0) * M**2
    expected = T_ratio ** (1.0 / (GAMMA - 1.0))  # 8.2^2.5
    assert rho0_over_rho(M, GAMMA) == pytest.approx(expected, rel=1e-12)


def test_stagnation_ratios_consistent_with_ideal_gas_law():
    # p0/(rho0 R) = T0 must equal T0 computed from T0/T * T, i.e.
    # (p0/p)/(rho0/rho) == T0/T identically for the perfect-gas isentropic
    # relations (independent cross-check of internal consistency).
    M = 3.7
    ratio_check = p0_over_p(M, GAMMA) / rho0_over_rho(M, GAMMA)
    assert ratio_check == pytest.approx(T0_over_T(M, GAMMA), rel=1e-10)


# 23. Scalar/vector consistency
def test_scalar_vector_consistency():
    Ms = np.array([1.5, 2.0, 3.0, 6.0, 8.0])
    vec_result = T0_over_T(Ms, GAMMA)
    scalar_results = np.array([T0_over_T(float(m), GAMMA) for m in Ms])
    assert np.allclose(vec_result, scalar_results, rtol=1e-12)


# 24. Invalid-input rejection
@pytest.mark.parametrize("gamma", [1.0, 0.9, 0.0, -1.4])
def test_invalid_gamma_rejected(gamma):
    with pytest.raises(ValueError):
        cp_from_gamma_R(gamma, DEFAULT_R)


@pytest.mark.parametrize("R_bad", [0.0, -287.05])
def test_invalid_R_rejected(R_bad):
    with pytest.raises(ValueError):
        cp_from_gamma_R(DEFAULT_GAMMA, R_bad)


def test_negative_mach_rejected_everywhere():
    for fn in (T0_over_T, p0_over_p, rho0_over_rho):
        with pytest.raises(ValueError):
            fn(-0.5, GAMMA)


# 26. No NaN/inf over representative sweep
def test_no_nan_inf_over_mach_sweep():
    Ms = np.linspace(0.0, 8.0, 200)
    for fn in (T0_over_T, p0_over_p, rho0_over_rho):
        vals = fn(Ms, GAMMA)
        assert np.all(np.isfinite(vals))
