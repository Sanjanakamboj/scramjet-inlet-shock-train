"""Independent verification of normal_shock.py.

Reference values are computed from hand-written reference formulas (not by
re-calling the production `normal_shock` function), per the project's
independent-verification requirement.
"""

import numpy as np
import pytest

from scramjet_inlet.gas_dynamics import density_ideal_gas, p0_over_p, speed_of_sound
from scramjet_inlet.normal_shock import normal_shock

GAMMA = 1.4


def _ref_M2(M1, g=GAMMA):
    return np.sqrt((1.0 + 0.5 * (g - 1.0) * M1**2) / (g * M1**2 - 0.5 * (g - 1.0)))


def _ref_p2_p1(M1, g=GAMMA):
    return 1.0 + (2.0 * g / (g + 1.0)) * (M1**2 - 1.0)


def _ref_rho2_rho1(M1, g=GAMMA):
    return ((g + 1.0) * M1**2) / ((g - 1.0) * M1**2 + 2.0)


def _ref_T2_T1(M1, g=GAMMA):
    return _ref_p2_p1(M1, g) / _ref_rho2_rho1(M1, g)


def _ref_p02_p01(M1, g=GAMMA):
    # Standard exact perfect-gas normal-shock stagnation-pressure ratio,
    # written independently from the production formula.
    term1 = (((g + 1.0) * M1**2) / ((g - 1.0) * M1**2 + 2.0)) ** (g / (g - 1.0))
    term2 = ((g + 1.0) / (2.0 * g * M1**2 - (g - 1.0))) ** (1.0 / (g - 1.0))
    return term1 * term2


M1_SAMPLE = 6.0


# 8. Normal-shock M2
def test_M2_hand_formula():
    result = normal_shock(M1_SAMPLE, GAMMA)
    assert result.M2 == pytest.approx(_ref_M2(M1_SAMPLE), rel=1e-12)


# 9. p2/p1
def test_p2_p1_hand_formula():
    result = normal_shock(M1_SAMPLE, GAMMA)
    assert result.p2_p1 == pytest.approx(_ref_p2_p1(M1_SAMPLE), rel=1e-12)


# 10. rho2/rho1
def test_rho2_rho1_hand_formula():
    result = normal_shock(M1_SAMPLE, GAMMA)
    assert result.rho2_rho1 == pytest.approx(_ref_rho2_rho1(M1_SAMPLE), rel=1e-12)


# 11. T2/T1
def test_T2_T1_hand_formula():
    result = normal_shock(M1_SAMPLE, GAMMA)
    assert result.T2_T1 == pytest.approx(_ref_T2_T1(M1_SAMPLE), rel=1e-12)


# 12. Identity p2/p1 = (rho2/rho1)(T2/T1)
def test_pressure_ratio_identity():
    result = normal_shock(M1_SAMPLE, GAMMA)
    assert result.p2_p1 == pytest.approx(result.rho2_rho1 * result.T2_T1, rel=1e-10)


# 13. M2 < 1 for representative supersonic inputs
@pytest.mark.parametrize("M1", [1.05, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0])
def test_M2_subsonic(M1):
    result = normal_shock(M1, GAMMA)
    assert result.M2 < 1.0


# 14, 15, 16: p2>p1, T2>T1, rho2>rho1
@pytest.mark.parametrize("M1", [1.05, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0])
def test_static_quantities_increase(M1):
    result = normal_shock(M1, GAMMA)
    assert result.p2_p1 > 1.0
    assert result.T2_T1 > 1.0
    assert result.rho2_rho1 > 1.0


# 17. 0 < p02/p01 < 1
@pytest.mark.parametrize("M1", [1.05, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0])
def test_stagnation_pressure_recovery_bounds(M1):
    result = normal_shock(M1, GAMMA)
    assert 0.0 < result.p02_p01 < 1.0


def test_p02_p01_hand_formula():
    result = normal_shock(M1_SAMPLE, GAMMA)
    assert result.p02_p01 == pytest.approx(_ref_p02_p01(M1_SAMPLE), rel=1e-10)


# 18. T0,2/T0,1 approx 1 (stagnation temperature conservation, adiabatic shock)
def test_stagnation_temperature_conserved():
    M1 = M1_SAMPLE
    T1 = 220.0
    result = normal_shock(M1, GAMMA)
    T2 = result.T2_T1 * T1
    T0_1 = T1 * (1.0 + 0.5 * (GAMMA - 1.0) * M1**2)
    T0_2 = T2 * (1.0 + 0.5 * (GAMMA - 1.0) * result.M2**2)
    assert T0_2 / T0_1 == pytest.approx(1.0, rel=1e-9)


# 19. Direct stagnation-pressure recovery agrees with reconstruction from
# downstream static state and M2 (fully independent reconstruction path:
# upstream p1,M1 -> p0,1 via isentropic relation; downstream p2 (from shock
# static ratio) and M2 -> p0,2 via isentropic relation; compare ratio).
def test_reconstructed_recovery_matches_direct():
    M1 = M1_SAMPLE
    p1 = 2500.0

    result = normal_shock(M1, GAMMA)

    # Reconstruction path (independent of the production p02_p01 formula):
    p0_1 = p1 * p0_over_p(M1, GAMMA)
    p2 = p1 * result.p2_p1
    p0_2 = p2 * p0_over_p(result.M2, GAMMA)
    reconstructed_p02_p01 = p0_2 / p0_1

    assert reconstructed_p02_p01 == pytest.approx(result.p02_p01, rel=1e-9)


# 20. Weak-shock limit M1 -> 1+
def test_weak_shock_limit():
    M1 = 1.0 + 1e-6
    result = normal_shock(M1, GAMMA)
    assert result.M2 == pytest.approx(1.0, abs=1e-4)
    assert result.p2_p1 == pytest.approx(1.0, abs=1e-4)
    assert result.rho2_rho1 == pytest.approx(1.0, abs=1e-4)
    assert result.T2_T1 == pytest.approx(1.0, abs=1e-4)
    assert result.p02_p01 == pytest.approx(1.0, abs=1e-6)


# 21. Increasing shock strength monotonically worsens stagnation-pressure
# recovery over the tested range.
def test_recovery_monotonically_decreases_with_M1():
    M1_values = np.array([1.05, 1.2, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0])
    recoveries = np.array([normal_shock(m, GAMMA).p02_p01 for m in M1_values])
    assert np.all(np.diff(recoveries) < 0.0)


# 22. Static pressure rise increases with shock strength over tested range.
def test_static_pressure_rise_increases_with_M1():
    M1_values = np.array([1.05, 1.2, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0])
    p2_p1_values = np.array([normal_shock(m, GAMMA).p2_p1 for m in M1_values])
    assert np.all(np.diff(p2_p1_values) > 0.0)


# 23. Scalar/vector consistency
def test_scalar_vector_consistency():
    M1_values = np.array([1.5, 2.0, 3.0, 6.0, 8.0])
    vec_result = normal_shock(M1_values, GAMMA)
    scalar_M2 = np.array([normal_shock(float(m), GAMMA).M2 for m in M1_values])
    assert np.allclose(vec_result.M2, scalar_M2, rtol=1e-12)


# 24/25. Invalid-input rejection: M1 <= 1
@pytest.mark.parametrize("M1_bad", [1.0, 0.9, 0.5, 0.0, -2.0])
def test_normal_shock_rejects_subsonic_or_sonic_M1(M1_bad):
    with pytest.raises(ValueError):
        normal_shock(M1_bad, GAMMA)


# 26. No NaN/inf over M1 in [1.05, 8]
def test_no_nan_inf_over_shock_sweep():
    M1_values = np.linspace(1.05, 8.0, 300)
    result = normal_shock(M1_values, GAMMA)
    for arr in (result.M2, result.p2_p1, result.rho2_rho1, result.T2_T1, result.p02_p01):
        assert np.all(np.isfinite(arr))


# Strong-shock density-ratio limit: (gamma+1)/(gamma-1) = 6 for gamma=1.4
def test_strong_shock_density_limit_trend():
    limit = (GAMMA + 1.0) / (GAMMA - 1.0)
    assert limit == pytest.approx(6.0, rel=1e-12)
    r_at_8 = normal_shock(8.0, GAMMA).rho2_rho1
    r_at_20 = normal_shock(20.0, GAMMA).rho2_rho1
    # Density ratio must stay below the limiting value and approach it
    # monotonically as M1 grows.
    assert r_at_8 < limit
    assert r_at_20 < limit
    assert r_at_20 > r_at_8


def test_representative_M1_6_full_state_end_to_end():
    """Independent end-to-end reconstruction of the representative
    hypersonic operating point (M1=6, T1=220 K, p1=2500 Pa) using only
    speed_of_sound/density_ideal_gas/p0_over_p and the hand reference
    formulas above -- not the production normal_shock() internals."""
    M1, T1, p1 = 6.0, 220.0, 2500.0
    g, R = 1.4, 287.05

    a1 = speed_of_sound(T1, g, R)
    assert a1 == pytest.approx(np.sqrt(g * R * T1), rel=1e-12)

    rho1 = density_ideal_gas(p1, T1, R)
    assert rho1 == pytest.approx(p1 / (R * T1), rel=1e-12)

    M2 = _ref_M2(M1, g)
    p2 = p1 * _ref_p2_p1(M1, g)
    T2 = T1 * _ref_T2_T1(M1, g)
    rho2 = rho1 * _ref_rho2_rho1(M1, g)

    result = normal_shock(M1, g)
    assert M2 == pytest.approx(result.M2, rel=1e-10)
    assert p2 == pytest.approx(p1 * result.p2_p1, rel=1e-10)
    assert T2 == pytest.approx(T1 * result.T2_T1, rel=1e-10)
    assert rho2 == pytest.approx(rho1 * result.rho2_rho1, rel=1e-10)

    # Ideal-gas law must still hold downstream
    assert p2 == pytest.approx(rho2 * R * T2, rel=1e-9)
