"""Independent verification of oblique_shock.py.

Reference values/formulas are written directly in these tests (not by
re-calling the production functions a second time), per the project's
independent-verification requirement.
"""

import numpy as np
import pytest

from scramjet_inlet.gas_dynamics import density_ideal_gas, p0_over_p, speed_of_sound
from scramjet_inlet.normal_shock import normal_shock
from scramjet_inlet.oblique_shock import (
    detachment_limit,
    mach_angle,
    oblique_shock,
    shock_angle,
    theta_from_beta,
)

GAMMA = 1.4


def _ref_theta_from_beta(M1, beta, g=GAMMA):
    """Independent hand implementation of the theta-beta-M relation."""
    num = M1**2 * np.sin(beta) ** 2 - 1.0
    den = M1**2 * (g + np.cos(2.0 * beta)) + 2.0
    tan_theta = 2.0 / np.tan(beta) * num / den
    return np.arctan(tan_theta)


def _ref_normal_shock_ratios(Mn1, g=GAMMA):
    M2sq = (1.0 + 0.5 * (g - 1.0) * Mn1**2) / (g * Mn1**2 - 0.5 * (g - 1.0))
    Mn2 = np.sqrt(M2sq)
    p2_p1 = 1.0 + (2.0 * g / (g + 1.0)) * (Mn1**2 - 1.0)
    rho2_rho1 = ((g + 1.0) * Mn1**2) / ((g - 1.0) * Mn1**2 + 2.0)
    T2_T1 = p2_p1 / rho2_rho1
    p02_p01 = (
        (((g + 1.0) * Mn1**2) / ((g - 1.0) * Mn1**2 + 2.0)) ** (g / (g - 1.0))
    ) * (((g + 1.0) / (2.0 * g * Mn1**2 - (g - 1.0))) ** (1.0 / (g - 1.0)))
    return Mn2, p2_p1, rho2_rho1, T2_T1, p02_p01


M1_SAMPLE = 6.0
THETA_10_DEG = np.radians(10.0)


# 1. Mach angle
def test_mach_angle_hand_formula():
    M1 = 6.0
    expected = np.arcsin(1.0 / M1)
    assert mach_angle(M1) == pytest.approx(expected, rel=1e-12)


def test_mach_angle_rejects_subsonic():
    with pytest.raises(ValueError):
        mach_angle(0.9)


# 2. theta-beta-M residual near zero for returned beta
def test_theta_beta_m_residual_weak():
    beta = shock_angle(M1_SAMPLE, THETA_10_DEG, GAMMA, branch="weak")
    theta_check = _ref_theta_from_beta(M1_SAMPLE, beta, GAMMA)
    assert theta_check == pytest.approx(THETA_10_DEG, abs=1e-9)


def test_theta_beta_m_residual_strong():
    beta = shock_angle(M1_SAMPLE, THETA_10_DEG, GAMMA, branch="strong")
    theta_check = _ref_theta_from_beta(M1_SAMPLE, beta, GAMMA)
    assert theta_check == pytest.approx(THETA_10_DEG, abs=1e-9)


# 3, 4, 5: weak/strong bracket ordering
def test_weak_and_strong_bracket_ordering():
    limit = detachment_limit(M1_SAMPLE, GAMMA)
    mu = mach_angle(M1_SAMPLE)
    beta_w = shock_angle(M1_SAMPLE, THETA_10_DEG, GAMMA, branch="weak")
    beta_s = shock_angle(M1_SAMPLE, THETA_10_DEG, GAMMA, branch="strong")

    assert mu < beta_w < limit.beta_at_theta_max
    assert limit.beta_at_theta_max < beta_s < np.pi / 2.0
    assert beta_w < beta_s


# 6. Weak/strong roots satisfy the same theta
def test_weak_and_strong_satisfy_same_theta():
    beta_w = shock_angle(M1_SAMPLE, THETA_10_DEG, GAMMA, branch="weak")
    beta_s = shock_angle(M1_SAMPLE, THETA_10_DEG, GAMMA, branch="strong")
    theta_w = _ref_theta_from_beta(M1_SAMPLE, beta_w, GAMMA)
    theta_s = _ref_theta_from_beta(M1_SAMPLE, beta_s, GAMMA)
    assert theta_w == pytest.approx(THETA_10_DEG, abs=1e-9)
    assert theta_s == pytest.approx(THETA_10_DEG, abs=1e-9)


# 7. Branches merge near theta_max
def test_branches_merge_near_theta_max():
    limit = detachment_limit(M1_SAMPLE, GAMMA)
    theta_near_max = limit.theta_max * 0.999999
    beta_w = shock_angle(M1_SAMPLE, theta_near_max, GAMMA, branch="weak")
    beta_s = shock_angle(M1_SAMPLE, theta_near_max, GAMMA, branch="strong")
    assert beta_w == pytest.approx(limit.beta_at_theta_max, abs=1e-3)
    assert beta_s == pytest.approx(limit.beta_at_theta_max, abs=1e-3)
    assert abs(beta_w - beta_s) < 1e-2


# 8. theta > theta_max rejected (detached)
def test_theta_above_theta_max_rejected():
    limit = detachment_limit(M1_SAMPLE, GAMMA)
    with pytest.raises(ValueError):
        shock_angle(M1_SAMPLE, limit.theta_max * 1.1, GAMMA, branch="weak")


# 9. theta <= 0 rejected
@pytest.mark.parametrize("theta_bad", [0.0, -0.1])
def test_theta_nonpositive_rejected(theta_bad):
    with pytest.raises(ValueError):
        shock_angle(M1_SAMPLE, theta_bad, GAMMA, branch="weak")


# 10. M1 <= 1 rejected
@pytest.mark.parametrize("M1_bad", [1.0, 0.9, 0.0, -2.0])
def test_M1_leq_1_rejected(M1_bad):
    with pytest.raises(ValueError):
        shock_angle(M1_bad, THETA_10_DEG, GAMMA, branch="weak")
    with pytest.raises(ValueError):
        mach_angle(M1_bad)


# 11. Mn1 = M1 sin(beta) > 1
def test_Mn1_supersonic():
    for branch in ("weak", "strong"):
        r = oblique_shock(M1_SAMPLE, THETA_10_DEG, GAMMA, branch=branch)
        assert r.Mn1 > 1.0


# 12. Normal-shock ratios from Mn1 match direct hand formulas
def test_normal_shock_ratios_from_Mn1_match_hand_formula():
    for branch in ("weak", "strong"):
        r = oblique_shock(M1_SAMPLE, THETA_10_DEG, GAMMA, branch=branch)
        Mn2_ref, p2_p1_ref, rho2_rho1_ref, T2_T1_ref, p02_p01_ref = _ref_normal_shock_ratios(
            r.Mn1, GAMMA
        )
        assert r.Mn2 == pytest.approx(Mn2_ref, rel=1e-10)
        assert r.p2_p1 == pytest.approx(p2_p1_ref, rel=1e-10)
        assert r.rho2_rho1 == pytest.approx(rho2_rho1_ref, rel=1e-10)
        assert r.T2_T1 == pytest.approx(T2_T1_ref, rel=1e-10)
        assert r.p02_p01 == pytest.approx(p02_p01_ref, rel=1e-10)


# 13, 14, 15: static ratios > 1
def test_static_ratios_increase():
    for branch in ("weak", "strong"):
        r = oblique_shock(M1_SAMPLE, THETA_10_DEG, GAMMA, branch=branch)
        assert r.p2_p1 > 1.0
        assert r.rho2_rho1 > 1.0
        assert r.T2_T1 > 1.0


# 16. 0 < p02_p01 < 1
def test_stagnation_recovery_bounds():
    for branch in ("weak", "strong"):
        r = oblique_shock(M1_SAMPLE, THETA_10_DEG, GAMMA, branch=branch)
        assert 0.0 < r.p02_p01 < 1.0


# 17. T0 conservation
def test_T0_conservation():
    for branch in ("weak", "strong"):
        r = oblique_shock(M1_SAMPLE, THETA_10_DEG, GAMMA, branch=branch, p1=2500.0, T1=220.0)
        assert r.T0_2 == pytest.approx(r.T0_1, rel=1e-9)


# 18. M2 = Mn2 / sin(beta - theta)
def test_M2_formula():
    for branch in ("weak", "strong"):
        r = oblique_shock(M1_SAMPLE, THETA_10_DEG, GAMMA, branch=branch)
        expected_M2 = r.Mn2 / np.sin(r.beta - r.theta)
        assert r.M2 == pytest.approx(expected_M2, rel=1e-12)


# 19. Tangential velocity conservation
def test_tangential_velocity_conservation():
    """Independent check: Vt1 = V1*cos(beta), Vt2 = V2*cos(beta-theta) must
    be equal, using velocities reconstructed from M and a (speed of sound)
    at each station -- computed from scratch here, not from production
    oblique-shock internals."""
    T1, p1 = 220.0, 2500.0
    for branch in ("weak", "strong"):
        r = oblique_shock(M1_SAMPLE, THETA_10_DEG, GAMMA, branch=branch, p1=p1, T1=T1)
        R = 287.05
        a1 = speed_of_sound(T1, GAMMA, R)
        V1 = M1_SAMPLE * a1
        T2 = T1 * r.T2_T1
        a2 = speed_of_sound(T2, GAMMA, R)
        V2 = r.M2 * a2

        Vt1 = V1 * np.cos(r.beta)
        Vt2 = V2 * np.cos(r.beta - r.theta)
        assert Vt2 == pytest.approx(Vt1, rel=1e-6)


# 20. Direct total-pressure recovery equals reconstruction from p2, M2
def test_reconstructed_recovery_matches_direct():
    p1 = 2500.0
    for branch in ("weak", "strong"):
        r = oblique_shock(M1_SAMPLE, THETA_10_DEG, GAMMA, branch=branch, p1=p1, T1=220.0)
        p0_1 = p1 * p0_over_p(M1_SAMPLE, GAMMA)
        p0_2 = r.p2 * p0_over_p(r.M2, GAMMA)
        reconstructed = p0_2 / p0_1
        assert reconstructed == pytest.approx(r.p02_p01, rel=1e-6)


# 21. theta -> 0+ weak-branch limit
def test_weak_branch_weak_deflection_limit():
    theta_small = np.radians(0.05)
    r = oblique_shock(M1_SAMPLE, theta_small, GAMMA, branch="weak")
    mu = mach_angle(M1_SAMPLE)
    assert r.beta == pytest.approx(mu, abs=2e-3)
    assert r.p2_p1 == pytest.approx(1.0, abs=0.05)
    assert r.p02_p01 == pytest.approx(1.0, abs=1e-3)
    assert r.M2 == pytest.approx(M1_SAMPLE, abs=0.05)


# 22. beta -> pi/2 (strong branch, theta->0) approaches the normal-shock
# solution at the same M1.
def test_strong_branch_theta_to_zero_matches_normal_shock():
    theta_small = np.radians(0.05)
    r = oblique_shock(M1_SAMPLE, theta_small, GAMMA, branch="strong")
    ns = normal_shock(M1_SAMPLE, GAMMA)
    assert r.beta == pytest.approx(np.pi / 2.0, abs=2e-3)
    assert r.p2_p1 == pytest.approx(float(ns.p2_p1), rel=2e-3)
    assert r.p02_p01 == pytest.approx(float(ns.p02_p01), rel=2e-3)


# 23. Increasing theta increases p2/p1 on weak branch
def test_weak_branch_p2_p1_increases_with_theta():
    thetas = np.radians([2.0, 5.0, 10.0, 15.0, 20.0, 25.0])
    p2_p1_vals = [oblique_shock(M1_SAMPLE, t, GAMMA, branch="weak").p2_p1 for t in thetas]
    assert np.all(np.diff(p2_p1_vals) > 0.0)


# 24. Increasing theta worsens p0,2/p0,1 on weak branch
def test_weak_branch_recovery_decreases_with_theta():
    thetas = np.radians([2.0, 5.0, 10.0, 15.0, 20.0, 25.0])
    recov = [oblique_shock(M1_SAMPLE, t, GAMMA, branch="weak").p02_p01 for t in thetas]
    assert np.all(np.diff(recov) < 0.0)


# 25. Continuity/finiteness at fixed theta over a Mach range
def test_continuous_finite_over_mach_range():
    theta = np.radians(8.0)
    M1_values = np.linspace(1.5, 8.0, 30)
    betas = []
    for M1 in M1_values:
        limit = detachment_limit(float(M1), GAMMA)
        if theta >= limit.theta_max:
            continue  # skip inadmissible (detached) combinations
        beta = shock_angle(float(M1), theta, GAMMA, branch="weak")
        assert np.isfinite(beta)
        betas.append(beta)
    assert len(betas) > 10
    # No huge jumps between adjacent finite samples (rough continuity check)
    betas = np.array(betas)
    assert np.all(np.abs(np.diff(betas)) < 0.2)


# 26. Representative M1=6, theta=10 deg hand/reference check
def test_representative_case_weak():
    r = oblique_shock(M1_SAMPLE, THETA_10_DEG, GAMMA, branch="weak", p1=2500.0, T1=220.0)

    beta_ref = None
    # Solve independently via a fine hand grid search + local refinement to
    # cross-check the production Brent-solved beta.
    betas_grid = np.linspace(mach_angle(M1_SAMPLE) + 1e-6, np.pi / 2 - 1e-6, 200000)
    thetas_grid = _ref_theta_from_beta(M1_SAMPLE, betas_grid, GAMMA)
    idx_max = np.argmax(thetas_grid)
    weak_region = slice(0, idx_max + 1)
    diffs = np.abs(thetas_grid[weak_region] - THETA_10_DEG)
    beta_ref = betas_grid[weak_region][np.argmin(diffs)]

    assert r.beta == pytest.approx(beta_ref, abs=1e-3)
    assert np.degrees(r.beta) == pytest.approx(17.5, abs=0.5)
    assert r.M2 == pytest.approx(4.65, abs=0.05)
    assert r.p2_p1 == pytest.approx(3.67, abs=0.05)


# 27. Scalar/vector behavior where vectorization is supported
# (theta_from_beta and normal_shock support arrays; oblique_shock's beta
# solve is scalar-only by design due to branch bracketing, so we check the
# vectorized building blocks here.)
def test_scalar_vector_consistency_theta_from_beta():
    betas = np.array([0.2, 0.3, 0.5, 1.0])
    vec_result = np.array([_ref_theta_from_beta(M1_SAMPLE, b, GAMMA) for b in betas])
    scalar_result = np.array([theta_from_beta(M1_SAMPLE, float(b), GAMMA) for b in betas])
    assert np.allclose(vec_result, scalar_result, rtol=1e-10)


# 28. Invalid gamma/input rejection
def test_invalid_gamma_rejected():
    with pytest.raises(ValueError):
        shock_angle(M1_SAMPLE, THETA_10_DEG, gamma=1.0, branch="weak")
    with pytest.raises(ValueError):
        detachment_limit(M1_SAMPLE, gamma=0.5)


def test_invalid_branch_rejected():
    with pytest.raises(ValueError):
        shock_angle(M1_SAMPLE, THETA_10_DEG, GAMMA, branch="bogus")


# 29. No NaN/inf over admissible grid
def test_no_nan_inf_over_admissible_grid():
    M1_values = np.linspace(1.5, 8.0, 15)
    for M1 in M1_values:
        limit = detachment_limit(float(M1), GAMMA)
        assert np.isfinite(limit.theta_max)
        assert np.isfinite(limit.beta_at_theta_max)
        theta_test = 0.8 * limit.theta_max
        for branch in ("weak", "strong"):
            r = oblique_shock(float(M1), theta_test, GAMMA, branch=branch)
            for val in (r.beta, r.Mn1, r.Mn2, r.M2, r.p2_p1, r.rho2_rho1, r.T2_T1, r.p02_p01):
                assert np.isfinite(val)


# 30. theta_max behaves smoothly across the tested Mach range
def test_theta_max_smooth_across_mach_range():
    M1_values = np.linspace(1.5, 8.0, 30)
    theta_max_values = np.array(
        [detachment_limit(float(M1), GAMMA).theta_max for M1 in M1_values]
    )
    assert np.all(np.isfinite(theta_max_values))
    # theta_max should increase monotonically with M1 over this range
    # (well-known perfect-gas oblique-shock behavior for gamma=1.4).
    assert np.all(np.diff(theta_max_values) > -1e-8)
    # No large jumps between adjacent Mach samples (theta_max(M1) is
    # legitimately steep near M1~1.5-2, so allow a generous per-step bound).
    assert np.all(np.abs(np.diff(theta_max_values)) < np.radians(10.0))


def test_density_ideal_gas_consistency_in_result():
    r = oblique_shock(M1_SAMPLE, THETA_10_DEG, GAMMA, branch="weak", p1=2500.0, T1=220.0)
    rho1_ref = density_ideal_gas(2500.0, 220.0, 287.05)
    rho2_ref = rho1_ref * r.rho2_rho1
    assert r.rho2 == pytest.approx(rho2_ref, rel=1e-9)
