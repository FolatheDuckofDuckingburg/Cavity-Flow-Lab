import numpy as np
import pytest

from cavity_flow.solver import CavityConfig, CavitySolver


def test_projection_reduces_divergence_to_poisson_tolerance():
    rng = np.random.default_rng(7)
    solver = CavitySolver(CavityConfig(nx=8, ny=8, pressure_tolerance=1e-10,
                                       pressure_max_iterations=10000))
    u = rng.normal(size=solver.u.shape)
    v = rng.normal(size=solver.v.shape)
    solver.apply_boundaries_to(u, v)
    before = solver.divergence(u, v)
    phi, _, residual = solver._poisson(before)
    uc, vc = u.copy(), v.copy()
    uc[:, 1:-1] -= (phi[:, 1:] - phi[:, :-1]) / solver.dx
    vc[1:-1, :] -= (phi[1:, :] - phi[:-1, :]) / solver.dy
    solver.apply_boundaries_to(uc, vc)
    after = solver.divergence(uc, vc)
    assert residual <= solver.cfg.pressure_tolerance
    assert np.max(np.abs(after)) < 1e-8
    assert np.max(np.abs(after)) < np.max(np.abs(before)) * 1e-6


def test_poisson_solution_is_finite_and_zero_mean():
    solver = CavitySolver(CavityConfig(nx=8, ny=8, pressure_tolerance=1e-9,
                                       pressure_max_iterations=10000))
    rhs = np.zeros((8, 8))
    rhs[2, 3] = 1.0
    rhs[5, 4] = -1.0
    phi, _, residual = solver._poisson(rhs)
    assert np.isfinite(phi).all()
    assert abs(float(np.mean(phi))) < 1e-12
    assert residual <= solver.cfg.pressure_tolerance


def test_projection_rejects_wrong_rhs_shape():
    solver = CavitySolver(CavityConfig(nx=8, ny=8))
    with pytest.raises(ValueError):
        solver._poisson(np.zeros((7, 8)))
