import numpy as np

from cavity_flow.solver import CavityConfig, CavitySolver


def test_projection_reduces_divergence_to_poisson_tolerance():
    rng = np.random.default_rng(7)
    solver = CavitySolver(
        CavityConfig(
            nx=8,
            ny=8,
            pressure_tolerance=1e-10,
            pressure_max_iterations=10000,
        )
    )
    u = rng.normal(size=solver.u.shape)
    v = rng.normal(size=solver.v.shape)
    solver.apply_boundaries_to(u, v)

    divergence_before = solver.divergence(u, v)
    rhs = divergence_before.copy()  # dt = 1 for this operator-level test
    phi, _, residual = solver._poisson(rhs)

    u_corrected, v_corrected = u.copy(), v.copy()
    u_corrected[:, 1:-1] -= (phi[:, 1:] - phi[:, :-1]) / solver.dx
    v_corrected[1:-1, :] -= (phi[1:, :] - phi[:-1, :]) / solver.dy
    solver.apply_boundaries_to(u_corrected, v_corrected)

    divergence_after = solver.divergence(u_corrected, v_corrected)
    assert residual <= solver.cfg.pressure_tolerance
    assert np.max(np.abs(divergence_after)) < 1e-8
    assert np.max(np.abs(divergence_after)) < np.max(np.abs(divergence_before)) * 1e-6


def test_poisson_solution_is_gauge_invariant_and_finite():
    solver = CavitySolver(
        CavityConfig(nx=8, ny=8, pressure_tolerance=1e-9, pressure_max_iterations=10000)
    )
    rhs = np.zeros((8, 8))
    rhs[2, 3] = 1.0
    rhs[5, 4] = -1.0
    phi, _, residual = solver._poisson(rhs)

    assert np.isfinite(phi).all()
    assert abs(float(np.mean(phi))) < 1e-12
    assert residual <= solver.cfg.pressure_tolerance
