"""Software invariants only; passing these tests is not CFD validation."""
import numpy as np
import pytest

from cavity_flow import CavityConfig, CavitySolver


@pytest.mark.parametrize("nx,ny", [(4, 4), (8, 6), (16, 16)])
def test_staggered_array_shapes(nx, ny):
    solver = CavitySolver(CavityConfig(nx=nx, ny=ny))
    assert solver.u.shape == (ny, nx + 1)
    assert solver.v.shape == (ny + 1, nx)
    assert solver.p.shape == (ny, nx)


def test_initial_state_is_finite_and_divergence_free():
    solver = CavitySolver(CavityConfig(nx=8, ny=8))
    assert np.isfinite(solver.u).all()
    assert np.isfinite(solver.v).all()
    assert np.isfinite(solver.p).all()
    assert solver.divergence_linf() == pytest.approx(0.0, abs=1e-14)


def test_wall_velocity_values_excluding_corner_conflicts():
    cfg = CavityConfig(nx=8, ny=8, lid_speed=1.25)
    solver = CavitySolver(cfg)
    assert np.allclose(solver.u[1:-1, 0], 0.0)
    assert np.allclose(solver.u[1:-1, -1], 0.0)
    assert np.allclose(solver.u[0, :], 0.0)
    assert np.allclose(solver.u[-1, 1:-1], cfg.lid_speed)
    assert np.allclose(solver.v[0, :], 0.0)
    assert np.allclose(solver.v[-1, :], 0.0)
    assert np.allclose(solver.v[:, 0], 0.0)
    assert np.allclose(solver.v[:, -1], 0.0)


@pytest.mark.parametrize("kwargs", [
    {"nx": 3}, {"ny": 2}, {"reynolds": 0}, {"length": 0},
    {"lid_speed": 0}, {"dt": 0}, {"dt": -0.1}, {"cfl": 0},
])
def test_invalid_configuration_is_rejected(kwargs):
    with pytest.raises(ValueError):
        CavityConfig(**kwargs)


def test_one_step_reports_pre_and_post_projection_diagnostics():
    solver = CavitySolver(CavityConfig(nx=8, ny=8, pressure_max_iterations=10000))
    row = solver.step()
    assert solver.iteration == 1
    assert solver.time > 0
    assert row["dt"] > 0
    assert row["divergence_pre_projection_linf"] >= 0
    assert row["divergence_post_projection_linf"] == pytest.approx(solver.divergence_linf())
    assert np.isfinite(solver.u).all() and np.isfinite(solver.v).all()
    assert np.isfinite(solver.p).all()
    assert np.isfinite(float(row["pressure_residual_linf"]))


def test_run_history_includes_final_step():
    solver = CavitySolver(CavityConfig(nx=8, ny=8))
    history = solver.run(steps=3, log_every=2)
    assert history[-1]["iteration"] == 3


def test_run_until_steady_respects_max_steps():
    solver = CavitySolver(CavityConfig(nx=8, ny=8))
    history = solver.run_until_steady(max_steps=3, min_steps=2, tolerance=1e-30)
    assert len(history) == 3
