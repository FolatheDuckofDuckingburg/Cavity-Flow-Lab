"""Software invariants only; these tests do not prove CFD validity."""
import numpy as np
import pytest
from cavity_flow import CavityConfig, CavitySolver


@pytest.mark.parametrize("nx,ny", [(4, 4), (8, 6), (16, 16)])
def test_staggered_array_shapes(nx, ny):
    s = CavitySolver(CavityConfig(nx=nx, ny=ny))
    assert s.u.shape == (ny, nx + 1)
    assert s.v.shape == (ny + 1, nx)
    assert s.p.shape == (ny, nx)


def test_initial_state_is_finite_and_divergence_free():
    s = CavitySolver(CavityConfig(nx=8, ny=8))
    assert np.isfinite(s.u).all() and np.isfinite(s.v).all()
    assert np.isfinite(s.p).all()
    assert s.divergence_linf() == pytest.approx(0.0, abs=1e-14)


def test_wall_velocity_values():
    cfg = CavityConfig(nx=8, ny=8, lid_speed=1.25)
    s = CavitySolver(cfg)
    assert np.allclose(s.u[:, 0], 0.0)
    assert np.allclose(s.u[:, -1], 0.0)
    assert np.allclose(s.u[0, :], 0.0)
    assert np.allclose(s.u[-1, :], cfg.lid_speed)
    assert np.allclose(s.v[0, :], 0.0) and np.allclose(s.v[-1, :], 0.0)
    assert np.allclose(s.v[:, 0], 0.0) and np.allclose(s.v[:, -1], 0.0)


def test_one_step_is_finite_and_advances_time():
    s = CavitySolver(CavityConfig(nx=8, ny=8))
    d = s.step()
    assert s.iteration == 1 and s.time > 0.0 and d["dt"] > 0.0
    assert np.isfinite(s.u).all() and np.isfinite(s.v).all() and np.isfinite(s.p).all()
    assert np.isfinite(d["divergence_linf"])
    assert np.isfinite(d["pressure_residual_linf"])
