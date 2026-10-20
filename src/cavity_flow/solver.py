"""Experimental MAC-grid lid-driven cavity solver.

This implementation is educational and not yet a validated production CFD code.
In particular, its wall treatment and spatial discretization require comparison
against independently verified reference solutions.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class CavityConfig:
    nx: int = 64
    ny: int = 64
    reynolds: float = 100.0
    lid_speed: float = 1.0
    length: float = 1.0
    dt: float | None = None
    cfl: float = 0.25
    pressure_tolerance: float = 1e-8
    pressure_max_iterations: int = 2000

    def __post_init__(self) -> None:
        if self.nx < 4 or self.ny < 4:
            raise ValueError("nx and ny must be >= 4")
        if not math.isfinite(self.reynolds) or self.reynolds <= 0:
            raise ValueError("reynolds must be finite and positive")
        if not math.isfinite(self.length) or self.length <= 0:
            raise ValueError("length must be finite and positive")
        if not math.isfinite(self.lid_speed) or self.lid_speed <= 0:
            raise ValueError("lid_speed must be finite and positive")
        if not math.isfinite(self.cfl) or not 0 < self.cfl <= 1:
            raise ValueError("cfl must be finite and in (0, 1]")
        if self.dt is not None and (not math.isfinite(self.dt) or self.dt <= 0):
            raise ValueError("dt must be finite and positive when provided")
        if (not math.isfinite(self.pressure_tolerance)
                or self.pressure_tolerance <= 0
                or self.pressure_max_iterations < 1):
            raise ValueError("pressure tolerance and iteration limit must be positive")


class CavitySolver:
    """2D cavity solver with staggered velocities and cell-centred pressure.

    Array convention:
      u: (ny, nx + 1), x-face velocity
      v: (ny + 1, nx), y-face velocity
      p: (ny, nx), cell-centred pressure
    """

    def __init__(self, config: CavityConfig):
        self.cfg = config
        self.dx = config.length / config.nx
        self.dy = config.length / config.ny
        self.nu = config.lid_speed * config.length / config.reynolds
        self.u = np.zeros((config.ny, config.nx + 1), dtype=float)
        self.v = np.zeros((config.ny + 1, config.nx), dtype=float)
        self.p = np.zeros((config.ny, config.nx), dtype=float)
        self.time = 0.0
        self.iteration = 0
        self.apply_velocity_boundaries()

    def apply_velocity_boundaries(self) -> None:
        self.apply_boundaries_to(self.u, self.v)

    def apply_boundaries_to(self, u: FloatArray, v: FloatArray) -> None:
        """Apply the current code's wall-value convention.

        The lid/side-wall corner convention is explicit: the lid assignment
        overwrites u's top-row corner values. This simplified treatment remains
        a known numerical limitation, not a claim of high-order wall accuracy.
        """
        u[:, 0] = u[:, -1] = 0.0
        u[0, :] = 0.0
        u[-1, :] = self.cfg.lid_speed
        v[0, :] = v[-1, :] = 0.0
        v[:, 0] = v[:, -1] = 0.0

    def divergence(self, u: FloatArray | None = None, v: FloatArray | None = None) -> FloatArray:
        """Finite-volume divergence at cell centres."""
        uu = self.u if u is None else u
        vv = self.v if v is None else v
        if uu.shape != self.u.shape or vv.shape != self.v.shape:
            raise ValueError("velocity arrays have shapes incompatible with this grid")
        return (uu[:, 1:] - uu[:, :-1]) / self.dx + (vv[1:, :] - vv[:-1, :]) / self.dy

    def divergence_linf(self) -> float:
        return float(np.max(np.abs(self.divergence())))

    def kinetic_energy(self) -> float:
        """Mean cell-centred kinetic energy per unit mass."""
        uc = 0.5 * (self.u[:, :-1] + self.u[:, 1:])
        vc = 0.5 * (self.v[:-1, :] + self.v[1:, :])
        return float(0.5 * np.mean(uc**2 + vc**2))

    def stable_dt(self) -> float:
        """Conservative explicit advection/diffusion time-step estimate."""
        umax = max(float(np.max(np.abs(self.u))), abs(self.cfg.lid_speed), 1e-12)
        vmax = max(float(np.max(np.abs(self.v))), 1e-12)
        dt_adv = self.cfg.cfl / (umax / self.dx + vmax / self.dy)
        dt_diff = 0.24 / (self.nu * (1.0 / self.dx**2 + 1.0 / self.dy**2))
        estimate = min(dt_adv, dt_diff)
        return min(self.cfg.dt, estimate) if self.cfg.dt is not None else estimate

    def _tendencies(self) -> tuple[FloatArray, FloatArray]:
        """Compute centred advection/diffusion predictor tendencies."""
        u, v, dx, dy, nu = self.u, self.v, self.dx, self.dy, self.nu
        du_dx = np.zeros_like(u)
        du_dy = np.zeros_like(u)
        lap_u = np.zeros_like(u)
        du_dx[:, 1:-1] = (u[:, 2:] - u[:, :-2]) / (2 * dx)
        du_dy[1:-1, :] = (u[2:, :] - u[:-2, :]) / (2 * dy)
        lap_u[:, 1:-1] += (u[:, 2:] - 2*u[:, 1:-1] + u[:, :-2]) / dx**2
        lap_u[1:-1, :] += (u[2:, :] - 2*u[1:-1, :] + u[:-2, :]) / dy**2
        vcc = 0.5 * (v[:-1, :] + v[1:, :])
        v_at_u = np.zeros_like(u)
        v_at_u[:, 1:-1] = 0.5 * (vcc[:, :-1] + vcc[:, 1:])
        fu = -(u * du_dx + v_at_u * du_dy) + nu * lap_u

        dv_dx = np.zeros_like(v)
        dv_dy = np.zeros_like(v)
        lap_v = np.zeros_like(v)
        dv_dx[:, 1:-1] = (v[:, 2:] - v[:, :-2]) / (2 * dx)
        dv_dy[1:-1, :] = (v[2:, :] - v[:-2, :]) / (2 * dy)
        lap_v[:, 1:-1] += (v[:, 2:] - 2*v[:, 1:-1] + v[:, :-2]) / dx**2
        lap_v[1:-1, :] += (v[2:, :] - 2*v[1:-1, :] + v[:-2, :]) / dy**2
        ucc = 0.5 * (u[:, :-1] + u[:, 1:])
        u_at_v = np.zeros_like(v)
        u_at_v[1:-1, :] = 0.5 * (ucc[:-1, :] + ucc[1:, :])
        fv = -(u_at_v * dv_dx + v * dv_dy) + nu * lap_v
        return fu, fv

    def _apply_pressure_laplacian(self, a: FloatArray) -> FloatArray:
        """Apply the cell-centred Laplacian with zero-normal-gradient ghosts."""
        east = np.empty_like(a)
        west = np.empty_like(a)
        north = np.empty_like(a)
        south = np.empty_like(a)
        east[:, :-1], east[:, -1] = a[:, 1:], a[:, -1]
        west[:, 1:], west[:, 0] = a[:, :-1], a[:, 0]
        north[:-1, :], north[-1, :] = a[1:, :], a[-1, :]
        south[1:, :], south[0, :] = a[:-1, :], a[0, :]
        return (east - 2*a + west) / self.dx**2 + (north - 2*a + south) / self.dy**2

    def _poisson(self, rhs: FloatArray) -> tuple[FloatArray, int, float]:
        """Jacobi solve for Laplacian(phi)=rhs with homogeneous Neumann walls.

        The right-hand side is projected to zero mean to satisfy the Neumann
        compatibility condition. The returned residual is computed from the
        discrete operator used by the solver.
        """
        if rhs.shape != (self.cfg.ny, self.cfg.nx):
            raise ValueError(f"rhs must have shape {(self.cfg.ny, self.cfg.nx)}, got {rhs.shape}")
        if not np.isfinite(rhs).all():
            raise ValueError("Poisson RHS must contain only finite values")
        rhs = np.asarray(rhs, dtype=float) - float(np.mean(rhs))
        phi = np.zeros_like(rhs)
        dx2, dy2 = self.dx**2, self.dy**2
        diagonal = np.full_like(rhs, 2.0 / dx2 + 2.0 / dy2)
        diagonal[:, 0] -= 1.0 / dx2
        diagonal[:, -1] -= 1.0 / dx2
        diagonal[0, :] -= 1.0 / dy2
        diagonal[-1, :] -= 1.0 / dy2

        residual = float("inf")
        for it in range(1, self.cfg.pressure_max_iterations + 1):
            weighted = np.zeros_like(phi)
            weighted[:, :-1] += phi[:, 1:] / dx2
            weighted[:, 1:] += phi[:, :-1] / dx2
            weighted[:-1, :] += phi[1:, :] / dy2
            weighted[1:, :] += phi[:-1, :] / dy2
            phi = (weighted - rhs) / diagonal
            phi -= float(np.mean(phi))
            if it == 1 or it % 10 == 0 or it == self.cfg.pressure_max_iterations:
                residual = float(np.max(np.abs(self._apply_pressure_laplacian(phi) - rhs)))
                if residual <= self.cfg.pressure_tolerance:
                    return phi, it, residual
        return phi, self.cfg.pressure_max_iterations, residual

    def step(self) -> dict[str, float | bool]:
        """Advance one predictor/projection step and return diagnostics."""
        dt = self.stable_dt()
        old_u, old_v = self.u.copy(), self.v.copy()
        fu, fv = self._tendencies()
        us, vs = self.u + dt*fu, self.v + dt*fv
        self.apply_boundaries_to(us, vs)
        divergence_before = self.divergence(us, vs)
        rhs = divergence_before / dt
        phi, pressure_iterations, pressure_residual = self._poisson(rhs)

        un, vn = us.copy(), vs.copy()
        un[:, 1:-1] -= dt * (phi[:, 1:] - phi[:, :-1]) / self.dx
        vn[1:-1, :] -= dt * (phi[1:, :] - phi[:-1, :]) / self.dy
        self.u, self.v = un, vn
        self.apply_velocity_boundaries()
        self.p += phi
        self.time += dt
        self.iteration += 1
        divergence_after = self.divergence_linf()
        return {
            "iteration": float(self.iteration),
            "time": self.time,
            "dt": dt,
            "divergence_pre_projection_linf": float(np.max(np.abs(divergence_before))),
            "divergence_post_projection_linf": divergence_after,
            "pressure_iterations": float(pressure_iterations),
            "pressure_residual_linf": pressure_residual,
            "pressure_converged": pressure_residual <= self.cfg.pressure_tolerance,
            "du_linf": float(np.max(np.abs(self.u-old_u))),
            "dv_linf": float(np.max(np.abs(self.v-old_v))),
            "kinetic_energy": self.kinetic_energy(),
        }

    def run(self, steps: int, log_every: int = 100) -> list[dict[str, float | bool]]:
        """Run a fixed number of steps; return snapshots including the final step."""
        if steps < 1 or log_every < 1:
            raise ValueError("steps and log_every must be >= 1")
        history: list[dict[str, float | bool]] = []
        for _ in range(steps):
            row = self.step()
            if self.iteration == 1 or self.iteration % log_every == 0 or _ == steps - 1:
                history.append(row.copy())
        return history

    def run_until_steady(
        self, max_steps: int, tolerance: float = 1e-7,
        min_steps: int = 100, consecutive: int = 5,
    ) -> list[dict[str, float | bool]]:
        """Run until velocity updates stay small for several consecutive steps.

        This is a numerical stopping heuristic, not proof of physical steady state.
        """
        if max_steps < 1 or min_steps < 1 or consecutive < 1 or tolerance <= 0:
            raise ValueError("max_steps, min_steps, consecutive, and tolerance must be positive")
        history: list[dict[str, float | bool]] = []
        stable_count = 0
        for _ in range(max_steps):
            row = self.step()
            history.append(row.copy())
            update = max(float(row["du_linf"]), float(row["dv_linf"]))
            if self.iteration >= min_steps and update <= tolerance:
                stable_count += 1
            else:
                stable_count = 0
            if stable_count >= consecutive:
                break
        return history

    def cell_centered_fields(self) -> tuple[FloatArray, FloatArray, FloatArray]:
        """Return cell-centred u, v, and a copy of pressure."""
        uc = 0.5 * (self.u[:, :-1] + self.u[:, 1:])
        vc = 0.5 * (self.v[:-1, :] + self.v[1:, :])
        return uc, vc, self.p.copy()
