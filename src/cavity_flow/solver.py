"""Readable MAC-grid lid-driven cavity solver starter.

The discretization is deliberately kept explicit for study. It needs independent
numerical verification before benchmark-accuracy claims should be made.
"""
from __future__ import annotations
from dataclasses import dataclass
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
        if self.reynolds <= 0 or self.length <= 0:
            raise ValueError("reynolds and length must be positive")
        if not 0 < self.cfl <= 1:
            raise ValueError("cfl must be in (0, 1]")
        if self.pressure_tolerance <= 0 or self.pressure_max_iterations < 1:
            raise ValueError("pressure tolerance and iteration limit must be positive")


class CavitySolver:
    """2D incompressible cavity with staggered velocity and cell pressure arrays."""

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
        # Horizontal velocity: stationary side walls; stationary bottom; moving top lid.
        u[:, 0] = u[:, -1] = 0.0
        u[0, :] = 0.0
        u[-1, :] = self.cfg.lid_speed
        # Vertical velocity is zero at all impermeable walls.
        v[0, :] = v[-1, :] = 0.0
        v[:, 0] = v[:, -1] = 0.0

    def divergence(self, u: FloatArray | None = None, v: FloatArray | None = None) -> FloatArray:
        """Finite-volume divergence at cell centres."""
        uu = self.u if u is None else u
        vv = self.v if v is None else v
        return (uu[:, 1:] - uu[:, :-1]) / self.dx + (vv[1:, :] - vv[:-1, :]) / self.dy

    def divergence_linf(self) -> float:
        return float(np.max(np.abs(self.divergence())))

    def kinetic_energy(self) -> float:
        uc = 0.5 * (self.u[:, :-1] + self.u[:, 1:])
        vc = 0.5 * (self.v[:-1, :] + self.v[1:, :])
        return float(0.5 * np.mean(uc**2 + vc**2))

    def stable_dt(self) -> float:
        umax = max(float(np.max(np.abs(self.u))), abs(self.cfg.lid_speed), 1e-12)
        vmax = max(float(np.max(np.abs(self.v))), 1e-12)
        dt_adv = self.cfg.cfl / (umax / self.dx + vmax / self.dy)
        dt_diff = 0.24 / (self.nu * (1.0 / self.dx**2 + 1.0 / self.dy**2))
        estimate = min(dt_adv, dt_diff)
        return min(self.cfg.dt, estimate) if self.cfg.dt is not None else estimate

    def _tendencies(self) -> tuple[FloatArray, FloatArray]:
        """Centred advection/diffusion tendencies, interpolating cross velocities."""
        u, v, dx, dy, nu = self.u, self.v, self.dx, self.dy, self.nu
        du_dx = np.zeros_like(u); du_dy = np.zeros_like(u); lap_u = np.zeros_like(u)
        du_dx[:, 1:-1] = (u[:, 2:] - u[:, :-2]) / (2 * dx)
        du_dy[1:-1, :] = (u[2:, :] - u[:-2, :]) / (2 * dy)
        lap_u[:, 1:-1] += (u[:, 2:] - 2*u[:, 1:-1] + u[:, :-2]) / dx**2
        lap_u[1:-1, :] += (u[2:, :] - 2*u[1:-1, :] + u[:-2, :]) / dy**2
        vcc = 0.5 * (v[:-1, :] + v[1:, :])
        v_at_u = np.zeros_like(u)
        v_at_u[:, 1:-1] = 0.5 * (vcc[:, :-1] + vcc[:, 1:])
        fu = -(u * du_dx + v_at_u * du_dy) + nu * lap_u

        dv_dx = np.zeros_like(v); dv_dy = np.zeros_like(v); lap_v = np.zeros_like(v)
        dv_dx[:, 1:-1] = (v[:, 2:] - v[:, :-2]) / (2 * dx)
        dv_dy[1:-1, :] = (v[2:, :] - v[:-2, :]) / (2 * dy)
        lap_v[:, 1:-1] += (v[:, 2:] - 2*v[:, 1:-1] + v[:, :-2]) / dx**2
        lap_v[1:-1, :] += (v[2:, :] - 2*v[1:-1, :] + v[:-2, :]) / dy**2
        ucc = 0.5 * (u[:, :-1] + u[:, 1:])
        u_at_v = np.zeros_like(v)
        u_at_v[1:-1, :] = 0.5 * (ucc[:-1, :] + ucc[1:, :])
        fv = -(u_at_v * dv_dx + v * dv_dy) + nu * lap_v
        return fu, fv

    def _poisson(self, rhs: FloatArray) -> tuple[FloatArray, int, float]:
        """Solve the cell-centred pressure correction with homogeneous Neumann walls.

        The Jacobi diagonal counts only real neighboring cells. The residual uses
        mirrored ghost values, which represent zero-normal-gradient boundary faces.
        This matches the Laplacian induced by divergence(gradient(phi)) on this MAC grid.
        """
        nx, ny = self.cfg.nx, self.cfg.ny
        dx2, dy2 = self.dx**2, self.dy**2
        rhs = np.asarray(rhs, dtype=float)
        if rhs.shape != (ny, nx):
            raise ValueError(f"rhs must have shape {(ny, nx)}, got {rhs.shape}")
        # A Neumann Poisson problem is solvable only when the RHS has zero mean.
        rhs = rhs - np.mean(rhs)
        phi = np.zeros_like(rhs)

        diagonal = np.full_like(rhs, 2.0 / dx2 + 2.0 / dy2)
        diagonal[:, 0] -= 1.0 / dx2
        diagonal[:, -1] -= 1.0 / dx2
        diagonal[0, :] -= 1.0 / dy2
        diagonal[-1, :] -= 1.0 / dy2

        def apply_laplacian(a: FloatArray) -> FloatArray:
            east = np.empty_like(a); west = np.empty_like(a)
            north = np.empty_like(a); south = np.empty_like(a)
            east[:, :-1], east[:, -1] = a[:, 1:], a[:, -1]
            west[:, 1:], west[:, 0] = a[:, :-1], a[:, 0]
            north[:-1, :], north[-1, :] = a[1:, :], a[-1, :]
            south[1:, :], south[0, :] = a[:-1, :], a[0, :]
            return (east - 2.0*a + west) / dx2 + (north - 2.0*a + south) / dy2

        residual = float("inf")
        for it in range(1, self.cfg.pressure_max_iterations + 1):
            # Sum only actual neighbors; boundary ghost values are not unknowns.
            weighted = np.zeros_like(phi)
            weighted[:, :-1] += phi[:, 1:] / dx2
            weighted[:, 1:] += phi[:, :-1] / dx2
            weighted[:-1, :] += phi[1:, :] / dy2
            weighted[1:, :] += phi[:-1, :] / dy2
            phi = (weighted - rhs) / diagonal
            phi -= np.mean(phi)  # fix the arbitrary pressure gauge
            if it == 1 or it % 10 == 0:
                residual = float(np.max(np.abs(apply_laplacian(phi) - rhs)))
                if residual <= self.cfg.pressure_tolerance:
                    return phi, it, residual
        residual = float(np.max(np.abs(apply_laplacian(phi) - rhs)))
        return phi, self.cfg.pressure_max_iterations, residual

    def step(self) -> dict[str, float]:
        dt = self.stable_dt()
        old_u, old_v = self.u.copy(), self.v.copy()
        fu, fv = self._tendencies()
        us, vs = self.u + dt*fu, self.v + dt*fv
        self.apply_boundaries_to(us, vs)
        rhs = self.divergence(us, vs) / dt
        phi, pit, pres = self._poisson(rhs)
        un, vn = us.copy(), vs.copy()
        un[:, 1:-1] -= dt * (phi[:, 1:] - phi[:, :-1]) / self.dx
        vn[1:-1, :] -= dt * (phi[1:, :] - phi[:-1, :]) / self.dy
        self.u, self.v = un, vn
        self.apply_velocity_boundaries()
        self.p += phi
        self.time += dt
        self.iteration += 1
        return {
            "iteration": float(self.iteration), "time": self.time, "dt": dt,
            "divergence_linf": self.divergence_linf(), "pressure_iterations": float(pit),
            "pressure_residual_linf": pres,
            "du_linf": float(np.max(np.abs(self.u-old_u))),
            "dv_linf": float(np.max(np.abs(self.v-old_v))),
            "kinetic_energy": self.kinetic_energy(),
        }

    def run(self, steps: int, log_every: int = 100) -> list[dict[str, float]]:
        if steps < 1 or log_every < 1:
            raise ValueError("steps and log_every must be >= 1")
        history = []
        for _ in range(steps):
            row = self.step()
            if self.iteration == 1 or self.iteration % log_every == 0:
                history.append(row.copy())
        return history

    def cell_centered_fields(self) -> tuple[FloatArray, FloatArray, FloatArray]:
        uc = 0.5 * (self.u[:, :-1] + self.u[:, 1:])
        vc = 0.5 * (self.v[:-1, :] + self.v[1:, :])
        return uc, vc, self.p.copy()
