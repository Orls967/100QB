"""Fungsi karakteristik Heston (formulasi Little Heston Trap, Albrecher et al. 2007)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class HestonParams:
    """Parameter model Heston: ``(v0, kappa, theta, xi, rho)``."""

    v0: float
    kappa: float
    theta: float
    xi: float
    rho: float

    def as_tuple(self) -> tuple[float, float, float, float, float]:
        """Kembalikan tuple 5 parameter."""
        return (self.v0, self.kappa, self.theta, self.xi, self.rho)

    @property
    def feller_lhs(self) -> float:
        """Nilai ``2 * kappa * theta`` pada kondisi Feller."""
        return 2.0 * self.kappa * self.theta

    @property
    def feller_rhs(self) -> float:
        """Nilai ``xi^2`` pada kondisi Feller."""
        return self.xi * self.xi

    @property
    def feller_satisfied(self) -> bool:
        """Apakah kondisi Feller ``2 * kappa * theta > xi^2`` terpenuhi."""
        return self.feller_lhs > self.feller_rhs


DEFAULT_PARAMS = HestonParams(v0=0.04, kappa=2.0, theta=0.04, xi=0.5, rho=-0.7)


def heston_cf(
    u: np.ndarray | complex | float,
    tau: float,
    s0: float,
    r: float,
    params: HestonParams,
) -> np.ndarray | complex:
    """Fungsi karakteristik ``E[exp(i * u * ln S_T)]`` dengan Little Heston Trap."""
    u_arr = np.asarray(u, dtype=np.complex128)
    v0, kappa, theta, xi, rho = params.as_tuple()

    x = np.log(s0) + r * tau
    beta = kappa - rho * xi * 1j * u_arr
    d = np.sqrt(beta * beta + (xi * xi) * (u_arr * u_arr + 1j * u_arr))
    g = (beta - d) / (beta + d)

    exp_dt = np.exp(-d * tau)
    c_term = (kappa * theta / (xi * xi)) * (
        (beta - d) * tau - 2.0 * np.log((1.0 - g * exp_dt) / (1.0 - g))
    )
    d_term = ((beta - d) / (xi * xi)) * ((1.0 - exp_dt) / (1.0 - g * exp_dt))
    phi = np.exp(1j * u_arr * x + c_term + d_term * v0)
    return complex(phi) if np.ndim(u) == 0 else phi


def heston_cumulants(tau: float, s0: float, r: float, params: HestonParams) -> tuple[float, float, float]:
    """Kumulan ke-1, ke-2, dan ke-4 dari ``ln(S_T)`` secara numerik dari ``ln phi(u)``."""
    h = 1e-3
    ks = np.array([-2, -1, 0, 1, 2], dtype=float)
    u_pts = -1j * ks * h  # MGF M(s) = phi(-i*s)
    psi = np.real(np.log(heston_cf(u_pts, tau, s0, r, params)))
    c1 = float((-psi[4] + 8.0 * psi[3] - 8.0 * psi[1] + psi[0]) / (12.0 * h))
    c2 = float((-psi[4] + 16.0 * psi[3] - 30.0 * psi[2] + 16.0 * psi[1] - psi[0]) / (12.0 * h * h))
    c4 = float((psi[4] - 4.0 * psi[3] + 6.0 * psi[2] - 4.0 * psi[1] + psi[0]) / (h**4))
    return c1, max(c2, 1e-12), max(c4, 0.0)
