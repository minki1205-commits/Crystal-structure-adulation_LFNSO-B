"""Row-vector convention: A_super = S @ A_parent; reciprocal units Å⁻¹, no 2π."""
import numpy as np


def matrix(s):
    s = np.asarray(s, dtype=float)
    if s.shape != (3, 3) or not np.isfinite(s).all() or abs(np.linalg.det(s)) < 1e-10:
        raise ValueError('Expected nonsingular finite 3×3 transformation')
    return s


def reflection_to_super(h, s):
    return np.asarray(h) @ matrix(s).T


def zone_to_super(u, s):
    return np.asarray(u) @ np.linalg.inv(matrix(s))


def reciprocal_cartesian(h, lattice):
    return np.asarray(h) @ np.linalg.inv(matrix(lattice)).T


def zone_law(h, u, tolerance=1e-9):
    return bool(abs(np.dot(h, u)) < tolerance)


def structure_factors(frac, hkl, weights):
    """Explicit supplied scattering weights; not a dynamical ED intensity model."""
    frac, hkl, weights = map(np.asarray, (frac, hkl, weights))
    return np.exp(2j * np.pi * (hkl @ frac.T)) @ weights
