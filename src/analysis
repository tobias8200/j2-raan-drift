import numpy as np
from constants import MU, R_EARTH, J2


# Specific Keplerian energy (two-body reference)
def specific_energy(r, v):
    r_norm = np.linalg.norm(r, axis=0)
    v_norm = np.linalg.norm(v, axis=0)
    return v_norm**2 / 2 - MU / r_norm


# Linear least-squares slope
def fit_slope(t, y):
    m, b = np.polyfit(t, y, 1)
    return m, b


# Wrap angles to [-pi, pi]
def angle_wrap_pi(x):
    return (x + np.pi) % (2*np.pi) - np.pi


# Analytical J2 secular RAAN drift rate (rad/s)
def theoretical_raan_rate(a, e, inc):
    p = a * (1.0 - e**2)
    n = np.sqrt(MU / a**3)
    return -(3.0/2.0) * J2 * (R_EARTH**2 / p**2) * n * np.cos(inc)


# Total gravitational potential including J2 correction
def j2_potential(r):
    r = np.asarray(r)

    if r.ndim == 1:
        r = r.reshape(3, 1)

    x, y, z = r
    r_norm = np.linalg.norm(r, axis=0)
    r2 = r_norm**2
    z2 = z**2

    factor = J2 * (R_EARTH**2 / r2) * 0.5 * (3.0 * z2 / r2 - 1.0)
    U = -MU / r_norm * (1.0 - factor)

    return U[0] if U.size == 1 else U


# Specific total mechanical energy (including J2 potential)
def total_specific_energy(r, v):
    v_norm = np.linalg.norm(v, axis=0)
    return 0.5 * v_norm**2 + j2_potential(r)


# Specific angular momentum vector
def angular_momentum(r, v):
    return np.cross(r.T, v.T).T
