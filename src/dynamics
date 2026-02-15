import numpy as np
from constants import MU, R_EARTH, J2


def kepler_acceleration(r):
    r_norm = np.linalg.norm(r)
    return -MU * r / r_norm**3


def j2_acceleration(r):
    x, y, z = r
    r_norm = np.linalg.norm(r)

    factor = (3/2) * J2 * MU * R_EARTH**2 / r_norm**5
    zx_ratio = (z**2) / r_norm**2

    ax = factor * x * (5*zx_ratio - 1)
    ay = factor * y * (5*zx_ratio - 1)
    az = factor * z * (5*zx_ratio - 3)

    return np.array([ax, ay, az])


def total_acceleration(r):
    return kepler_acceleration(r) + j2_acceleration(r)


def state_derivative(t, y):
    r = y[:3]
    v = y[3:]

    a = total_acceleration(r)

    return np.concatenate((v, a))
