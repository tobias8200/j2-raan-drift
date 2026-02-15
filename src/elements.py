import numpy as np
from constants import MU


def orbital_elements(r, v):
    r_norm = np.linalg.norm(r)
    v_norm = np.linalg.norm(v)

    h = np.cross(r, v)
    h_norm = np.linalg.norm(h)

    k = np.array([0, 0, 1])
    n = np.cross(k, h)
    n_norm = np.linalg.norm(n)

    e_vec = (1/MU) * (np.cross(v, h) - MU * r / r_norm)
    e = np.linalg.norm(e_vec)

    energy = v_norm**2 / 2 - MU / r_norm
    a = -MU / (2 * energy)

    i = np.arccos(h[2] / h_norm)

    if n_norm != 0:
        RAAN = np.arccos(n[0] / n_norm)
        if n[1] < 0:
            RAAN = 2*np.pi - RAAN
    else:
        RAAN = 0.0

    if n_norm != 0 and e > 1e-10:
        omega = np.arccos(np.dot(n, e_vec) / (n_norm * e))
        if e_vec[2] < 0:
            omega = 2*np.pi - omega
    else:
        omega = 0.0

    return {
        "a": a,
        "e": e,
        "i": i,
        "RAAN": RAAN,
        "omega": omega
    }
