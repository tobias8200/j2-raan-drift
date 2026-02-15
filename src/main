import numpy as np
import matplotlib.pyplot as plt

from constants import MU, R_EARTH
from propagation import propagate_orbit
from elements import orbital_elements
from analysis import (
    fit_slope,
    angle_wrap_pi,
    theoretical_raan_rate,
    total_specific_energy,
    angular_momentum,
)

# Initial conditions (inclined LEO)
altitude = 700e3
i_deg = 50.0
inc = np.deg2rad(i_deg)

r0 = np.array([R_EARTH + altitude, 0.0, 0.0])
v_circ = np.sqrt(MU / np.linalg.norm(r0))

# rotate velocity out of equatorial plane to set inclination about x-axis
v0 = np.array([0.0, v_circ * np.cos(inc), v_circ * np.sin(inc)])

y0 = np.concatenate((r0, v0))

# Propagation setup
t_final = 5.0 * 24.0 * 3600.0  # 5 days
t_eval = np.linspace(0.0, t_final, 5000)

sol = propagate_orbit(y0, (0.0, t_final), t_eval)

t = sol.t                       # seconds
r = sol.y[:3, :]                # (3, N)
v = sol.y[3:, :]                # (3, N)

# Orbital elements over time
a_list = np.empty(r.shape[1])
e_list = np.empty(r.shape[1])
i_list = np.empty(r.shape[1])
raan_list = np.empty(r.shape[1])
omega_list = np.empty(r.shape[1])

for k in range(r.shape[1]):
    elems = orbital_elements(r[:, k], v[:, k])
    a_list[k] = elems["a"]
    e_list[k] = elems["e"]
    i_list[k] = elems["i"]
    raan_list[k] = elems["RAAN"]
    omega_list[k] = elems["omega"]

# unwrap RAAN for a clean linear trend
raan_unwrap = np.unwrap(raan_list)

# (1) Numerical vs analytical RAAN drift comparison
m_num, b_num = fit_slope(t, raan_unwrap)  # rad/s
elems0 = orbital_elements(r[:, 0], v[:, 0])
m_th = theoretical_raan_rate(elems0["a"], elems0["e"], elems0["i"])  # rad/s

rad2deg = 180.0 / np.pi
sec_per_day = 86400.0

m_num_deg_day = m_num * rad2deg * sec_per_day
m_th_deg_day  = m_th  * rad2deg * sec_per_day
rel_err = (m_num - m_th) / m_th * 100.0

print("RAAN drift comparison (J2)")
print(f"Initial elements: a={elems0['a']:.1f} m, e={elems0['e']:.6f}, i={np.rad2deg(elems0['i']):.3f} deg")
print(f"Numerical slope:   {m_num:.6e} rad/s  ({m_num_deg_day:.4f} deg/day)")
print(f"Theoretical slope: {m_th:.6e} rad/s  ({m_th_deg_day:.4f} deg/day)")
print(f"Relative error:    {rel_err:.3f} %")

# Plot RAAN drift
plt.figure()
plt.plot(t/3600.0, raan_unwrap)
plt.xlabel("Time [hours]")
plt.ylabel("RAAN [rad]")
plt.title("RAAN Drift due to J2 (numerical)")
plt.show()

# Residual after removing linear fit
raan_fit = m_num * t + b_num
raan_res = angle_wrap_pi(raan_unwrap - raan_fit)

plt.figure()
plt.plot(t/3600.0, raan_res)
plt.xlabel("Time [hours]")
plt.ylabel("RAAN residual [rad]")
plt.title("RAAN residual after removing linear drift")
plt.show()

# (2) Energy + angular momentum checks
E = total_specific_energy(r, v)
E0 = E[0]
E_rel = (E - E0) / abs(E0)

h = angular_momentum(r, v)
hz = h[2, :]
hz0 = hz[0]
hz_rel = (hz - hz0) / abs(hz0)

h_norm = np.linalg.norm(h, axis=0)
h_norm_rel = (h_norm - h_norm[0]) / abs(h_norm[0])

plt.figure()
plt.plot(t/3600.0, E_rel)
plt.xlabel("Time [hours]")
plt.ylabel("Relative energy error [-]")
plt.title("Energy conservation (including J2 potential)")
plt.show()

plt.figure()
plt.plot(t/3600.0, hz_rel)
plt.xlabel("Time [hours]")
plt.ylabel("Relative error in $h_z$ [-]")
plt.title("$h_z$ conservation (axisymmetry check)")
plt.show()

plt.figure()
plt.plot(t/3600.0, h_norm_rel)
plt.xlabel("Time [hours]")
plt.ylabel("Relative change in |h| [-]")
plt.title("|h| is not conserved under J2 (expected)")
plt.show()

# show a(t), omega(t)
plt.figure()
plt.plot(t/3600.0, a_list - a_list[0])
plt.xlabel("Time [hours]")
plt.ylabel("Δa [m]")
plt.title("Semi-major axis variation")
plt.show()

plt.figure()
plt.plot(t/3600.0, np.unwrap(omega_list))
plt.xlabel("Time [hours]")
plt.ylabel("Argument of perigee ω [rad]")
plt.title("Argument of perigee evolution (includes J2 effects)")
plt.show()
