# J2 RAAN Drift – Numerical Orbit Propagation (Cowell Formulation)

This project demonstrates high-accuracy numerical orbit propagation of an inclined Low Earth Orbit (LEO) including Earth's J2 perturbation.  
The secular Right Ascension of the Ascending Node (RAAN) drift obtained from numerical integration is validated against the analytical J2 theory.

The implementation is fully Cartesian and uses Cowell’s formulation of the equations of motion.

---

## Physical Model

The spacecraft state vector is propagated in Cartesian coordinates:

r̈ = a_kepler + a_J2

### Force Model
- Central gravity (two-body term)
- First zonal harmonic J2

### Numerical Integration
- Cowell formulation
- DOP853 (explicit Runge–Kutta 8(5,3))
- Relative tolerance: 1e-10  
- Absolute tolerance: 1e-12  

---

## Simulation Setup

Initial orbit:
- Altitude: 700 km  
- Inclination: 50°  
- Eccentricity: 0 (circular LEO)  

Propagation duration:
- 5 days  

The orbit is initialized in the equatorial plane and rotated about the x-axis to achieve the desired inclination.

---

## Analytical J2 RAAN Drift

The secular RAAN drift due to J2 is approximated by:

Ω̇ = -(3/2) J2 (R_E / p)^2 n cos(i)

where  

p = a(1 − e²)  
n = √(μ / a³)

---

## Numerical vs. Analytical Validation

```
RAAN drift comparison (J2)
Initial elements: a=7078137.0 m, e=0.000000, i=50.000 deg
Numerical slope:   -9.021868e-07 rad/s  (-4.4661 deg/day)
Theoretical slope: -8.986261e-07 rad/s  (-4.4485 deg/day)
Relative error:    0.396 %
```

The numerical secular drift matches the analytical prediction within **0.4%**, confirming the correct implementation of the J2 perturbation model and the integration scheme.

---

## Additional Consistency Checks

The simulation verifies several physical properties:

- Conservation of total mechanical energy including J2 potential  
- Conservation of h_z due to axisymmetry of the J2 gravity field  
- Expected non-conservation of |h| under J2 perturbation  
- Secular evolution of the argument of perigee  
- No secular drift in semi-major axis  

These checks confirm both physical correctness and numerical stability.

---

## Project Structure

```
src/
 ├── constants.py       # Physical constants
 ├── dynamics.py        # Acceleration model (Kepler + J2)
 ├── elements.py        # Orbital element extraction
 ├── propagation.py     # Cowell propagation (DOP853)
 ├── analysis.py        # Validation & diagnostics
 └── main.py            # Simulation entry point
```
## How to Run

Clone the repository and install the required dependencies:

```
pip install -r requirements.txt
```

From the repository root, execute:

```
python src/main.py
```

### Requirements

- numpy
- scipy
- matplotlib

---

## Purpose

This project serves as a compact demonstration of:

- Orbit propagation using Cowell’s method  
- Implementation of the J2 perturbation  
- Quantitative validation against analytical orbital mechanics theory  
- Numerical stability and conservation analysis  

It is intended as a technical portfolio project in the context of astrodynamics and flight dynamics engineering.
