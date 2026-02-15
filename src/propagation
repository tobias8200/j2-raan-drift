from scipy.integrate import solve_ivp
from dynamics import state_derivative

# Propagates the orbit using Cowell's formulation
# Integration is performed with the high-order explicit Runge-Kutta method DOP853.
def propagate_orbit(y0, t_span, t_eval=None):
    sol = solve_ivp(
        state_derivative,
        t_span,
        y0,
        method="DOP853",
        rtol=1e-10,
        atol=1e-12,
        t_eval=t_eval
    )

    return sol
