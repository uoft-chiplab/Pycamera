"""
make_polylog_table.py

Offline generator for the polylogarithm lookup table used by the degenerate
Fermi-gas ("polylog") fit in PyCamera.  Run this ONCE (on any machine with
scipy - Python 2 or Python 3) to produce 'polylog_table.mat'; the fit itself
then only interpolates the table and never evaluates a polylog directly.

The table stores, on a grid of the reduced variable

    u = log(fugacity) - (x-x0)^2 / (2 sigma^2)   [the fit's argument]

the three polylog orders needed by the app:

    Li2  = Li_2(-e^u)      (2D column density; used by the mock Fermi cloud)
    Li52 = Li_{5/2}(-e^u)  (1D ROI-summed profile; the fit model)
    Li3  = Li_3(-e^u)      (atom number and T/T_F conversions)

These are computed from the complete Fermi-Dirac integral

    F_j(u) = 1/Gamma(j+1) * Integral_0^inf t^j / (exp(t-u) + 1) dt

using the identity   Li_{j+1}(-e^u) = -F_j(u).  So
    Li_2  (-e^u) = -F_1  (u)
    Li_{5/2}(-e^u) = -F_{3/2}(u)
    Li_3  (-e^u) = -F_2  (u)

Usage:
    python make_polylog_table.py [output.mat] [u_min] [u_max] [n_points]
"""

import sys

from scipy.integrate import quad
from scipy.io import savemat
from scipy.special import gamma
from numpy import exp, linspace, zeros


def fermi_dirac_integral(j, u):
    """ Complete Fermi-Dirac integral F_j(u) (order j, argument u). """
    # The integrand decays like t^j * exp(-(t-u)) for t >> u, so an upper
    # limit a comfortable distance above u captures it to machine precision.
    upper = max(60.0, u + 60.0)
    integrand = lambda t: t ** j / (exp(t - u) + 1.0)
    val, _err = quad(integrand, 0.0, upper, limit=200)
    return val / gamma(j + 1.0)


def polylog_neg_exp(s, u):
    """ Li_s(-e^u) via  Li_{j+1}(-e^u) = -F_j(u)  with j = s-1. """
    return -fermi_dirac_integral(s - 1.0, u)


def build_table(u_min=-30.0, u_max=30.0, n_points=4000):
    """ Tabulate Li_2, Li_{5/2}, Li_3 of -e^u over a uniform u grid. """
    u = linspace(u_min, u_max, n_points)
    Li2 = zeros(n_points)
    Li52 = zeros(n_points)
    Li3 = zeros(n_points)
    for i in range(n_points):
        Li2[i] = polylog_neg_exp(2.0, u[i])
        Li52[i] = polylog_neg_exp(2.5, u[i])
        Li3[i] = polylog_neg_exp(3.0, u[i])
    return u, Li2, Li52, Li3


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "polylog_table.mat"
    u_min = float(sys.argv[2]) if len(sys.argv) > 2 else -30.0
    u_max = float(sys.argv[3]) if len(sys.argv) > 3 else 30.0
    n_points = int(sys.argv[4]) if len(sys.argv) > 4 else 4000

    print("Building polylog table: u in [%g, %g], %d points ..."
          % (u_min, u_max, n_points))
    u, Li2, Li52, Li3 = build_table(u_min, u_max, n_points)
    savemat(out, {"u": u, "Li2": Li2, "Li52": Li52, "Li3": Li3})
    print("Wrote %s" % out)
    # Quick sanity spot-check at u = 0:  Li_s(-1) closed forms.
    #   Li_2(-1) = -pi^2/12,  Li_3(-1) = -3/4 * zeta(3)
    print("spot-check at u=0:  Li2=%.6f (exact -0.822467)  "
          "Li3=%.6f (exact -0.901543)"
          % (polylog_neg_exp(2.0, 0.0), polylog_neg_exp(3.0, 0.0)))


if __name__ == "__main__":
    main()
