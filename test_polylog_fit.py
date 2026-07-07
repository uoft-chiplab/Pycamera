"""
test_polylog_fit.py

Hardware-free numeric self-test for the degenerate-Fermi-gas (polylog) fit.

The production fit lives in data_analysis2.py / experiment2.py, which are
Python-2.5 / Enthought-Traits / wxPython modules that cannot be imported on a
modern dev machine.  The *numerics* of the fit (the polylog model + a leastsq
fit + the lookup table), however, are plain scipy and run fine under Python 3.
This script therefore mirrors 'polylog1D' / 'polylog_fit1D' exactly and checks:

  1. the committed 'polylog_table.mat' against known closed-form values, and
  2. that a leastsq fit recovers the log-fugacity 'q' (and hence T/T_F) from
     synthetic Fermi profiles built with the same model.

Run:  python test_polylog_fit.py [polylog_table.mat]
Exits non-zero on any failure.
"""

import sys

from numpy import (arange, exp, sqrt, pi, argmax, array, mean, zeros_like,
                   random)
from scipy.io import loadmat
from scipy.interpolate import interp1d
from scipy.optimize import leastsq
from scipy.special import zeta


# --- Mirror of data_analysis2.polylog1D / polylog_fit1D ---------------------

def polylog1D(x, params, F):
    A = params[0]
    sig2 = abs(params[1]) + 1e-12
    x0 = params[2]
    bg = params[3]
    q = params[4]
    u = q - (x - x0) ** 2 / (2 * sig2)
    return A * F(u) / F(q) + bg


def polylog_fit1D(xdata, gdata, p0, F):
    errfunc = lambda p, x, val: polylog1D(x, p, F) - val
    fit, _ = leastsq(errfunc, array(p0, dtype=float), args=(xdata, gdata))
    return fit


def load_table(path):
    tbl = loadmat(path)
    u = tbl['u'].ravel()
    F = interp1d(u, tbl['Li52'].ravel(), bounds_error=False, fill_value=0.0)
    G = interp1d(u, tbl['Li3'].ravel(), bounds_error=False, fill_value=0.0)
    return u, F, G


def ttf_from_q(q, G):
    """ T/T_F = ( -6 Li_3(-e^q) )^(-1/3). """
    return (-6.0 * float(G(q))) ** (-1.0 / 3.0)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "polylog_table.mat"
    u, F, G = load_table(path)
    failures = []

    # --- 1. Table spot-checks at u=0 (Li_s(-1) closed forms) ----------------
    # Li_3(-1)   = -3/4 zeta(3)
    # Li_{5/2}(-1) = -eta(5/2) = -(1 - 2^{-3/2}) zeta(5/2)
    li3_exact = -0.75 * zeta(3)
    li52_exact = -(1.0 - 2.0 ** (-1.5)) * zeta(2.5)
    for name, got, exact in (("Li3(-1)", float(G(0.0)), li3_exact),
                             ("Li52(-1)", float(F(0.0)), li52_exact)):
        err = abs(got - exact)
        ok = err < 1e-3
        print("table %-9s: got %.6f  exact %.6f  |err|=%.2e  %s"
              % (name, got, exact, err, "OK" if ok else "FAIL"))
        if not ok:
            failures.append(name)

    # --- 2. Fit recovery of q from synthetic Fermi profiles -----------------
    random.seed(1)
    A_true, sig2_true, x0_true, bg_true = 1000.0, 100.0, 250.0, 50.0
    xdata = arange(150, 350).astype(float)
    print("\nfit recovery (seed q=1.0 always):")
    for q_true in (0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0):
        clean = polylog1D(xdata, [A_true, sig2_true, x0_true, bg_true, q_true], F)
        noise = 0.003 * A_true * random.standard_normal(clean.shape)
        gdata = clean + noise
        # Seed the way the app does: gaussian-style guesses + q0 = 1.0.
        seed = [gdata.max() - gdata.min(), 120.0,
                xdata[argmax(gdata)], gdata.min(), 1.0]
        fit = polylog_fit1D(xdata, gdata, seed, F)
        q_fit = fit[4]
        resid = polylog1D(xdata, fit, F) - gdata
        rmse_rel = sqrt(mean(resid ** 2)) / A_true
        ttf_t = ttf_from_q(q_true, G)
        ttf_f = ttf_from_q(q_fit, G)
        ttf_relerr = abs(ttf_f - ttf_t) / ttf_t
        # q is well constrained once the gas is at least mildly degenerate.
        q_ok = (q_true < 1.0) or (abs(q_fit - q_true) < 0.5)
        shape_ok = rmse_rel < 0.01
        ttf_ok = (q_true < 1.0) or (ttf_relerr < 0.15)
        ok = q_ok and shape_ok and ttf_ok
        print("  q_true=%4.1f -> q_fit=%6.3f  T/TF: true=%.4f fit=%.4f "
              "(relerr %.1f%%)  rmse/A=%.4f  %s"
              % (q_true, q_fit, ttf_t, ttf_f, 100 * ttf_relerr, rmse_rel,
                 "OK" if ok else "FAIL"))
        if not ok:
            failures.append("q_true=%.1f" % q_true)

    print()
    if failures:
        print("FAILURES: %s" % ", ".join(failures))
        sys.exit(1)
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
