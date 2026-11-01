"""
Ionization and recombination rate coefficients of deuterium from OpenADAS (through CHERAB),
loaded the same way as in atomic/coronal_eq_solver.py.

Units: n_e in m^-3, T_e in eV, rates in m^3/s.
Needs the CHERAB module (e.g. `module load CHERAB/1.5.0-gfbf-2023b`).
"""
import numpy as np
from cherab.core.atomic import deuterium
from cherab.openadas import OpenADAS

_adas = OpenADAS(permit_extrapolation=False)
_rate_ion   = _adas.ionisation_rate(deuterium, 0)      # D0 -> D+
_rate_recom = _adas.recombination_rate(deuterium, 1)   # D+ -> D0

# Relative step of the finite differences used for the derivatives
_REL_STEP = 1.0e-3

# Relative margin kept from the edges of the ADAS tables
_EDGE_MARGIN = 1.0e-4


def _inner_range(value_range):
    # the table limits are not always evaluable exactly at the edge (rounding), so stay slightly inside
    return value_range[0] * (1.0 + _EDGE_MARGIN), value_range[1] * (1.0 - _EDGE_MARGIN)


def _clamp(value, value_range):
    return min(max(value, value_range[0]), value_range[1])


def _rate_and_derivatives(rate, n_e, t_e, zero_below_range=False):
    """
    Rate coefficient and its derivatives with respect to n_e and T_e.

    n_e and T_e are clamped to the range of the ADAS table (the rates are held constant
    outside it, so the derivatives are zero there). With zero_below_range, the rate and its
    derivatives are zero instead when T_e is below the lower limit of the table
    (above the upper limit, and for any n_e, the rate is still held constant).
    """
    n_lo, n_hi = _inner_range(rate.density_range)
    t_lo, t_hi = _inner_range(rate.temperature_range)

    if zero_below_range and t_e < t_lo:
        return 0.0, 0.0, 0.0

    n_c = _clamp(n_e, (n_lo, n_hi))
    t_c = _clamp(t_e, (t_lo, t_hi))
    value = rate(n_c, t_c)

    d_n = 0.0
    if n_lo < n_e < n_hi:
        n_p = min(n_c * (1.0 + _REL_STEP), n_hi)
        n_m = max(n_c * (1.0 - _REL_STEP), n_lo)
        d_n = (rate(n_p, t_c) - rate(n_m, t_c)) / (n_p - n_m)

    d_t = 0.0
    if t_lo < t_e < t_hi:
        t_p = min(t_c * (1.0 + _REL_STEP), t_hi)
        t_m = max(t_c * (1.0 - _REL_STEP), t_lo)
        d_t = (rate(n_c, t_p) - rate(n_c, t_m)) / (t_p - t_m)

    return value, d_n, d_t


def ionisation_rate(n_e, t_e):
    """S(n_e, T_e) [m^3/s] and dS/dn_e, dS/dT_e. Zero below the lower T_e limit of the ADAS table."""
    return _rate_and_derivatives(_rate_ion, n_e, t_e, zero_below_range=True)


def recombination_rate(n_e, t_e):
    """alpha(n_e, T_e) [m^3/s] and dalpha/dn_e, dalpha/dT_e. Held constant outside the ADAS table."""
    return _rate_and_derivatives(_rate_recom, n_e, t_e)
