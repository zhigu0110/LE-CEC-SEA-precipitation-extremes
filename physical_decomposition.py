"""Historical-baseline decomposition corresponding to manuscript Equation 2.

Adapted from script 04 in the author's LE_CEC_code package.
Inputs must already be period means on the same grid and spatial support,
with one ensemble mean per SMILE. Humidity must be validated independently.
"""
import numpy as np


def decompose(rx_his, rx_future, q_his, q_future, tas_his, tas_future):
    """Return thermodynamic, dynamic, residual, and supporting fields.

    Historical period: 1985-2014. Future period: 2070-2099.
    Rx1day units: mm/day. Q units: kg/kg. Temperature units: K or Celsius,
    consistently across both periods. All input arrays must have equal shape.
    The final axis layout is chosen by the caller, for example model/lat/lon.
    Sensitivity and component units are mm/day/K.
    """
    arrays = [np.asarray(a, dtype=float) for a in
              (rx_his, rx_future, q_his, q_future, tas_his, tas_future)]
    if len({a.shape for a in arrays}) != 1:
        raise ValueError("Inputs must have identical shapes and aligned coordinates.")
    rh, rf, qh, qf, th, tf = arrays
    valid = (np.isfinite(rh) & np.isfinite(rf) & np.isfinite(qh)
             & np.isfinite(qf) & np.isfinite(th) & np.isfinite(tf)
             & (qh > 0) & (qf > 0) & (tf != th))
    dt = np.where(valid, tf - th, np.nan)
    qh_safe = np.where(valid, qh, np.nan)
    qf_safe = np.where(valid, qf, np.nan)
    beta_his = rh / qh_safe
    beta_future = rf / qf_safe
    delta_beta = beta_future - beta_his
    sens_q = (qf - qh) / dt
    return {
        "thermo_term": beta_his * sens_q,
        "dynamic_term": qh_safe * delta_beta / dt,
        "residual_term": delta_beta * sens_q,
        "total_sensitivity": (rf - rh) / dt,
        "beta_bar": beta_his,
        "sens_q": sens_q,
        "rx1day_his": np.where(valid, rh, np.nan),
        "delta_t": dt,
    }
