"""Weighted spatial evaluation metrics for historical Rx1day.

Adapted from the calculation in script 09 of the author's LE_CEC_code package.
Inputs must be 1985-2014 climatologies in mm/day on the same grid and mask.
No observation merging, regridding, or substitute data are performed here.
"""
import numpy as np


def spatial_metrics(model, observation, weights, mask):
    """Return weighted COR, RMSE, RSD, PBIAS (%), and valid cell count.

    Weights must have the full grid shape; for a regular lat/lon grid they
    can be the broadcast cosine of latitude. Positive PBIAS denotes wet bias.
    RSD is the model spatial SD divided by the observation spatial SD.
    """
    model, observation, weights = [np.asarray(a, float) for a in (model, observation, weights)]
    mask = np.asarray(mask, bool)
    if not (model.shape == observation.shape == weights.shape == mask.shape):
        raise ValueError("All arrays must share the same grid shape.")
    valid = (mask & np.isfinite(model) & np.isfinite(observation)
             & np.isfinite(weights) & (weights > 0))
    if not np.any(valid):
        raise ValueError("No valid paired grid cells.")
    m, o, w = model[valid], observation[valid], weights[valid]
    mean_m, mean_o = np.average(m, weights=w), np.average(o, weights=w)
    am, ao = m - mean_m, o - mean_o
    var_m, var_o = np.average(am**2, weights=w), np.average(ao**2, weights=w)
    covariance = np.average(am * ao, weights=w)
    cor = covariance / np.sqrt(var_m * var_o) if var_m > 0 and var_o > 0 else np.nan
    rsd = np.sqrt(var_m / var_o) if var_m > 0 and var_o > 0 else np.nan
    return {"COR": float(cor), "RMSE": float(np.sqrt(np.average((m-o)**2, weights=w))),
            "RSD": float(rsd), "PBIAS_percent": float(100*(mean_m-mean_o)/mean_o) if mean_o != 0 else np.nan,
            "n_grid_cells": int(valid.sum())}
