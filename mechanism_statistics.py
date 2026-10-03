"""Mechanism statistics corresponding to manuscript Figures 3 and 4.

Adapted from scripts 05 and 06 in the author's LE_CEC_code package.
This module computes statistics only; it does not reproduce final artwork.
"""
import numpy as np
from scipy.stats import pearsonr


def paired_correlation(x, y, axis=0):
    """Calculate Pearson r, two-sided p, and paired valid sample counts.

    Use one equally weighted SMILE ensemble mean per sample. For spatial
    correlations, inputs can have shape (model, latitude, longitude).
    Constant series or fewer than three finite pairs yield NaN statistics.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    if x.shape != y.shape:
        raise ValueError("Paired inputs must have identical shapes.")
    x, y = np.moveaxis(x, axis, 0), np.moveaxis(y, axis, 0)
    shape = x.shape[1:]
    xf, yf = x.reshape(x.shape[0], -1), y.reshape(y.shape[0], -1)
    r, p = np.full(xf.shape[1], np.nan), np.full(xf.shape[1], np.nan)
    n = np.zeros(xf.shape[1], dtype=int)
    for i in range(xf.shape[1]):
        finite = np.isfinite(xf[:, i]) & np.isfinite(yf[:, i])
        a, b = xf[finite, i], yf[finite, i]
        n[i] = len(a)
        if n[i] >= 3 and np.ptp(a) > 0 and np.ptp(b) > 0:
            r[i], p[i] = pearsonr(a, b)
    return {"r": r.reshape(shape), "p": p.reshape(shape),
            "n": n.reshape(shape)}


def area_mean(field, latitude, mask):
    """Cosine-latitude mean over the last two axes on a regular lat/lon grid.

    Supply the common study-region mask explicitly. For paired comparisons,
    callers must use the same finite support for both variables.
    """
    field, latitude, mask = np.asarray(field, float), np.asarray(latitude, float), np.asarray(mask, bool)
    if field.ndim < 2 or mask.shape != field.shape[-2:]:
        raise ValueError("The mask must match the final latitude/longitude axes.")
    if latitude.ndim != 1 or latitude.size != field.shape[-2]:
        raise ValueError("Latitude must match the penultimate axis.")
    weights = np.cos(np.deg2rad(latitude))[:, None]
    valid = np.isfinite(field) & mask
    denominator = np.sum(np.where(valid, weights, 0), axis=(-2, -1))
    numerator = np.sum(np.where(valid, field * weights, 0), axis=(-2, -1))
    return np.divide(numerator, denominator,
                     out=np.full_like(numerator, np.nan), where=denominator > 0)
