from dataclasses import dataclass
import math
import numpy as np
from scipy.signal import hilbert
from .data import Dataset
from .travel import probe_positions, propagation_cache


def interpolate(signal, one_based_index):
    """Preserve MATLAB's 1 <= floor(index) < Nt, including zero at Nt."""
    base = np.floor(one_based_index)
    valid = np.isfinite(base) & (base >= 1) & (base < len(signal))
    out = np.zeros_like(one_based_index, dtype=np.float32)
    i = base[valid].astype(np.int64) - 1
    alpha = one_based_index[valid] - base[valid]
    out[valid] = (np.float32(1) - alpha) * signal[i] + alpha * signal[i + 1]
    return out


def time_zero(peaks, probes, dt, c):
    direct = np.abs(probes[:, None] - probes[None, :]) * np.float32(1e-3) / np.float32(c.v1)
    candidates = peaks - direct / dt
    valid = candidates[np.isfinite(candidates) & (candidates > 0)]
    if not valid.size:
        raise ValueError("no valid time-zero candidates")
    return np.float32(np.median(valid)) + np.float32(c.t0_fine)


@dataclass
class Result:
    image: np.ndarray
    envelope: np.ndarray
    t0_index: np.float32
    dt: np.float32
    probes: np.ndarray


def reconstruct(directory, config, progress=None):
    ds = Dataset(directory)
    probes = probe_positions(config, ds.n)
    peaks = ds.peak_indices()
    t0 = time_zero(peaks, probes, ds.dt, config)
    if not np.isfinite(t0):
        raise ValueError("time zero must be finite")
    dist, travel2, travel3 = propagation_cache(config, probes)
    image1 = np.zeros((config.h1 - 1, config.nx), np.float32)
    image2 = np.zeros((config.h2, config.nx), np.float32)
    image3 = np.zeros((config.h3 + 1, config.nx), np.float32)
    z = np.arange(1, config.h1, dtype=np.float32)[:, None]
    mm_to_m = np.float32(1e-3)
    for k in range(ds.n):
        signals = ds.read(k)
        for j in range(ds.n):
            signal = signals[:, j].copy()
            cut = max(0, min(ds.nt, math.floor(float(peaks[k, j]) + config.cut_point)))
            signal[:cut] = 0
            index1 = (dist[j] + dist[k]) * mm_to_m / np.float32(config.v1) / ds.dt + t0
            weight = (z / dist[j]) * (z / dist[k])
            image1 += interpolate(signal, index1) * weight
            index2 = (travel2[j] + travel2[k]) * mm_to_m / ds.dt + t0
            image2 += interpolate(signal, index2)
            index3 = (travel3[j] + travel3[k]) * mm_to_m / ds.dt + t0
            image3 += interpolate(signal, index3)
        if progress:
            progress(k + 1, ds.n)
    image = np.concatenate((image1, image2, image3), axis=0)
    # MATLAB hilbert(matrix) acts down columns; scipy defaults to the last axis.
    envelope = np.abs(hilbert(image, axis=0)).astype(np.float32)
    if not np.isfinite(envelope).all() or envelope.max() <= 0:
        raise ValueError("reconstruction has no finite nonzero signal")
    return Result(image, envelope, t0, ds.dt, probes)
