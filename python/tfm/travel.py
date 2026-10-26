import numpy as np


def probe_positions(config, n):
    positions = config.x_loc + np.arange(n, dtype=np.int64) * config.gap
    if np.any(positions < 1) or np.any(positions > config.nx):
        raise ValueError("probe coordinates must lie in 1..nx mm")
    return positions.astype(np.float32)


def travel_tables(c, max_dx):
    """One-way times in mm/(m/s); conversion to seconds occurs at sampling."""
    v1, v2, v3 = map(np.float32, (c.v1, c.v2, c.v3))
    dz2 = np.arange(c.h2, dtype=np.float32)[:, None]
    dz3 = np.arange(c.h3 + 1, dtype=np.float32)[:, None]
    t12 = np.empty(max_dx + 1, np.float32)
    lut2 = np.empty((c.h2, max_dx + 1), np.float32)
    lut3 = np.empty((c.h3 + 1, max_dx + 1), np.float32)
    for d in range(max_dx + 1):
        u = np.arange(d + 1, dtype=np.float32)
        rem = np.float32(d) - u
        first = np.sqrt(u**2 + np.float32(c.h1)**2) / v1
        t12[d] = np.min(first + np.sqrt(rem**2 + np.float32(c.h2)**2) / v2)
        lut2[:, d] = np.min(first + np.sqrt(rem**2 + dz2**2) / v2, axis=1)
        lut3[:, d] = np.min(t12[:d + 1] + np.sqrt(rem**2 + dz3**2) / v3, axis=1)
    return lut2, lut3


def propagation_cache(c, probes):
    x = np.arange(1, c.nx + 1, dtype=np.float32)
    z = np.arange(1, c.h1, dtype=np.float32)[:, None]
    dx = np.abs(x[None, :] - probes[:, None]).astype(np.int64)
    lut2, lut3 = travel_tables(c, int(dx.max()))
    dist1 = np.sqrt((x[None, None, :] - probes[:, None, None])**2 + z**2)
    return dist1, lut2[:, dx].transpose(1, 0, 2).copy(), lut3[:, dx].transpose(1, 0, 2).copy()
