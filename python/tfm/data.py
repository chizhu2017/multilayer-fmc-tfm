from pathlib import Path
import numpy as np


def read_matrix(path):
    """Whitespace/comma delimited numeric TXT; skip COMSOL %/# comments."""
    rows = []
    with Path(path).open(encoding="utf-8-sig") as stream:
        for line in stream:
            line = line.strip()
            if not line or line.startswith(("%", "#")):
                continue
            rows.append([float(x) for x in line.replace(",", " ").split()])
    try:
        matrix = np.asarray(rows, dtype=np.float64)
    except ValueError as exc:
        raise ValueError(f"{path}: inconsistent row widths") from exc
    if matrix.ndim != 2 or matrix.shape[0] < 2 or matrix.shape[1] < 2:
        raise ValueError(f"{path}: need >=2 samples and >=1 channel")
    if not np.isfinite(matrix).all():
        raise ValueError(f"{path}: all values must be finite")
    with np.errstate(over="ignore"):
        result = matrix.astype(np.float32)
    if not np.isfinite(result).all():
        raise ValueError(f"{path}: values exceed float32 range")
    return result


class Dataset:
    """Two passes over numbered files; keep only one transmit file in memory."""
    def __init__(self, directory):
        directory = Path(directory)
        files = list(directory.glob("*.txt"))
        self.n = len(files)
        if not self.n:
            raise ValueError("no TXT files found")
        self.paths = [directory / f"{i}.txt" for i in range(1, self.n + 1)]
        if set(files) != set(self.paths):
            raise ValueError("input must contain exactly 1.txt through N.txt")
        first = read_matrix(self.paths[0])
        self.nt = len(first)
        self.time = first[:, 0].copy()
        if first.shape[1] != self.n + 1:
            raise ValueError("file count must match receiver count")
        delta = np.diff(self.time)
        if np.any(delta <= 0):
            raise ValueError("time must strictly increase after float32 conversion")
        # MATLAB round(mean(1./diff(single(time)))); positive half rounds up.
        # Stable reduction; MATLAB's single-precision mean can differ by a few ULPs.
        mean = np.float32(np.mean(np.float32(1) / delta, dtype=np.float64))
        self.fs = np.float32(np.floor(float(mean) + .5))
        self.dt = np.float32(1) / self.fs
        if not np.isfinite(self.dt) or self.dt <= 0:
            raise ValueError("invalid sampling interval")
        # A uniform sample index is the algorithm's acquisition model.
        if np.max(np.abs(delta - self.dt)) > np.float32(.01) * self.dt:
            raise ValueError("time step must be uniform within 1% of dt")

    def read(self, index):
        matrix = read_matrix(self.paths[index])
        if matrix.shape != (self.nt, self.n + 1):
            raise ValueError(f"{self.paths[index]}: inconsistent dimensions")
        if np.max(np.abs(matrix[:, 0] - self.time)) > np.float32(.01) * self.dt:
            raise ValueError(f"{self.paths[index]}: inconsistent time axis")
        return matrix[:, 1:].copy()

    def peak_indices(self):
        # MATLAB max selects the first maximum of the SIGNED trace.
        return np.stack([
            (np.argmax(self.read(k), axis=0) + 1).astype(np.float32)
            for k in range(self.n)
        ])
