"""Generate deterministic synthetic FMC, with independently enumerated ray paths."""
import argparse
import json
from pathlib import Path
import numpy as np


def one_way(x, z, probe, c):
    # Exhaustively enumerate absolute interface positions, independent of the LUT.
    q = np.arange(1, c["nx"] + 1, dtype=np.float64)
    first = np.hypot(q - probe, c["h1"]) * 1e-3 / c["v1"]
    if z < c["h1"]:
        return np.hypot(x-probe, z) * 1e-3 / c["v1"]
    if z < c["h1"] + c["h2"]:
        return np.min(first + np.hypot(x-q, z-c["h1"]) * 1e-3 / c["v2"])
    second = np.min(first[:, None] + np.hypot(q[None, :] - q[:, None], c["h2"]) * 1e-3 / c["v2"], axis=0)
    return np.min(second + np.hypot(x-q, z-c["h1"]-c["h2"]) * 1e-3 / c["v3"])


def generate(output, config_path):
    c = json.loads(Path(config_path).read_text(encoding="utf-8"))
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    if any(output.glob("*.txt")):
        raise ValueError("demo output already contains TXT files; choose an empty directory")
    probes = [c["x_loc"] + i*c["gap"] for i in range(4)]
    if probes[-1] > c["nx"]:
        raise ValueError("demo requires space for four probes")
    dt = 5e-8
    time = np.arange(1600) * dt
    delay = 1.5e-6
    targets = [(10, 7, 1.), (19, 17, .8), (16, 27, 1.2)]
    for k, tx in enumerate(probes, 1):
        traces = []
        for rx in probes:
            direct = abs(tx-rx)*1e-3/c["v1"] + delay
            signal = 20*np.exp(-((time-direct)/1e-7)**2)
            for x, z, amplitude in targets:
                arrival = one_way(x, z, tx, c) + one_way(x, z, rx, c) + delay
                tau = (time-arrival)/1.5e-7
                signal += amplitude*(1-2*tau**2)*np.exp(-tau**2)
            traces.append(signal)
        np.savetxt(output/f"{k}.txt", np.column_stack([time, *traces]), fmt="%.17e", delimiter="\t")
    (output/"targets.json").write_text(json.dumps({"targets_x_z_amplitude": targets, "delay_s": delay}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output")
    parser.add_argument("--config", default="config/demo.json")
    args = parser.parse_args()
    generate(args.output, args.config)
