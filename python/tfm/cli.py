import argparse
import json
from pathlib import Path
from time import perf_counter
import numpy as np
from scipy.io import savemat
from .config import Config
from .reconstruction import reconstruct


def main():
    parser = argparse.ArgumentParser(description="Three-layer FMC-TFM reconstruction")
    parser.add_argument("input", type=Path, help="directory containing 1.txt..N.txt")
    parser.add_argument("--config", type=Path)
    parser.add_argument("--output", type=Path, default=Path("results/python"))
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()
    try:
        config = Config.load(args.config) if args.config else Config()
        start = perf_counter()
        result = reconstruct(args.input, config, lambda k, n: print(f"FMC {k}/{n}", flush=True))
        args.output.mkdir(parents=True, exist_ok=True)
        np.savetxt(args.output / "image.csv", result.image, delimiter=",", fmt="%.9g")
        np.savetxt(args.output / "envelope.csv", result.envelope, delimiter=",", fmt="%.9g")
        metadata = {"config": config.to_dict(), "shape": list(result.image.shape),
                    "t0_index": float(result.t0_index), "dt": float(result.dt),
                    "probe_x_mm": result.probes.tolist(), "elapsed_seconds": perf_counter() - start,
                    "units": {"coordinates": "mm", "dt": "s", "t0_index": "samples (MATLAB 1-based convention)"}}
        (args.output / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
        np.savez_compressed(args.output / "result.npz", image=result.image, envelope=result.envelope,
                            t0_index=result.t0_index, dt=result.dt, probeX=result.probes)
        savemat(args.output / "result.mat", {"TFM_img": result.image, "TFM_env": result.envelope,
                "t0_index": result.t0_index, "dt": result.dt, "probeX": result.probes,
                **{k: getattr(config, k) for k in ("h1", "h2", "h3", "v1", "v2", "v3")}})
        if not args.no_plot:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(9, 7.5))
            plot = ax.imshow(result.envelope, origin="lower", cmap="jet", aspect="equal",
                             extent=(.5, config.nx + .5, .5, config.nz + .5),
                             vmin=0, vmax=float(result.envelope.max()) * config.color_ratio)
            if config.show_interface:
                for depth in (config.h1, config.h1 + config.h2):
                    ax.axhline(depth, color="white", linestyle="--", linewidth=1.2)
            ax.set(xlabel="X (mm)", ylabel="Depth (mm)", title="Multi-layer TFM Reconstruction")
            fig.colorbar(plot, ax=ax)
            fig.tight_layout()
            fig.savefig(args.output / "reconstruction.png", dpi=160)
            plt.close(fig)
        print(json.dumps(metadata, indent=2))
    except (ValueError, OSError) as exc:
        parser.exit(1, f"error: {exc}\n")


if __name__ == "__main__":
    main()
