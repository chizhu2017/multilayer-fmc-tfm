"""Compare common CSV outputs, optionally against an original MATLAB MAT file."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.io import loadmat


def metric(reference, actual):
    reference = np.asarray(reference, dtype=np.float64)
    actual = np.asarray(actual, dtype=np.float64)
    if reference.shape != actual.shape:
        raise ValueError(f"shape mismatch: {reference.shape} vs {actual.shape}")
    if not np.isfinite(actual).all() or not np.isfinite(reference).all():
        raise ValueError("nonfinite output")
    delta = actual - reference
    scale = max(float(np.max(np.abs(reference))), np.finfo(float).tiny)
    norm = max(float(np.linalg.norm(reference)), np.finfo(float).tiny)
    return {"max_abs_error": float(np.max(np.abs(delta))),
            "max_error_over_peak": float(np.max(np.abs(delta)) / scale),
            "relative_l2_error": float(np.linalg.norm(delta) / norm)}


def compare(left, right, matlab=None, tolerance=5e-4):
    report = {"tolerance_relative_l2_and_peak": tolerance, "python_vs_rust": {}}
    py_meta = json.loads((Path(left)/"metadata.json").read_text())
    rs_meta = json.loads((Path(right)/"metadata.json").read_text())
    report["metadata"] = {
        "python": {k: py_meta[k] for k in ("shape", "dt", "t0_index", "probe_x_mm")},
        "rust": {k: rs_meta[k] for k in ("shape", "dt", "t0_index", "probe_x_mm")},
    }
    for key in ("shape", "probe_x_mm", "dt", "t0_index"):
        if not np.allclose(py_meta[key], rs_meta[key], rtol=1e-6, atol=0):
            raise ValueError(f"metadata mismatch: {key}")
    reference = loadmat(matlab) if matlab else None
    if reference is not None:
        report["metadata"]["matlab"] = {
            "dt": float(reference["dt"].item()), "t0_index": float(reference["t0_index"].item())}
    passed = True
    for name, field in (("image", "TFM_img"), ("envelope", "TFM_env")):
        py = np.loadtxt(Path(left)/f"{name}.csv", delimiter=",", ndmin=2)
        rs = np.loadtxt(Path(right)/f"{name}.csv", delimiter=",", ndmin=2)
        metrics = metric(py, rs)
        report["python_vs_rust"][name] = metrics
        passed &= metrics["relative_l2_error"] < tolerance and metrics["max_error_over_peak"] < tolerance
        if reference is not None:
            for label, array in (("matlab_vs_python", py), ("matlab_vs_rust", rs)):
                metrics = metric(reference[field], array)
                report.setdefault(label, {})[name] = metrics
                passed &= metrics["relative_l2_error"] < tolerance and metrics["max_error_over_peak"] < tolerance
    if reference is not None:
        for key in ("dt", "t0_index"):
            atol = 5e-4 if key == "t0_index" else 0
            if not np.allclose(reference[key].item(), py_meta[key], rtol=5e-6, atol=atol):
                raise ValueError(f"MATLAB metadata mismatch: {key}")
    report["passed"] = bool(passed)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("python_output")
    parser.add_argument("rust_output")
    parser.add_argument("--matlab")
    parser.add_argument("--report", type=Path)
    parser.add_argument("--tolerance", type=float, default=5e-4)
    args = parser.parse_args()
    report = compare(args.python_output, args.rust_output, args.matlab, args.tolerance)
    encoded = json.dumps(report, indent=2)
    print(encoded)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(encoded + "\n", encoding="utf-8")
    raise SystemExit(0 if report["passed"] else 1)
