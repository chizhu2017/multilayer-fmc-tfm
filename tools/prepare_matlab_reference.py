"""Create a runnable reference copy; modify only parameters, plotting and output."""
import argparse
import json
import re
from pathlib import Path


def prepare(source, data, config, output):
    source = Path(source).read_text(encoding="utf-8-sig")
    c = json.loads(Path(config).read_text(encoding="utf-8"))
    source = source.replace("clc;", "").replace("clear;", "").replace("close all;", "")
    data = str(Path(data).resolve()).replace("\\", "/").replace("'", "''")
    source = re.sub(r"file_path\s*=\s*\.\.\.\s*'[^']*';", lambda _: f"file_path = '{data}';", source, count=1)
    for key, matlab in (("nx", "Nx"), ("h1", "h1"), ("h2", "h2"), ("h3", "h3"),
                        ("v1", "v1"), ("v2", "v2"), ("v3", "v3"), ("x_loc", "x_loc"),
                        ("gap", "gap"), ("cut_point", "cut_point"), ("t0_fine", "t0_fine")):
        value = f"single({c[key]})" if key in ("v1", "v2", "v3", "t0_fine") else str(c[key])
        source = re.sub(rf"^{matlab}\s*=.*?;", f"{matlab} = {value};", source, count=1, flags=re.M)
    source = source.replace("save_result = true;", "save_result = false;")
    source = source.replace("figure('Color'", "figure('Visible','off','Color'")
    output = Path(output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    mat = str(output.with_suffix(".mat")).replace("\\", "/").replace("'", "''")
    marker = "%% ================================================================\n% 局部函数："
    source = source.replace(marker, f"save('{mat}', 'TFM_img','TFM_env','t0_index','dt','probeX','-v7');\nclose all;\n\n" + marker)
    if f"save('{mat}'" not in source:
        raise ValueError("reference insertion point not found")
    output.write_text(source, encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data")
    parser.add_argument("--source", default="matlab/TMF_original.m")
    parser.add_argument("--config", default="config/default.json")
    parser.add_argument("--output", default="validation/local/reference.m")
    args = parser.parse_args()
    prepare(args.source, args.data, args.config, args.output)
