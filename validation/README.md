# Validation record / 验证记录

Run date / 执行日期：2026-10-26, Windows x86_64.

## Cases / 数据工况

| Case / 工况 | Probes / 阵元 | Channels / 通道 | Samples / 样本数 | Image / 网格 |
|---|---:|---:|---:|---|
| Deterministic synthetic / 确定性合成数据 | 4 | 16 | 1600 | 32×32 |
| Local original acquisition / 本地原始采集数据 | 11 | 121 | 1422 | 600×600 |

The real acquisition is the input directory named by the supplied `TMF.m`. Raw traces, images derived from those traces and absolute local paths are excluded from this repository. Only the aggregate error metrics are committed. The synthetic example can be reproduced entirely from the repository.

实际数据来自用户提供的 `TMF.m` 所指向的采集目录。原始波形、该数据的重建图及本地绝对路径均未上传，只提交汇总误差；合成示例可完全在仓库内复现。

## Measured agreement / 实测一致程度

Relative L2 error is `||actual-reference||₂ / ||reference||₂`, computed in float64 after loading outputs. Peak-normalized error is `max(abs(actual-reference))/max(abs(reference))`. Each reported image and envelope must satisfy **both** measures below `5×10⁻⁴`. Metadata checks require matching shapes/probes, Python/Rust agreement within `rtol=10⁻⁶`, and MATLAB comparison within `dt rtol=5×10⁻⁶`, `t0 rtol=5×10⁻⁶, atol=5×10⁻⁴ samples`.

相对 L2 误差按 `||待测值-参考值||₂/||参考值||₂` 计算；峰值归一化误差为 `最大绝对误差/参考值绝对峰值`。读回输出后用 float64 计算，图像与包络的两项误差均须小于 **5×10⁻⁴**。元数据还检查网格、阵元、采样间隔和校时编号；Python/Rust 的相对容差为 `10⁻⁶`，对 MATLAB 的 dt 相对容差为 `5×10⁻⁶`，t0 允许 `5×10⁻⁴` 样本的绝对差及 `5×10⁻⁶` 的相对差。

| Case / 工况 | Comparison / 对比 | Signed image relative L2 / 原始图像 | Envelope relative L2 / 包络 |
|---|---|---:|---:|
| Synthetic / 合成 | Python vs Rust | 1.10×10⁻⁹ | 8.43×10⁻⁸ |
| Synthetic / 合成 | MATLAB vs Python | 4.16×10⁻⁴ | 2.82×10⁻⁴ |
| Synthetic / 合成 | MATLAB vs Rust | 4.16×10⁻⁴ | 2.82×10⁻⁴ |
| Real / 实际 | Python vs Rust | 1.32×10⁻⁹ | 1.71×10⁻⁷ |
| Real / 实际 | MATLAB vs Python | 4.17×10⁻⁵ | 1.12×10⁻⁵ |
| Real / 实际 | MATLAB vs Rust | 4.17×10⁻⁵ | 1.12×10⁻⁵ |

Full metrics, including maximum absolute error and calibration metadata:

完整误差及校时元数据：

- [Synthetic comparison / 合成数据对比](demo-parity.json)
- [Real-acquisition comparison / 实际数据对比](real-data-parity.json)

The signed images match between Python/Rust to CSV formatting precision on these cases. The larger differences from MATLAB include single-precision mean-reduction effects on `fs`/`dt` and `t0`; Hilbert FFT arithmetic also differs. These values describe the tested cases, not an error bound for arbitrary acquisitions. The validation threshold was set before the Rust parity runs; it is an implementation acceptance criterion, not a physical accuracy claim.

这两组数据的 Python/Rust 原始图像差异处于 CSV 格式化精度范围；与 MATLAB 的较大差异涉及单精度求均值对 `fs/dt/t0` 的影响，包络还受 FFT 运算差异影响。上述值是实测结果，并非任意输入的误差上界。阈值在 Rust 一致性对比前已设定，用于验收移植实现；不能用来代表物理成像精度。

## Reproduction / 复现

Run the demo commands in the main README to create `data/demo`, `results/demo-python` and `results/demo-rust`. Python package versions used in this validation are in [python-requirements.txt](python-requirements.txt); Rust dependency versions and checksums are in [Cargo.lock](../rust/Cargo.lock).

先执行主 README 中的合成示例命令，生成输入及两种语言输出。验证环境 Python 顶层依赖版本见 [python-requirements.txt](python-requirements.txt)，Rust 依赖版本及校验值见 [Cargo.lock](../rust/Cargo.lock)。

For MATLAB R2023b with Signal Processing Toolbox:

使用带 Signal Processing Toolbox 的 MATLAB R2023b：

```bash
python tools/prepare_matlab_reference.py data/demo --config config/demo.json --output validation/local/demo_reference.m
matlab -batch "run('validation/local/demo_reference.m')"
python tools/compare.py results/demo-python results/demo-rust --matlab validation/local/demo_reference.mat
```

The helper copies the reference, overrides input/configuration parameters, suppresses visible plotting, and saves a MAT result before local functions. It does not change reconstruction equations. Generated scripts and MAT files stay under ignored `validation/local/`. To reproduce the real case locally, substitute the original input directory and `config/default.json`.

辅助脚本生成参考脚本副本，覆盖输入和配置、隐藏绘图，并在局部函数前保存 MAT，不改重建公式。生成脚本和 MAT 保存在被忽略的 `validation/local/`。实际工况的本地复现只需换成原输入目录及 `config/default.json`。

Local checks / 本地检查：

- Python 3.13.2: 10 tests passed / 10 项测试通过。
- Rust 1.99.1 (`x86_64-pc-windows-gnu`): 3 tests passed / 3 项测试通过。
- Rustfmt and Clippy (`--all-targets -- -D warnings`) passed / 格式与静态检查通过。
- Synthetic and real-acquisition MATLAB/Python/Rust comparisons passed / 两组数据三种语言对比通过。

GitHub Actions reruns the Python/Rust tests and synthetic comparison on pushes and pull requests. It does not execute MATLAB or access local measurement data.

GitHub Actions 在 push 和 pull request 时重新执行 Python/Rust 测试和合成数据对比，不运行 MATLAB，也不访问本地测量数据。
