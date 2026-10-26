# Algorithm details / 算法细节

## 1. Calibration / 校时

Let `p[k,j]` be the **1-based** index of the first maximum of the signed trace. Coordinates are in mm, speeds in m/s, and sample interval `dt` in seconds:

```text
fs = round(mean(1 / diff(float32(time))))
dt = float32(1 / fs)
q[k,j] = p[k,j] - (abs(x[j]-x[k]) * 1e-3 / v1) / dt
t0 = median(q[finite(q) and q>0]) + t0_fine
```

`p` 是有符号回波第一个最大值的 MATLAB 1 起始编号。只在中位数计算前过滤无效和非正候选，`t0_fine` 的单位是样本，不是秒。移植版采用 float64 对倒数步长求和，求均值后转回 float32；其他校时运算保持 float32。原 MATLAB 的单精度归约可能使 `fs` 略有变化。

The ports use a float64 sum for the reciprocal-step mean, cast the mean back to float32, and retain float32 for the other calibration operations. This is a deliberate reduction-stability choice; it can yield a slightly different `fs` than MATLAB's single-precision reduction.

## 2. Discrete rays / 离散射线

The following `L` tables use **mm/(m/s)**, not seconds. Multiply by `1e-3` when converting to physical travel time.

以下表的单位是 **mm/(m/s)**。必须乘 `1e-3` 才是秒，避免与已经换算成米的第一层路径混淆。

For a horizontal integer offset `d`, define `u=0…d`:

```text
T12[d] = min_u sqrt(u²+h1²)/v1 + sqrt((d-u)²+h2²)/v2
L2[dz,d] = min_u sqrt(u²+h1²)/v1 + sqrt((d-u)²+dz²)/v2
L3[dz,d] = min_u T12[u] + sqrt((d-u)²+dz²)/v3
```

`L2` covers `dz=0…h2-1`; `L3` covers `dz=0…h3`. `T12` must be available through offset `d` before computing `L3[:,d]`. The implementation computes offsets in increasing order, preserving this dependency.

`L2` 的相对深度为 `0…h2-1`，`L3` 为 `0…h3`。第三层使用前两层的最优前缀，属于动态规划；按 d 递增计算，确保所需的 `T12[0…d]` 均已就绪。

For the flat isotropic layers in this model, an optimal path can be taken with monotone horizontal coordinates between source and target. Clamping any crossing point to that horizontal interval cannot increase its adjacent segment lengths; sorting crossing coordinates removes backtracking. Thus integer-offset enumeration can reproduce the absolute interface search within the image domain. Python tests compare the tables against unrestricted absolute crossing-point enumeration on small grids.

在此水平、各向同性层状模型里，可以选取水平坐标从源点到目标点单调的最优路径：将界面交点限制在两端点区间不会增加相邻线段长度，去掉水平回折也不会增加长度。因此可按整数水平偏移搜索，再通过平移恢复绝对坐标。测试使用小网格的绝对界面交点穷举作独立对照。

## 3. Sampling and weighting / 采样与加权

```text
Layer 1:
  r_a = sqrt((x-x[a])²+z²)
  r_b = sqrt((x-x[b])²+z²)
  it = ((r_a+r_b)*1e-3/v1)/dt+t0
  weight = (z/r_a)*(z/r_b)

Layers 2,3:
  it = (L[dz,abs(x-x[a])]+L[dz,abs(x-x[b])])*1e-3/dt+t0
  weight = 1

i = floor(it), alpha = it-i
S = (1-alpha)*signal[i-1] + alpha*signal[i]   # 0-based storage
valid only if 1 <= i < Nt; otherwise S=0
image[z,x] += S*weight
```

Muting uses `cut=floor(float64(p)+cut_point)`, clips it to `[0,Nt]`, then sets the first `cut` samples to zero. All channel sums use the original loop order (file index outer, receiver index inner). The image is not divided by the number of channels.

截断采用 `floor(float64(p)+cut_point)`，限制到 `[0,Nt]` 后将前 cut 个样本置零。通道叠加保持原脚本的文件外循环、通道内循环顺序，不除以通道数。

## 4. Envelope / 包络

```text
image = concatenate(layer1, layer2, layer3)   # shape [nz,nx]
envelope = abs(analytic_signal(image, axis=depth))
```

The FFT Hilbert multiplier is 1 at DC, 2 at strictly positive frequencies, 1 at the Nyquist bin for even lengths, and 0 at negative frequencies. The inverse RustFFT result is divided by `nz`. No extra padding is used. Depth partitions are joined **before** the transform, as in MATLAB; transforming each layer separately would alter interface behavior.

FFT 乘子为 DC=1、正频率=2、偶数长度 Nyquist=1、负频率=0，RustFFT 的逆变换还需除以 nz。没有额外补零。必须先合并三层再变换；逐层分别求包络会改变界面附近的数值。

## 5. Cost and layout / 复杂度与布局

Let `N` be probe count, `Nt` samples per trace, `P=nx*nz` image pixels, `D` maximum horizontal offset and `H=h2+h3+1`.

| Component / 环节 | Cost / 规模 |
|---|---|
| Two input passes / 两遍读取 | O(N² Nt) work, O(N Nt) working signal memory |
| Discrete tables / 走时表 | O(H D²) work, O(H D) memory |
| Probe caches / 阵元缓存 | O(N P) work and memory |
| Coherent sum / 相干叠加 | O(N² P) work |
| Hilbert / 包络 | O(nx nz log(nz)) work |

Propagation caches occupy about `4*N*nx*nz` bytes, excluding temporary arrays. The original 11-probe, 600×600 case uses approximately **15.11 MiB** for these caches. TXT parsing also temporarily holds text/numeric conversion buffers. Peak process memory is therefore larger than the cache size.

传播缓存约占 `4*N*nx*nz` 字节，不含临时数组。11 阵元、600×600 网格的缓存约 **15.11 MiB**。TXT 解析还会持有文本和数值转换缓冲，因此缓存大小不等于进程峰值内存。

Python caches use `[probe,depth,x]`; Rust stores the same order in flat row-major vectors. Both output `[depth,x]` CSV matrices. No full `Nt×N×N` FMC cube or interface-search matrix is retained.

Python 缓存使用 `[probe,depth,x]`；Rust 按相同顺序展平，x 连续。两者 CSV 都是深度行、横向列，不保留完整 FMC 三维数组或大型界面搜索矩阵。
