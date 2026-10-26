# DirectSum 直接相加法

`DirectSum` 求解器对每个目标粒子遍历全部源粒子，直接累加重力加速度，
复杂度为 ``\mathcal{O}(N^2)``，没有多极展开截断误差，并且可以任意划分到
多线程、多进程或 GPU 上执行。

## 适用场景

- **中小规模问题**：论文中的基准算例在 ``N \lesssim 10^5`` 时使用直接求和；
  单核 CPU 上每步的实用规模约为 ``N \sim 10^4``，多线程、分布式内存或 GPU
  可以把这个范围继续向上扩展。
- **参考解**：作为 Tree、FFT、FDM 求解器的精度基准。
- **新核函数原型**：潮汐张量、后牛顿修正等可以先在 DirectSum 上验证。

## 算法

对目标粒子 ``i``（软化长度 ``h_i``）和源粒子 ``j``：

```math
\mathbf{a}_i = G \sum_{j \neq i} m_j \,
               \frac{\mathbf{r}_j - \mathbf{r}_i}{h_i^3} \,
               g\!\left(\frac{\lVert \mathbf{r}_j - \mathbf{r}_i \rVert}{h_i}\right),
```

其中 ``g(u)`` 是 Gadget-2 使用的三次样条软化核
（Springel, Yoshida & White 2001；Springel 2005）：

```math
g(u) =
\begin{cases}
\dfrac{32}{3} - 38.4\,u^2 + 32\,u^3, & 0 \le u < \tfrac{1}{2}, \\
\dfrac{64}{3} - 48\,u + 38.4\,u^2 - \dfrac{32}{3}\,u^3
 - \dfrac{1}{15\,u^3}, & \tfrac{1}{2} \le u < 1, \\
\dfrac{1}{u^3}, & u \ge 1.
\end{cases}
```

该核函数在 ``u=1/2`` 与 ``u=1`` 处为 ``\mathcal{C}^2``；当 ``u \ge 1`` 时退化为
通常的 Newton 引力 ``G m_j (\mathbf{r}_j-\mathbf{r}_i)/r^3``。

> 注意：[`compute_force_pairwise`](@ref) 使用的是 Plummer 软化
> ``(\lVert\Delta\mathbf{r}\rVert^2+h^2)^{-3/2}``，与主循环中的三次样条核不同。
> 论文正文将力软化统称为 Plummer 软化，但 `compute_force_from_particle` 中实际
> 使用的是上面的 Gadget-2 三次样条核。

## 软化长度

`ForceSofteningTable` 为 6 个元素，分别对应 gas、halo、disk、bulge、star、
black hole 六类粒子（`Collection` 枚举顺序）。`SimConfig` 默认取 `1.0e-4` 长度单位（`uAstro` 下为 kpc），
`GravityConfig` 默认取 `0.01u"kpc"`。

一个常用的估计是粒子平均间距（论文采用的做法）

```math
h \sim \left(\frac{V}{N}\right)^{1/3},
```

这正是 [`suggest_softlen`](@ref) 的返回值。可用 [`set_softlen!`](@ref) 或
[`suggest_softlen!`](@ref) 显式设置；软化长度必须严格为正，并同时用于力核与
自适应时间步判据。

## 并行方式

| 方式 | 实现 |
| ---- | ---- |
| 多线程 | 本地目标粒子按 `Threads.nthreads()` 分块，用 `Threads.@threads` 并行。 |
| 分布式内存 | 各进程先计算本地源粒子对其他进程目标的贡献；目标暂存在 `sim.buffer.gravToEval`，结果通过 `sim.buffer.gravResult` 返回，由 `combine_force` 累加。通信使用 `Distributed.jl`（TCP/IP），未实现 MPI 后端。 |
| GPU | 单个 all-pairs CUDA 核，每个 GPU 线程对应一个目标粒子并遍历全部源粒子（`GPU_Algorithm = AllPairs()`）。`Tiled` 仅有类型定义，尚未实现。 |

## 时间积分

粒子求解器共用 Gadget-2 风格的蛙跳（kick–drift–kick）格式：半步 kick、
drift、重新计算力、再半步 kick。默认 `TimeIntegrationAlgorithm = Leapfrog()`，
也可选 `Euler()`。

- `TimeStep > 0`：固定步长。
- `TimeStep == 0`（`SimConfig` 默认）：自适应步长，在
  `TimeBase = 1 << 30` 的整数时间轴上量化，原始步长为

  ```math
  \Delta t = \sqrt{\frac{2\,\eta\,h_i}{\lVert \mathbf{a}_i \rVert}},
  ```

  其中 `ErrTolTimestep = 0.025`，并受 `MinStep`/`MaxStep` 限制。

## 典型配置

```julia
ds = Simulation(particles;
                GravitySolver = DirectSum(),
                device        = CPU(),
                ForceSofteningTable = fill(0.1u"kpc", 6),
                TimeEnd       = 0.1u"Gyr",
                TimeStep      = 0.0u"Gyr")  # 0 表示自适应
run(ds)
```

英文文档见 [DirectSum gravity solver](@ref "DirectSum gravity solver")。
