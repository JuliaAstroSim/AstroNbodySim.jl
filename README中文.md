# AstroNbodySim.jl

[![codecov](https://codecov.io/gh/JuliaAstroSim/AstroNbodySim.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaAstroSim/AstroNbodySim.jl)
[![][docs-dev-img]][docs-dev-url]

天体物理数值模拟计算项目。

请遵守`GPL 3.0`协议。

## 安装

```julia
]add AstroNbodySim
```
or
```julia
]add https://github.com/JuliaAstroSim/AstroNbodySim.jl
```

**GPU加速需要提前安装 [NVIDIA CUDA toolkit](https://developer.nvidia.com/cuda-toolkit)**

## 文档

- [**Dev**][docs-dev-url] &mdash; *开发版本文档.*

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://juliaastrosim.github.io/AstroNbodySim.jl/dev

## 特性

- 带单位计算
- 用户友好的
  - 较为完善的文档
  - 易读的代码结构
  - 广播的数组操作
  - 支持 `Float16`, `Float32`, `Float64`, `Int128`, `BigFloat`, `Measurement` 等多种类型
- 跨平台，易于部署：`Linux`, `Windows`, `MacOS`
- 多种形式的并行：多线程，分布式并行，GPU加速
- 模块化和可扩展性：10+ packages，通用化设计
- 实时可视化
- 高覆盖的自动测试

## 特性预览

### GPU 模拟加速及实时可视化
![Realtime visualization on GPU](https://github.com/JuliaAstroSim/AstroNbodySim.jl/blob/main/docs/src/examples/pics/readme/Plummer.gif)

### 星系碰撞
![Galactic collision](https://github.com/JuliaAstroSim/AstroNbodySim.jl/blob/main/docs/src/examples/pics/readme/GalacticCollision.gif)

### 误差传递
![Uncertainty propagation](https://github.com/JuliaAstroSim/AstroNbodySim.jl/blob/main/docs/src/examples/pics/examples/01-binary/Uncertainty%20of%20elliptic%20orbit.png)

### 背景引力势场的自动微分
![Autodiff of background potential field](https://github.com/JuliaAstroSim/AstroNbodySim.jl/blob/main/docs/src/examples/pics/examples/AutodiffBackground.png)

### 用户自定义模拟流程：潮汐瓦解（TDE）

位置可视化：

![User-difined pipeline: Tidal disruption event (TDE)](https://github.com/JuliaAstroSim/AstroNbodySim.jl/blob/main/docs/src/examples/pics/examples/07-TDEcluster/TDE-elliptic-mosaic.png)

吸积历史：

![User-difined pipeline: TDE accretion history](https://github.com/JuliaAstroSim/AstroNbodySim.jl/blob/main/docs/src/examples/pics/examples/07-TDEcluster/TDE-elliptic-AccretionHistory.png)

### 拉格朗日半径和标度半径（scale radius）
![Lagrange radii and scale radius](https://github.com/JuliaAstroSim/AstroNbodySim.jl/blob/main/docs/src/examples/pics/examples/03-plummer/Plummer-LagrangianRadii.png)

![Lagrange radii and scale radius](https://github.com/JuliaAstroSim/AstroNbodySim.jl/blob/main/docs/src/examples/pics/examples/03-plummer/Plummer-ScaleRadius.png)

### 太阳系可视化
![Solar System](https://github.com/JuliaAstroSim/AstroNbodySim.jl/blob/main/docs/src/examples/pics/readme/SolarSystem.gif)

### 轨道回溯
![Trajectory lookback](https://github.com/JuliaAstroSim/AstroNbodySim.jl/blob/main/docs/src/examples/pics/readme/traj_lookback_all_arrow_axes_20250721v2_2Gyr.gif)

## 支持与引用

本项目用于学术研究。欢迎 Star 来支持我们。如果你在研究、教学和其他活动中用到了我们的代码，请引用如下文章：
```tex
%TODO WIP
```

## 代码架构

```mermaid
graph LR
    A[AstroSimBase.jl] --> C[PhysicalParticles.jl]
    B[ParallelOperations.jl] --> C
    C --> D[AstroIO.jl]
    C --> E[AstroIC.jl]
    C --> F[PhysicalFDM.jl]
    C --> G[PhysicalFFT.jl]
    B --> H[PhysicalTrees.jl]
    C --> H
    B --> I[PhysicalMeshes.jl]
    C --> I
    D --> J[AstroPlot.jl]
    H --> J
    I --> J
    J --> K[AstroNbodySim.jl]
    D --> K
    E --> K
    F --> K
    G --> K
    H --> K
    I --> K
    J --> K
    K --> L[WaveDM.jl]
    F --> L
    G --> L
    I --> L
    N[GalacticDynamics.jl] --> E
```


- Basic data structure: [PhysicalParticles.jl](https://github.com/JuliaAstroSim/PhysicalParticles.jl)
- File I/O: [AstroIO.jl](https://github.com/JuliaAstroSim/AstroIO.jl)
- Initial Condition: [AstroIC.jl](https://github.com/JuliaAstroSim/AstroIC.jl)
- Parallelism: [ParallelOperations.jl](https://github.com/JuliaAstroSim/ParallelOperations.jl)
- Trees: [PhysicalTrees.jl](https://github.com/JuliaAstroSim/PhysicalTrees.jl)
- Meshes: [PhysicalMeshes.jl](https://github.com/JuliaAstroSim/PhysicalMeshes.jl)
- Plotting: [AstroPlot.jl](https://github.com/JuliaAstroSim/AstroPlot.jl)
- Simulation: [AstroNbodySim](https://github.com/JuliaAstroSim/AstroNbodySim.jl)
- Benchmark: [BenchmarkPlots](https://github.com/JuliaAstroSim/BenchmarkPlots.jl)
- Parameter space exploration: [ParameterSpace](https://github.com/JuliaAstroSim/ParameterSpace.jl)