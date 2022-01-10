# AstroNbodySim.jl

[![codecov](https://codecov.io/gh/JuliaAstroSim/AstroNbodySim.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaAstroSim/AstroNbodySim.jl)
[![][docs-dev-img]][docs-dev-url]

天体物理数值模拟计算项目。

请遵守`GPL 3.0`协议。

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

## 开发教程

1. 强烈建议使用相关包的最新`master`分支进行开发，首先 clone:
   - [ParallelOperations.jl](https://github.com/JuliaAstroSim/ParallelOperations.jl)
   - [PhysicalParticles.jl](https://github.com/JuliaAstroSim/PhysicalParticles.jl)
   - [AstroIO.jl](https://github.com/JuliaAstroSim/AstroIO.jl)
   - [AstroIC.jl](https://github.com/JuliaAstroSim/AstroIC.jl)
   - [PhysicalTrees.jl](https://github.com/JuliaAstroSim/PhysicalTrees.jl)
   - [PhysicalMeshes.jl](https://github.com/JuliaAstroSim/PhysicalMeshes.jl)
   - [AstroPlot.jl](https://github.com/JuliaAstroSim/AstroPlot.jl)
2. `dev --local [absolute path]` 安装上面的依赖包，比如
   ```jl
   pkg> dev --local /home/user/work/AstroNbodySim
   ```
3. 建议使用`VS Code`及其`Julia`扩展进行开发，`Revise.jl`可以热更新改动
4. 欢迎提交 issue 和 PR

## 新用户教程

0. 学习`Julia`基础（建议啃英文文档）：https://docs.julialang.org/en/v1/
1. 阅读依赖包的文档
2. 尝试运行`AstroNbodySim/examples`里的示例。首先使用`AstroNbodySim/examples/install_pkgs.jl`安装示例中用到的依赖包。
   默认输出函数为当前终端目录下的`./output`文件夹。如果用户定义的输出文件夹不存在，程序会自动创建一个以避免报错。
3. 如果想要终止正在运行的模拟，在输出目录下创建一个名为`stop`的文件即可（与`Gadget2`一样）
   ```
   echo > output/stop
   ```
4. 在`REPL`中使用问号`help?`可以快速查看函数支持的参数和关键字，比如
   ```jl
   help?> prepare
   search: prepare preprocessdata

     function prepare(simulation::Simulation)

     Do the following operations:
     1. Say hello
     2. Preprocess data
     3. Check the output directory, make a new one if not exist
     4. Remove "stop" file
     5. Set the global preferred units
     6. Set up logging, timing, profiling and analyzing log files
   ```

### 默认选项和注意事项

1. `Gadget2`的默认输出精度为单浮点`float`，在`Julia`里为`Float32`。但我们在整个项目里都使用双浮点`Float64`
2. 粒子类型独立存储在`Type`数组中，目前共6种：`Gas`, `Halo`, `Disk`, `Bulge`, `Star`, `BlackHole`（在`Gadget2`中是`bndry`类型）
3. 注意区分`AstroNbodySim.Plot`与`Plots`模块

## 支持与引用

本项目用于学术研究。欢迎 Star 来支持我们。如果你在研究、教学和其他活动中用到了我们的代码，请引用如下文章：
```tex
% arxiv
```

## 常见问答

## 代码开发约定

1. 命名约定（函数或变量）
   1. 尽量不要使用有可能令人迷惑的专业名词的缩写，比如黑洞和`Bames-Hut`算法都可以缩写为`BH`
   2. 为避免误解，最好将英文单词完整打出来，比如`Particle`和`Partial`都可以简写为`Part`（在读文档前，你能一眼就猜出来`Gadget2`里的`NumPart`意思是什么吗？）
   3. 给那些会修改传入参数的值的函数加上`!`后缀
   4. 变量和模块采用单词首字母大写的Pascal命名法，比如`UpperCase`；函数命名全小写加下划线，比如`lower_case`
2. 文档和注释
   1. 常用函数要有`Julia`风格的注释，包含简介、示例、参数说明
   2. 不要啰嗦
   3. 代码块和复杂算法都要有说明
3. Pull requests
   1. 合并前一定要进行过完整测试
   2. 提供更新内容简要说明
   3. 提供示例代码
4. 测试
5. 变量作用域
   1. 没有必要设置全局变量
   2. -
6. 不要上传测试用的参数文件和数据。请在`.gitignore`里标记它们
7. 适当使用`@inbounds`以提高循环效率
8. `Numerics`模块仅包含最基础的功能函数，且都有返回值或修改其参数值。高级函数都写在其他模块里
9. 每次赋值变量的时候都要小心`=`与`copy()`的不同：**`Julia`对于`A=B`形式的赋值不会创建新的变量**
   但这样并不意味着`a = 1.0; b = a`之后就能通过修改`a`的值来改变`b`的值
   `b = a`所做的事情是将常量`1.0`的引用拷贝给了`b`，而不是`a`的引用。因此`Pos = ParticleData.Pos`之类的语句有时并不适用