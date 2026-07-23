# DirectSum 直接相加法（中文版）

`DirectSum` 求解器对每一对粒子直接计算 Newton 引力 —— ``\mathcal{O}(N^2)``，
对于中小规模问题（``N \lesssim 10^4``）以及作为其他求解器的参考解非常实用。

## 算法公式

对所有有序对 ``(i, j)``（``i \neq j``），softening 长度 ``h_i``：

```math
\mathbf{F}_i \mathrel{+}=
G \frac{m_i m_j \, \Delta \mathbf{r}}
     {(\lVert \Delta\mathbf{r}\rVert^2 + h_i^2)^{3/2}}
\qquad
\Delta \mathbf{r} = \mathbf{r}_i - \mathbf{r}_j .
```

AstroNbodySim.jl 在
[`src/directsum/gravity.jl`](https://github.com/JuliaAstroSim/AstroNbodySim.jl/blob/main/src/directsum/gravity.jl)
中实现了三个分支：

| 区间                | 函数形式       | 光滑性         | 备注           |
| ------------------- | -------------- | -------------- | -------------- |
| ``r \ge h``         | 纯 Newton      | 无             | 精确 ``1/r^2`` |
| ``0.5h \le r < h``  | 三次样条       | ``\mathcal{C}^3`` | Aarseth Plummer 等价核 |
| ``r < 0.5h``        | 四次样条       | ``\mathcal{C}^4`` | Aarseth CHICK   |

详细示例请参考英文文档
[DirectSum](@ref) 与 `examples/05-plummer.md`。
