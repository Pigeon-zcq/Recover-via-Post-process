# Gegenbauer 重构方法 — Python 实现

本文件夹是对 Chen & Shu (2014, 2015) 两篇论文中 Gegenbauer 后处理重构算法的高精度 Python 复现，所有源 Fortran 程序均已翻译为 Python。

## 文件总览

```
Gegenbauer_Python/
├── gegenbauer_reconstruction.py          # 核心算法库（全部数值方法）
├── __init__.py                           # 包入口，导出核心库的所有符号
│
├── chen_shu_2015_fourier_examples.ipynb  # 复现 2015 论文 Examples 6.1 / 6.2
├── chen_shu_2015_transport_example.ipynb # 复现 2015 论文 变系数输运方程
├── chen_shu_2014_collocation_examples.ipynb  # 复现 2014 论文 Examples 5.1 / 5.2
├── logarithmic_gegenbauer_reconstruction_test.ipynb  # 对数奇异性的可行性实验
│
├── VariableCoeffsLinearTransport.f90     # 原始 Fortran：输运方程求解器
├── VariableCoeffsLinearTransportOptimize.f90  # 原始 Fortran：参数扫描
├── doublescollofindbestparameters.f      # 原始 Fortran：配点法最优参数搜索
├── COMMON1                               # Fortran 公共块声明文件
│
├── reproductions/                        # 自动生成的输出（CSV、Markdown 表、PNG 图）
│   ├── chen_shu_2015/                    # 快速精度复现结果
│   └── chen_shu_2015_full/               # 高精度完整复现结果
│
├── Chen和Shu_2015_*.md                   # 2015 论文阅读笔记
├── Chen和Shu_2014_*.md                   # 2014 论文阅读笔记
└── README_Python.md                      # Fortran → Python 函数名对照表
```

---

## 核心库：`gegenbauer_reconstruction.py`

整个项目的数值计算核心，使用 `mpmath` 实现任意精度算术（默认 80 位十进制精度）。
所有模块均从此文件导入。

### 功能分组

#### 1. 精度与基础工具
| 函数 | 对应 Fortran | 说明 |
|------|-------------|------|
| `set_precision(dps)` | — | 设置 `mpmath` 精度（十进制位数） |
| `as_mpf(value)` | — | 将数值安全转换为 `mp.mpf`，避免浮点截断 |
| `endpoint_power(base, exp)` | — | 端点处的幂函数，处理零点特殊情况 |

#### 2. 测试函数
| 函数 | 对应 Fortran | 说明 |
|------|-------------|------|
| `fs_periodic_singularity(z, p, q)` | `fs` | 周期化端点奇异函数 `(1+x)^(p/q)`，周期为 2 |
| `fv_transport_exact(z, p, q)` | `fv` | 变系数输运算例的精确解 |
| `paper_example_function(p, q)` | — | 返回 `cos(x) + sin(x)·(1+x)^(p/q)`（2015 论文用函数） |

#### 3. 奇异性消除映射
| 函数 | 对应 Fortran | 说明 |
|------|-------------|------|
| `transform_y_to_x(y, bet, a, b)` | `transf`（flag=1）| y 空间 → x 空间（奇异性消除映射） |
| `transform_x_to_y(x, bet, a, b)` | `transf`（flag=2）| x 空间 → y 空间（逆映射） |

映射公式：$x(y) = \left(\frac{(b-a)^\beta}{2}(1+y)\right)^{1/\beta} + a$，其中 $\beta = 1/q$。

#### 4. Gegenbauer 多项式
| 函数 | 对应 Fortran | 说明 |
|------|-------------|------|
| `rising_over_factorial(lam, max_order)` | — | 预计算系数 $(\lambda)_j/j!$ |
| `gegenbauer_grid_values(order, nn, lam)` | `gegen` | 在 Chebyshev 网格上批量计算 Gegenbauer 值 |
| `gegenbauer_norm(order, lam)` | — | 计算归一化常数 $h_l^\lambda$ |

#### 5. Fourier 系数计算
| 函数 | 对应 Fortran | 说明 |
|------|-------------|------|
| `composite_gauss_fourier_coefficients(func, n)` | `getPointValuefromFC` | 复合三点 Gauss 积分，计算前 `2N+1` 个 Fourier 系数 |
| `quad_fourier_coefficients(func, n)` | — | 高精度自适应积分（`mp.quad`），更慢但精度更高 |
| `fourier_coefficients(func, n, method)` | — | 统一接口，支持 `"gauss3"` 和 `"quad"` 两种方法 |
| `evaluate_fourier_series(coeffs, x, n)` | — | 在 x 处求值 Fourier 部分和 |

#### 6. Gegenbauer 重构流水线
| 函数 | 对应 Fortran | 说明 |
|------|-------------|------|
| `gegenbauer_coefficients_from_fourier(coeffs, n, m, lam, q, a, b)` | `gcoe` | 从 Fourier 系数计算 Gegenbauer 系数（Chebyshev-Gauss 积分） |
| `reconstruct_on_y_grid(hg, lam, q, nzn, a, b)` | — | 在 Chebyshev-y 网格上批量重构 |

#### 7. 论文算例接口
| 函数 | 说明 |
|------|------|
| `find_best_parameters(func, n, q, ...)` | 扫描 `(lambda, m)` 参数空间，寻找最小误差对，对应 Fortran `Optimize` 程序 |

#### 8. 配点法（2014 论文）
| 函数 | 对应 Fortran | 说明 |
|------|-------------|------|
| `collocation_points(nfc, npt, n, ...)` | — | 生成 Fourier 或 Chebyshev 配点网格 |
| `collocation_spectral_coefficients(values, x_nodes, nfc, n, npt)` | `ccoe` | 从配点函数值计算谱系数（支持 Fourier 和 Chebyshev 基） |
| `old_collocation_gegenbauer_coefficients(hf, hhf, nfc, n, m, lam, q)` | `gcoe`（配点版）| 2014 论文的配点法 Gegenbauer 系数计算 |

#### 9. 变系数输运方程（2015 论文 Section 6）
| 函数 | 对应 Fortran | 说明 |
|------|-------------|------|
| `fourier_galerkin_transport_coefficients(n, ...)` | `FourierGalerkinTransport` | RK4 时间推进求解 $u_t - xu_x = 0$，含指数滤波器 |
| `evaluate_transport_fourier_solution(n, ...)` | — | 评估后处理前的 Fourier-Galerkin 解 |
| `reconstruct_variable_transport(n, m, lam, ...)` | — | 对输运方程解施加 Gegenbauer 后处理 |
| `find_best_transport_parameters(n, ...)` | — | 扫描输运算例的最优后处理参数 |

---

## Notebooks

### `chen_shu_2015_fourier_examples.ipynb`
**目标**：复现 Chen & Shu (2015) 论文 Tables 1–2 和 Fig. 1（Examples 6.1, 6.2）。

- 使用 `gegenbauer_reconstruction.py` 中的 `fourier_coefficients`、`gegenbauer_coefficients_from_fourier`、`reconstruct_on_y_grid` 等函数
- 同时计算论文线性选取参数和文章给出的最优对 `(lambda, m)` 的误差
- 支持 `SEARCH_OPTIMAL = True` 模式扫描最优参数
- 以 `PIL.Image` 绘制逐点误差对数图（不依赖 matplotlib）

### `chen_shu_2015_transport_example.ipynb`
**目标**：复现 Chen & Shu (2015) 论文 Table 7（变系数输运方程后处理）。

- 内嵌 `TransportFGDouble` 类：**双精度 NumPy** 实现，与 Fortran 浮点路径完全一致，比 mpmath 快 100 倍以上
- 同时对比论文线性选取参数与扫描最优参数的误差
- 提供可选的 `RUN_MPMATH_CROSS_CHECK` 交叉验证开关（调用 `gegenbauer_reconstruction.py`）
- 绘制后处理前/后的误差对比图

### `chen_shu_2014_collocation_examples.ipynb`
**目标**：复现 Chen & Shu (2014) 论文 Examples 5.1 和 5.2（配点值 Gegenbauer 重构）。

- 使用 `collocation_points`、`collocation_spectral_coefficients`、`old_collocation_gegenbauer_coefficients` 等函数
- 同时测试 **Fourier 均匀配点**（`nfc=1`）和 **Chebyshev-Gauss 配点**（`nfc=2`）两种方案
- 实现配点值的加权处理：$f^\lambda(x_i) = (1-y_i^2)^{\lambda-1/2}(1+y_i)^{-(q-1)/2} f(x_i)$

### `logarithmic_gegenbauer_reconstruction_test.ipynb`
**目标**：探索对数奇异性 $u(x) = x\log x$ 的 Gegenbauer 重构可行性（实验性）。

- 使用对数映射 $x(y) = 2\exp\!\left(-\frac{2\log 2}{1+y}\right)$ 代替幂次映射
- 自定义 `fourier_coefficients_unit_interval`（在 $[0,1]$ 上的 Fourier 系数）
- 调用 `gegenbauer_reconstruction.py` 中的 `gegenbauer_grid_values`、`gegenbauer_norm` 等基础函数
- 扫描 `(lambda, m)` 候选网格，选出最优参数，绘制逐点误差图

---

## Fortran 原始文件

| 文件 | 功能 |
|------|------|
| `VariableCoeffsLinearTransport.f90` | 主求解程序：计算 Fourier 系数、RK4 时间推进、Gegenbauer 后处理 |
| `VariableCoeffsLinearTransportOptimize.f90` | 参数扫描版本：对所有 `(m, lambda)` 组合求最小误差 |
| `doublescollofindbestparameters.f` | 配点法参数搜索（双精度 Fortran 77） |
| `COMMON1` | Fortran 公共块声明，被上述程序 `INCLUDE` |

Python 与 Fortran 的函数名对应关系详见 `README_Python.md`。

---

## 依赖关系图

```
Fortran 源文件
    │  (翻译)
    ▼
gegenbauer_reconstruction.py   ← 所有数值算法
    │
    ├── __init__.py             ← 包导出
    │
    ├── chen_shu_2015_fourier_examples.ipynb
    │       (fourier_coefficients, gegenbauer_coefficients_from_fourier,
    │        reconstruct_on_y_grid, find_best_parameters)
    │
    ├── chen_shu_2015_transport_example.ipynb
    │       (内嵌 TransportFGDouble 双精度实现 + 可选 mpmath 交叉验证)
    │
    ├── chen_shu_2014_collocation_examples.ipynb
    │       (collocation_points, collocation_spectral_coefficients,
    │        old_collocation_gegenbauer_coefficients, reconstruct_on_y_grid)
    │
    └── logarithmic_gegenbauer_reconstruction_test.ipynb
            (gegenbauer_grid_values, gegenbauer_norm, rising_over_factorial
             + 自定义对数映射和 Fourier 系数函数)
```

---

## 快速上手

```bash
# 安装依赖
pip install mpmath numpy matplotlib pillow

# 运行配点法算例（2014 论文）
jupyter notebook Gegenbauer_Python/chen_shu_2014_collocation_examples.ipynb

# 运行 Fourier 算例（2015 论文 Examples 6.1 / 6.2）
jupyter notebook Gegenbauer_Python/chen_shu_2015_fourier_examples.ipynb

# 运行输运方程算例（2015 论文 Table 7）
jupyter notebook Gegenbauer_Python/chen_shu_2015_transport_example.ipynb

# 在 Python 中直接调用核心函数
python3 -c "
from Gegenbauer_Python import (
    fourier_coefficients,
    gegenbauer_coefficients_from_fourier,
    paper_example_function,
    reconstruct_on_y_grid,
    set_precision,
)
set_precision(60)
N, q, lam, m = 40, 2, 2, 1
func = paper_example_function(1, q)
coeffs = fourier_coefficients(func, N)
hg = gegenbauer_coefficients_from_fourier(coeffs, N, m, lam, q)
xs, _, values = reconstruct_on_y_grid(hg, lam, q, N)
errors = [abs(u - func(x)) for x, u in zip(xs, values)]
print('max error =', max(errors))
"
```

---

## 参考文献

- Chen, Z. & Shu, C.-W. (2014). *Recovering Exponential Accuracy from Collocation Point Values of Smooth Functions with End-Point Singularities*. Journal of Computational and Applied Mathematics.
- Chen, Z. & Shu, C.-W. (2015). *Recovering Exponential Accuracy in Fourier Spectral Methods Involving Piecewise Smooth Functions with Unbounded Derivative Singularities*. Journal of Scientific Computing.
