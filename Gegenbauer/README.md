# Gegenbauer Reconstruction for Functions with End-Point Singularities

从 Fourier 系数或配点值重构带端点奇异性函数的指数精度近似，基于坐标变换 + Gegenbauer 多项式展开。

## 数学背景

### 问题设置

考虑形如

$$
f(x) = a(x) + b(x)(1+x)^s, \quad x \in [a, b], \quad 0 < s = \frac{p}{q} < 1
$$

的函数，其中 $a(x)$, $b(x)$ 解析，但 $f$ 在左端点 $x=a$ 处导数无界。目标是从其 Fourier 系数（或配点值）重构出 $f$ 的指数精度逼近。

### 算法流程

**Step 1 — 坐标变换**（将奇异性消除）

引入一对一变换 $y \leftrightarrow x$（$\beta = 1/q$）：

$$
x = \Bigl[\tfrac{1}{2}(b-a)^\beta (1+y)\Bigr]^{1/\beta} + a, \qquad
y = -1 + \frac{(x-a)^\beta}{\tfrac{1}{2}(b-a)^\beta}
$$

变换后 $F(y) = f(x(y))$ 在 $[-1,1]$ 上解析。

**Step 2 — 近似 Gegenbauer 系数**（Chen & Shu 2015, eq. 3.1）

$$
\hat{g}^\lambda(l) = \frac{1}{h_l^\lambda} \int_{-1}^{1} (1-y^2)^{\lambda - \frac{1}{2}}\, f_N(x(y))\, C_l^\lambda(y)\, dy
$$

其中 $f_N(x) = \sum_{|k|\le N} \hat{f}_k e^{ik\pi x}$ 为截断 Fourier 级数，$h_l^\lambda$ 为归一化常数，
积分用 Chebyshev-Gauss 求积（$N_t - 1$ 个内点）近似。

Gegenbauer 多项式通过升阶乘 $(\lambda)_j / j!$ 的三角级数表示：

$$
C_n^\lambda(\cos\theta) = \sum_{j=0}^{n} \frac{(\lambda)_j}{j!} \cdot \frac{(\lambda)_{n-j}}{(n-j)!} \cdot \cos\bigl((n-2j)\theta\bigr)
$$

**Step 3 — 重构**

$$
f^{m,\lambda}_N(x) = \sum_{l=0}^{m} \hat{g}^\lambda(l)\, C_l^\lambda\!\bigl(y(x)\bigr)
$$

当 $\lambda = \alpha N$, $m = \gamma N$（$\alpha, \gamma$ 为适当常数）时，最大模误差指数衰减。

## 文件结构

```
Gegenbauer/
│
├── test_GegenbaurReconstruct.m   ← 主测试脚本（从这里开始）
├── GegenbaurReconstruct.m        ← 入口（调用上面的测试脚本）
│
├── 核心算法
│   ├── compute_gt.m              Pochhammer 系数 (λ)_j/j! 预计算
│   ├── gegen.m                   Chebyshev 节点处 C_n^λ(cos θ) 求值
│   ├── gcoe.m                    Fourier 系数 → Gegenbauer 系数
│   ├── gegenbauer_reconstruct.m  用 Gegenbauer 系数重构 f(x)
│   ├── transf_y2x.m              坐标变换 y → x
│   └── transf_x2y.m              坐标变换 x → y
│
├── 辅助工具
│   ├── getPointValuefromFC.m     3 点 Gauss 求积计算 Fourier 系数
│   ├── run_reconstruction.m      单次重构驱动（Fourier → Gegenbauer → 误差）
│   ├── findBestParameters.m      网格搜索最优 (λ, m) 参数对
│   └── FourierGalerkinTransport.m  RK4 求解 u_t − xu_x = 0（Example 6.5）
│
├── 测试函数
│   ├── f_example1.m              cos(x) + sin(x)·(1+x)^(p/q)
│   ├── f_transport_ic.m          输运方程初始条件 (1+x)^(p/q)
│   └── f_transport_exact.m       输运方程精确解 (1+2x)^(p/q) @ T=log2
│
├── 原始 Fortran 代码（参考）
│   ├── COMMON1                   Fortran 公共块（全局变量）
│   ├── doublescollofindbestparameters.f   配点法 + 最优参数搜索
│   ├── VariableCoeffsLinearTransport.f90  Fourier 系数法 + 输运方程
│   └── VariableCoeffsLinearTransportOptimize.f90  带参数优化版本
│
└── 参考文献（Markdown 格式）
    ├── Chen和Shu - 2014 - ...Collocation....md
    └── Chen和Shu - 2015 - ...Fourier....md
```

## 快速开始

在 MATLAB 中，进入本目录后运行：

```matlab
>> test_GegenbaurReconstruct
```

脚本将依次输出：

1. **Example 6.1** — $u(x) = \cos x + \sin x \cdot \sqrt{1+x}$，线性参数 $\lambda = N/16$, $m = 3N/80$
2. **Example 6.1** — 最优参数搜索（$N = 40, 80$）
3. **Example 6.2** — $u(x) = \cos x + \sin x \cdot (1+x)^{1/3}$，线性参数 $\lambda = N/8$, $m = N/32$
4. **Example 6.5** — 输运方程后处理，在子区间 $[-0.5, 0]$ 上恢复 $\sqrt{1+2x}$
5. **误差图** — Example 6.1，$N=80$ 的逐点误差（对数坐标）

### 自定义调用示例

```matlab
% 定义待重构函数（p=1, q=2，即 sqrt 奇异性）
p = 1;  q = 2;
f = @(x) f_example1(x, p, q);

% 参数设置
N   = 80;
lam = 5;    % Gegenbauer 参数 lambda
m   = 3;    % 保留项数
bet = 1/q;  % 变换指数

% 计算 Fourier 系数
hf = getPointValuefromFC(N, f);

% 预计算 Pochhammer 系数
gt = compute_gt(lam, max(m, 3*N) + 1);

% 计算 Gegenbauer 系数
hg = gcoe(hf, N, m, lam, bet, -1, 1, gt, 3*N);

% 在 N+1 个 Chebyshev 节点上重构
[u, xx] = gegenbauer_reconstruct(hg, m, lam, bet, -1, 1, N, gt);

% 误差
err = max(abs(u - arrayfun(f, xx)));
fprintf('L_inf error = %.3e\n', err);
```

### 搜索最优参数

```matlab
p = 1;  q = 2;
f = @(x) f_example1(x, p, q);
[m_opt, lam_opt, err_opt] = findBestParameters(40, p, q, -1, 1, f, f);
fprintf('Best: lambda=%d, m=%d, err=%.3e\n', lam_opt, m_opt, err_opt);
```

## 参数选取指导

由 Theorem 4.3（Chen & Shu 2015），当 $\lambda = \alpha N$, $m = \gamma N$ 时误差指数小。
实用参考值：

| 函数类型 | $\lambda$ | $m$ |
|---|---|---|
| $s = 1/2$（$q=2$） | $N/16$ | $3N/80$ |
| $s = 1/3$（$q=3$） | $N/8$  | $N/32$  |
| 输运方程子区间 | $N/80$ | $N/4$ |

最优参数可通过 `findBestParameters` 在 $[1, N/2]^2$ 上穷举得到。

## 数值结果预期（Table 1, Chen & Shu 2015）

$u(x) = \cos x + \sin x \cdot \sqrt{1+x}$，线性参数选取：

| $N$  | $\lambda$ | $m$ | $L^\infty$ 误差 |
|------|-----------|-----|----------------|
| 40   | 2         | 1   | ~5.9e-1        |
| 80   | 5         | 3   | ~2.9e-1        |
| 160  | 10        | 6   | ~1.6e-2        |
| 320  | 20        | 12  | ~1.3e-5        |

（误差随 $N$ 指数衰减。）

## 实现说明

- **大参数数值稳定性**：归一化常数 $h_l^\lambda$ 含 $\Gamma(l + 2\lambda)$，对大 $\lambda$ 会溢出。代码用 `gammaln` 在对数空间计算后取指数。
- **Pochhammer 系数**：`compute_gt` 用递推式 $g_j = g_{j-1}(\lambda+j-1)/j$ 避免直接调用 $\Gamma$ 函数。
- **Fourier 系数索引**：`hf(k+N+1)` 存储 $\hat{f}_k$，$k = -N, \ldots, N$。
- **Gegenbauer 系数索引**：`hg(l+1)` 存储 $\hat{g}^\lambda(l)$，$l = 0, \ldots, m$。
- **`gt` 数组索引**：`gt(j+1)` 存储 $(\lambda)_j / j!$，$j = 0, 1, \ldots$

## 参考文献

1. **Chen, Z. & Shu, C.-W. (2014)**. Recovering Exponential Accuracy from Collocation Point Values of Smooth Functions with End-Point Singularities. *J. Comput. Appl. Math.* **265**, 83–95.

2. **Chen, Z. & Shu, C.-W. (2015)**. Recovering Exponential Accuracy in Fourier Spectral Methods Involving Piecewise Smooth Functions with Unbounded Derivative Singularities. *J. Sci. Comput.* **65**, 1145–1165.
