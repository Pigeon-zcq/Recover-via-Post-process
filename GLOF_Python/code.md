# Fourier 系数到 GLOF 重构的算法流程

## 0. 目标与输入输出

本文档给出一个可编程实现的算法流程，用于从有限 Fourier 系数重构具有端点对数弱奇性的函数

$$
f(x)=x^r(-\log x)^k g(x),
\qquad x\in(0,1].
$$

其中：

- $r>0$；
- $k\in\mathbb N_0$；
- $g(x)$ 在 $[0,1]$ 附近解析或充分光滑；
- 奇性集中在 $x=0$；
- 已知有限 Fourier 系数 $\{\hat f_n\}_{|n|\le N_F}$。

本文采用如下 Fourier 约定：

$$
\hat f_n=\int_0^1 f(x)e^{-2\pi i n x}\,dx,
\qquad n\in\mathbb Z.
$$

Fourier 部分和定义为

$$
f_{N_F}(x)=\sum_{n=-N_F}^{N_F}\hat f_n e^{2\pi i n x}.
$$

目标是构造 GLOF 后处理重构

$$
R_{N_F,M}f(x)
=
\sum_{\ell=0}^{M}\tilde c_{\ell,N_F}
\mathcal S_\ell^{(\alpha,\beta,\lambda)}(x),
$$

其中 $\tilde c_{\ell,N_F}$ 由有限 Fourier 系数计算得到。

---

## 1. 符号约定

为避免编程实现中混淆，统一使用以下符号。

| 符号 | 含义 |
|---|---|
| $N_F$ | Fourier 截断半带宽，已知模式 $n=-N_F,\ldots,N_F$ |
| $M$ | GLOF 重构最高阶数 |
| $Q$ | Gauss-Laguerre / Gauss-GLOF 求积最高编号，共 $Q+1$ 个节点 |
| $\alpha$ | Laguerre / log 权重参数，要求 $\alpha>-1$ |
| $\beta$ | log 映射参数，要求 $\beta>-1$ |
| $\lambda$ | GLOF 权函数指数参数，$\lambda\in\mathbb R$ |
| $\nu$ | GLOF 幂因子指数，定义为 $\nu=(\beta-\lambda)/2$ |
| $r$ | 原函数端点奇性中的幂指数 |
| $k$ | 原函数端点奇性中的 log 幂次 |
| $x_j$ | 映射回 $(0,1)$ 上的求积节点 |
| $y_j$ | Laguerre-Gauss 节点 |
| $\omega_j^{(\alpha)}$ | Laguerre-Gauss 权重 |
| $\chi_j^{(\alpha,\beta,\lambda)}$ | GLOF 求积权重 |
| $\mathcal S_\ell^{(\alpha,\beta,\lambda)}$ | 第 $\ell$ 个 GLOF 基函数 |
| $\gamma_\ell^{(\alpha,\beta)}$ | GLOF 正交归一化常数 |
| $\tilde c_{\ell,N_F}$ | 由有限 Fourier 系数得到的近似 GLOF 系数 |

---

# 2. 参数选择

## 2.1 基本可积性条件

GLOF 权函数定义为

$$
\chi^{(\alpha,\lambda)}(x)
=
(-\log x)^\alpha x^\lambda.
$$

加权 $L^2$ 范数为

$$
\|u\|_{\chi^{(\alpha,\lambda)}}^2
=
\int_0^1 |u(x)|^2(-\log x)^\alpha x^\lambda\,dx.
$$

对模型函数

$$
f(x)=x^r(-\log x)^k
$$

要求

$$
\int_0^1 |f(x)|^2(-\log x)^\alpha x^\lambda\,dx<\infty.
$$

即

$$
\int_0^1 x^{2r+\lambda}(-\log x)^{2k+\alpha}\,dx<\infty.
$$

因此需要

$$
2r+\lambda>-1,
$$

以及

$$
2k+\alpha>-1.
$$

由于 $k\in\mathbb N_0$ 且通常取 $\alpha>-1$，第二个条件自然满足。

---

## 2.2 奇性匹配参数

GLOF 基函数含有幂因子

$$
x^{\frac{\beta-\lambda}{2}}.
$$

为了匹配原函数的端点幂奇性 $x^r$，令

$$
\frac{\beta-\lambda}{2}=r.
$$

因此选择

$$
\boxed{
\beta=\lambda+2r.
}
$$

此时

$$
\mathcal S_\ell^{(\alpha,\beta,\lambda)}(x)
=
x^r
L_\ell^{(\alpha)}
\left(-(\lambda+2r+1)\log x\right).
$$

这样，$x^r$ 被基函数精确吸收，而 $(-\log x)^k$ 变成 Laguerre 变量中的多项式因子。

---

## 2.3 推荐参数策略

### 2.3.1 直接 GLOF 投影型参数

若已有函数点值或可以直接积分，推荐固定参数：

$$
\lambda=0,
\qquad
\beta=2r,
\qquad
M\ \text{逐步增加}.
$$

此时

$$
\mathcal S_\ell^{(\alpha,2r,0)}(x)
=
x^r L_\ell^{(\alpha)}\left(-(2r+1)\log x\right).
$$

### 2.3.2 Fourier 后处理型参数

若从有限 Fourier 系数后处理，模仿 Chen-Shu 的参数思想，可令参数随 $N_F$ 增长：

$$
M=\lfloor \gamma N_F\rfloor,
$$

$$
\lambda=\lfloor \theta N_F\rfloor,
$$

$$
\beta=\lambda+2r.
$$

其中

$$
0<\gamma<1,
\qquad
0<\theta<1.
$$

实际计算中可先取

$$
\gamma\in[0.05,0.25],
\qquad
\theta\in[0.02,0.20].
$$

⚠️ 注意：$\lambda\sim N_F$ 有助于增强 $x=0$ 端点的权函数消失阶，但并不自动处理 $x=1$ 端点的分部积分边界项。若 Fourier 周期延拓在 $x=1$ 有跳跃，建议先做端点匹配、局部化或分段重构。

---

# 3. log 映射与 GLOF 基函数

## 3.1 log 映射

定义

$$
y=y(x)=-(\beta+1)\log x,
\qquad x\in(0,1).
$$

反变换为

$$
x=x(y)=e^{-y/(\beta+1)},
\qquad y\in(0,\infty).
$$

微分关系为

$$
\frac{dy}{dx}=-(\beta+1)x^{-1},
$$

$$
dx=-\frac{x}{\beta+1}\,dy.
$$

由于

$$
x:0^+\to 1
\quad\Longleftrightarrow\quad
y:\infty\to0,
$$

所以

$$
\int_0^1 F(x)\,dx
=
\frac1{\beta+1}
\int_0^\infty
F\left(e^{-y/(\beta+1)}\right)
e^{-y/(\beta+1)}\,dy.
$$

此外，

$$
-\log x=\frac{y}{\beta+1}.
$$

因此

$$
(-\log x)^\alpha x^\beta\,dx
=
\frac1{(\beta+1)^{\alpha+1}}y^\alpha e^{-y}\,dy.
$$

更一般地，

$$
(-\log x)^\alpha x^\lambda\,dx
=
x^{\lambda-\beta}
\frac1{(\beta+1)^{\alpha+1}}y^\alpha e^{-y}\,dy.
$$

---

## 3.2 Laguerre 多项式

令 $L_\ell^{(\alpha)}(y)$ 为广义 Laguerre 多项式，满足正交性

$$
\int_0^\infty
L_\ell^{(\alpha)}(y)
L_m^{(\alpha)}(y)
y^\alpha e^{-y}\,dy
=
\gamma_\ell^{(\alpha)}\delta_{\ell m},
$$

其中

$$
\gamma_\ell^{(\alpha)}
=
\frac{\Gamma(\ell+\alpha+1)}{\Gamma(\ell+1)}.
$$

用于数值计算的三项递推为

$$
L_0^{(\alpha)}(y)=1,
$$

$$
L_1^{(\alpha)}(y)=-y+\alpha+1,
$$

$$
L_{\ell+1}^{(\alpha)}(y)
=
\frac{2\ell+\alpha+1-y}{\ell+1}L_\ell^{(\alpha)}(y)
-
\frac{\ell+\alpha}{\ell+1}L_{\ell-1}^{(\alpha)}(y),
\qquad \ell\ge1.
$$

---

## 3.3 LOF 基函数

LOF 定义为

$$
\mathcal S_\ell^{(\alpha,\beta)}(x)
=
L_\ell^{(\alpha)}(y(x))
=
L_\ell^{(\alpha)}\left(-(\beta+1)\log x\right).
$$

其正交性为

$$
\int_0^1
\mathcal S_\ell^{(\alpha,\beta)}(x)
\mathcal S_m^{(\alpha,\beta)}(x)
(-\log x)^\alpha x^\beta\,dx
=
\gamma_\ell^{(\alpha,\beta)}\delta_{\ell m}.
$$

归一化常数为

$$
\gamma_\ell^{(\alpha,\beta)}
=
\frac{\Gamma(\ell+\alpha+1)}
{(\beta+1)^{\alpha+1}\Gamma(\ell+1)}.
$$

---

## 3.4 GLOF 基函数

GLOF 定义为

$$
\boxed{
\mathcal S_\ell^{(\alpha,\beta,\lambda)}(x)
=
x^{\frac{\beta-\lambda}{2}}
L_\ell^{(\alpha)}
\left(-(\beta+1)\log x\right).
}
$$

记

$$
\nu=\frac{\beta-\lambda}{2}.
$$

则

$$
\mathcal S_\ell^{(\alpha,\beta,\lambda)}(x)
=
x^\nu L_\ell^{(\alpha)}(y(x)).
$$

GLOF 权函数为

$$
\chi^{(\alpha,\lambda)}(x)
=
(-\log x)^\alpha x^\lambda.
$$

GLOF 正交性为

$$
\boxed{
\int_0^1
\mathcal S_\ell^{(\alpha,\beta,\lambda)}(x)
\mathcal S_m^{(\alpha,\beta,\lambda)}(x)
\chi^{(\alpha,\lambda)}(x)\,dx
=
\gamma_\ell^{(\alpha,\beta)}\delta_{\ell m}.
}
$$

其中

$$
\gamma_\ell^{(\alpha,\beta)}
=
\frac{\Gamma(\ell+\alpha+1)}
{(\beta+1)^{\alpha+1}\Gamma(\ell+1)}.
$$

---

# 4. Gauss-GLOF 求积节点与权重

## 4.1 Laguerre-Gauss 节点

令 $\{y_j^{(\alpha)}\}_{j=0}^{Q}$ 是 $L_{Q+1}^{(\alpha)}(y)$ 的 $Q+1$ 个正零点。

令 $\{\omega_j^{(\alpha)}\}_{j=0}^{Q}$ 为对应的 Laguerre-Gauss 权重，满足

$$
\int_0^\infty p(y)y^\alpha e^{-y}\,dy
=
\sum_{j=0}^{Q}p(y_j^{(\alpha)})\omega_j^{(\alpha)}
$$

对所有

$$
p\in\mathbb P_{2Q+1}
$$

精确。

编程实现建议直接调用标准库：

- Python: `scipy.special.roots_genlaguerre(Q+1, alpha)`
- MATLAB: 可使用 Golub-Welsch 算法或 Laguerre-Gauss 节点函数

---

## 4.2 映射到 $x$ 域的求积节点

由

$$
x=e^{-y/(\beta+1)}
$$

得到

$$
\boxed{
x_j
=
x_j^{(\alpha,\beta)}
=
\exp\left(-\frac{y_j^{(\alpha)}}{\beta+1}\right),
\qquad j=0,\ldots,Q.
}
$$

---

## 4.3 LOF 求积权重

由

$$
(-\log x)^\alpha x^\beta dx
=
\frac1{(\beta+1)^{\alpha+1}}y^\alpha e^{-y}dy
$$

得到 LOF 权重

$$
\boxed{
\chi_j^{(\alpha,\beta)}
=
(\beta+1)^{-\alpha-1}\omega_j^{(\alpha)}.
}
$$

因此

$$
\int_0^1 h(x)(-\log x)^\alpha x^\beta\,dx
\approx
\sum_{j=0}^{Q}h(x_j)\chi_j^{(\alpha,\beta)}.
$$

---

## 4.4 GLOF 求积权重

GLOF 权函数为

$$
\chi^{(\alpha,\lambda)}(x)=(-\log x)^\alpha x^\lambda.
$$

由于

$$
(-\log x)^\alpha x^\lambda dx
=
x^{\lambda-\beta}
(-\log x)^\alpha x^\beta dx,
$$

所以 GLOF 权重为

$$
\boxed{
\chi_j^{(\alpha,\beta,\lambda)}
=
x_j^{\lambda-\beta}\chi_j^{(\alpha,\beta)}.
}
$$

即

$$
\boxed{
\chi_j^{(\alpha,\beta,\lambda)}
=
x_j^{\lambda-\beta}
(\beta+1)^{-\alpha-1}
\omega_j^{(\alpha)}.
}
$$

于是

$$
\int_0^1 h(x)(-\log x)^\alpha x^\lambda\,dx
\approx
\sum_{j=0}^{Q}h(x_j)\chi_j^{(\alpha,\beta,\lambda)}.
$$

---

# 5. 从 Fourier 系数计算函数值

已知

$$
\{\hat f_n\}_{n=-N_F}^{N_F}.
$$

在任意节点 $x_j$ 上计算 Fourier 部分和

$$
\boxed{
f_{N_F}(x_j)
=
\sum_{n=-N_F}^{N_F}
\hat f_n e^{2\pi i n x_j}.
}
$$

矩阵形式为

$$
\mathbf f_N=E\hat{\mathbf f},
$$

其中

$$
\mathbf f_N
=
\begin{bmatrix}
f_{N_F}(x_0)\\
f_{N_F}(x_1)\\
\vdots\\
f_{N_F}(x_Q)
\end{bmatrix},
$$

$$
\hat{\mathbf f}
=
\begin{bmatrix}
\hat f_{-N_F}\\
\hat f_{-N_F+1}\\
\vdots\\
\hat f_{N_F}
\end{bmatrix},
$$

并且

$$
E_{j,n}
=
e^{2\pi i n x_j},
\qquad
j=0,\ldots,Q,
\qquad
n=-N_F,\ldots,N_F.
$$

⚠️ 由于 $x_j$ 是非均匀节点，不能直接用标准 FFT；可直接矩阵乘法，复杂度为 $O(QN_F)$。当 $N_F,Q$ 很大时，可考虑 NUFFT。

---

# 6. 构造 GLOF 基函数矩阵

在求积节点 $x_j$ 上计算

$$
y_j=-(\beta+1)\log x_j.
$$

由于 $x_j$ 本身由 $y_j^{(\alpha)}$ 映射得到，理论上有

$$
y_j=y_j^{(\alpha)}.
$$

对 $\ell=0,\ldots,M$，定义基函数矩阵

$$
\Phi_{j\ell}
=
\mathcal S_\ell^{(\alpha,\beta,\lambda)}(x_j).
$$

即

$$
\boxed{
\Phi_{j\ell}
=
x_j^\nu L_\ell^{(\alpha)}(y_j),
\qquad
\nu=\frac{\beta-\lambda}{2}.
}
$$

若采用奇性匹配参数 $\beta=\lambda+2r$，则

$$
\nu=r,
$$

从而

$$
\Phi_{j\ell}
=
x_j^r L_\ell^{(\alpha)}(y_j).
$$

---

# 7. 计算 GLOF 系数

## 7.1 精确 GLOF 系数

理想 GLOF 投影系数定义为

$$
c_\ell
=
\frac1{\gamma_\ell^{(\alpha,\beta)}}
\int_0^1
f(x)
\mathcal S_\ell^{(\alpha,\beta,\lambda)}(x)
(-\log x)^\alpha x^\lambda\,dx.
$$

其中

$$
\gamma_\ell^{(\alpha,\beta)}
=
\frac{\Gamma(\ell+\alpha+1)}
{(\beta+1)^{\alpha+1}\Gamma(\ell+1)}.
$$

---

## 7.2 由 Fourier 部分和得到的近似 GLOF 系数

用 $f_{N_F}$ 替代 $f$，定义

$$
\boxed{
\tilde c_{\ell,N_F}
=
\frac1{\gamma_\ell^{(\alpha,\beta)}}
\int_0^1
f_{N_F}(x)
\mathcal S_\ell^{(\alpha,\beta,\lambda)}(x)
(-\log x)^\alpha x^\lambda\,dx.
}
$$

用 Gauss-GLOF 求积近似：

$$
\boxed{
\tilde c_{\ell,N_F}^{(Q)}
=
\frac1{\gamma_\ell^{(\alpha,\beta)}}
\sum_{j=0}^{Q}
f_{N_F}(x_j)
\Phi_{j\ell}
\chi_j^{(\alpha,\beta,\lambda)}.
}
$$

其中

$$
\Phi_{j\ell}
=
\mathcal S_\ell^{(\alpha,\beta,\lambda)}(x_j),
$$

$$
\chi_j^{(\alpha,\beta,\lambda)}
=
x_j^{\lambda-\beta}
(\beta+1)^{-\alpha-1}
\omega_j^{(\alpha)}.
$$

---

## 7.3 矩阵形式

定义对角权重矩阵

$$
W=\operatorname{diag}\left(
\chi_0^{(\alpha,\beta,\lambda)},
\chi_1^{(\alpha,\beta,\lambda)},
\ldots,
\chi_Q^{(\alpha,\beta,\lambda)}
\right).
$$

定义归一化矩阵

$$
G=\operatorname{diag}\left(
\gamma_0^{(\alpha,\beta)},
\gamma_1^{(\alpha,\beta)},
\ldots,
\gamma_M^{(\alpha,\beta)}
\right).
$$

则

$$
\tilde{\mathbf c}
=
G^{-1}\Phi^T W\mathbf f_N.
$$

其中

$$
\tilde{\mathbf c}
=
\begin{bmatrix}
\tilde c_{0,N_F}\\
\tilde c_{1,N_F}\\
\vdots\\
\tilde c_{M,N_F}
\end{bmatrix}.
$$

如果直接由 Fourier 系数到 GLOF 系数，可以写成

$$
\tilde{\mathbf c}
=
B\hat{\mathbf f},
$$

其中

$$
B
=
G^{-1}\Phi^T W E.
$$

即

$$
\boxed{
B_{\ell n}
=
\frac1{\gamma_\ell^{(\alpha,\beta)}}
\sum_{j=0}^{Q}
e^{2\pi i n x_j}
\mathcal S_\ell^{(\alpha,\beta,\lambda)}(x_j)
\chi_j^{(\alpha,\beta,\lambda)}.
}
$$

于是

$$
\boxed{
\tilde c_{\ell,N_F}
=
\sum_{n=-N_F}^{N_F}
B_{\ell n}\hat f_n.
}
$$

---

# 8. 重构函数值

给定任意评价点

$$
x_{\mathrm{eval}}\in(0,1),
$$

先计算

$$
y_{\mathrm{eval}}
=
-(\beta+1)\log x_{\mathrm{eval}}.
$$

然后计算

$$
\Psi_{\ell}(x_{\mathrm{eval}})
=
\mathcal S_\ell^{(\alpha,\beta,\lambda)}(x_{\mathrm{eval}})
=
x_{\mathrm{eval}}^\nu
L_\ell^{(\alpha)}(y_{\mathrm{eval}}).
$$

最终重构值为

$$
\boxed{
R_{N_F,M}f(x_{\mathrm{eval}})
=
\sum_{\ell=0}^{M}
\tilde c_{\ell,N_F}
\Psi_\ell(x_{\mathrm{eval}}).
}
$$

若有多个评价点 $\{z_i\}_{i=0}^{N_{\mathrm{eval}}}$，定义评价矩阵

$$
\Psi_{i\ell}
=
z_i^\nu
L_\ell^{(\alpha)}
\left(-(\beta+1)\log z_i\right).
$$

则

$$
\boxed{
\mathbf R
=
\Psi\tilde{\mathbf c}.
}
$$

---

# 9. 完整算法流程

## Algorithm 1: Fourier-GLOF Reconstruction

### 输入

1. Fourier 系数：

$$
\hat f_n,\qquad n=-N_F,\ldots,N_F.
$$

2. 奇性参数：

$$
r>0,
\qquad
k\in\mathbb N_0.
$$

3. GLOF 参数：

$$
\alpha>-1,
\qquad
\lambda\in\mathbb R,
\qquad
\beta=\lambda+2r.
$$

4. 重构阶数：

$$
M.
$$

5. 求积阶数：

$$
Q.
$$

6. 评价点：

$$
z_i\in(0,1),
\qquad i=0,\ldots,N_{\mathrm{eval}}.
$$

---

### Step 1: 检查参数

检查

$$
\beta>-1,
$$

$$
2r+\lambda>-1,
$$

$$
M\le Q
$$

或更保守地取

$$
Q\ge 2M+N_{\mathrm{buffer}}.
$$

若由 Fourier 系数计算，建议

$$
Q\ge \max(2M,2N_F)
$$

作为初始选择。

⚠️ 如果 $f_{N_F}(x)$ 在求积节点附近振荡强烈，需增大 $Q$。

---

### Step 2: 生成 Laguerre-Gauss 节点与权重

计算

$$
\{y_j^{(\alpha)},\omega_j^{(\alpha)}\}_{j=0}^{Q},
$$

满足

$$
\int_0^\infty p(y)y^\alpha e^{-y}\,dy
=
\sum_{j=0}^{Q}p(y_j^{(\alpha)})\omega_j^{(\alpha)}
$$

对所有 $p\in\mathbb P_{2Q+1}$ 精确。

---

### Step 3: 映射到 $x$ 域节点

计算

$$
x_j
=
\exp\left(-\frac{y_j^{(\alpha)}}{\beta+1}\right),
\qquad j=0,\ldots,Q.
$$

---

### Step 4: 计算 GLOF 求积权重

先计算

$$
\chi_j^{(\alpha,\beta)}
=
(\beta+1)^{-\alpha-1}\omega_j^{(\alpha)}.
$$

再计算

$$
\chi_j^{(\alpha,\beta,\lambda)}
=
x_j^{\lambda-\beta}
\chi_j^{(\alpha,\beta)}.
$$

即

$$
\chi_j^{(\alpha,\beta,\lambda)}
=
x_j^{\lambda-\beta}
(\beta+1)^{-\alpha-1}
\omega_j^{(\alpha)}.
$$

---

### Step 5: 在求积节点上计算 Fourier 部分和

对每个 $j=0,\ldots,Q$，计算

$$
f_j
=
f_{N_F}(x_j)
=
\sum_{n=-N_F}^{N_F}
\hat f_n e^{2\pi i n x_j}.
$$

矩阵形式为

$$
f_j=\sum_{n=-N_F}^{N_F}E_{j n}\hat f_n,
$$

其中

$$
E_{j n}=e^{2\pi i n x_j}.
$$

---

### Step 6: 计算 GLOF 基函数矩阵

定义

$$
\nu=\frac{\beta-\lambda}{2}.
$$

在奇性匹配情形，

$$
\nu=r.
$$

对 $j=0,\ldots,Q$ 和 $\ell=0,\ldots,M$，计算

$$
\Phi_{j\ell}
=
x_j^\nu
L_\ell^{(\alpha)}(y_j^{(\alpha)}).
$$

其中 $L_\ell^{(\alpha)}$ 用三项递推计算：

$$
L_0^{(\alpha)}(y)=1,
$$

$$
L_1^{(\alpha)}(y)=-y+\alpha+1,
$$

$$
L_{\ell+1}^{(\alpha)}(y)
=
\frac{2\ell+\alpha+1-y}{\ell+1}L_\ell^{(\alpha)}(y)
-
\frac{\ell+\alpha}{\ell+1}L_{\ell-1}^{(\alpha)}(y).
$$

---

### Step 7: 计算归一化常数

对 $\ell=0,\ldots,M$，计算

$$
\gamma_\ell^{(\alpha,\beta)}
=
\frac{\Gamma(\ell+\alpha+1)}
{(\beta+1)^{\alpha+1}\Gamma(\ell+1)}.
$$

数值实现中建议使用 log-gamma 避免溢出：

$$
\log \gamma_\ell^{(\alpha,\beta)}
=
\log\Gamma(\ell+\alpha+1)
-
\log\Gamma(\ell+1)
-
(\alpha+1)\log(\beta+1).
$$

因此

$$
\gamma_\ell^{(\alpha,\beta)}
=
\exp\left(
\log \gamma_\ell^{(\alpha,\beta)}
\right).
$$

---

### Step 8: 计算近似 GLOF 系数

对每个 $\ell=0,\ldots,M$，计算

$$
\tilde c_{\ell,N_F}
=
\frac1{\gamma_\ell^{(\alpha,\beta)}}
\sum_{j=0}^{Q}
f_j
\Phi_{j\ell}
\chi_j^{(\alpha,\beta,\lambda)}.
$$

矩阵形式：

$$
\tilde{\mathbf c}
=
G^{-1}\Phi^T W\mathbf f.
$$

其中

$$
G_{\ell\ell}
=
\gamma_\ell^{(\alpha,\beta)},
$$

$$
W_{jj}
=
\chi_j^{(\alpha,\beta,\lambda)}.
$$

---

### Step 9: 在评价点上重构

对每个评价点 $z_i\in(0,1)$，计算

$$
Y_i=-(\beta+1)\log z_i.
$$

然后构造

$$
\Psi_{i\ell}
=
z_i^\nu L_\ell^{(\alpha)}(Y_i),
\qquad
\ell=0,\ldots,M.
$$

最终重构值为

$$
R_i
=
R_{N_F,M}f(z_i)
=
\sum_{\ell=0}^{M}
\tilde c_{\ell,N_F}\Psi_{i\ell}.
$$

矩阵形式：

$$
\mathbf R=\Psi\tilde{\mathbf c}.
$$

# 11. 可直接实现的矩阵版本

## 11.1 Fourier 评价矩阵

令 Fourier 模式集合为

$$
\mathcal N_F=\{-N_F,-N_F+1,\ldots,N_F\}.
$$

定义

$$
E_{j n}
=
e^{2\pi i n x_j},
\qquad
j=0,\ldots,Q,
\qquad
n\in\mathcal N_F.
$$

则

$$
\mathbf f=E\hat{\mathbf f}.
$$

---

## 11.2 GLOF 求积矩阵

定义

$$
\Phi_{j\ell}
=
x_j^\nu L_\ell^{(\alpha)}(y_j),
\qquad
j=0,\ldots,Q,
\qquad
\ell=0,\ldots,M.
$$

定义

$$
W=\operatorname{diag}(\chi_0,\ldots,\chi_Q),
$$

其中

$$
\chi_j=\chi_j^{(\alpha,\beta,\lambda)}.
$$

定义

$$
G=\operatorname{diag}(\gamma_0,\ldots,\gamma_M),
$$

其中

$$
\gamma_\ell=\gamma_\ell^{(\alpha,\beta)}.
$$

则从节点函数值到 GLOF 系数的矩阵为

$$
A_{\mathrm{GLOF}}
=
G^{-1}\Phi^T W.
$$

即

$$
\tilde{\mathbf c}
=
A_{\mathrm{GLOF}}\mathbf f.
$$

---

## 11.3 Fourier 系数到 GLOF 系数矩阵

由

$$
\mathbf f=E\hat{\mathbf f}
$$

可得

$$
\boxed{
\tilde{\mathbf c}
=
G^{-1}\Phi^TWE\hat{\mathbf f}.
}
$$

定义

$$
B=G^{-1}\Phi^TWE.
$$

则

$$
\boxed{
\tilde{\mathbf c}=B\hat{\mathbf f}.
}
$$

其中

$$
B_{\ell n}
=
\frac1{\gamma_\ell}
\sum_{j=0}^{Q}
x_j^\nu L_\ell^{(\alpha)}(y_j)
e^{2\pi i n x_j}
\chi_j.
$$

---

## 11.4 评价矩阵

对评价点 $\{z_i\}_{i=0}^{N_{\mathrm{eval}}}$，定义

$$
Y_i=-(\beta+1)\log z_i.
$$

定义

$$
\Psi_{i\ell}
=
z_i^\nu L_\ell^{(\alpha)}(Y_i).
$$

则最终重构为

$$
\boxed{
\mathbf R=\Psi B\hat{\mathbf f}.
}
$$

这也是最紧凑的程序实现形式。

---

# 12. 误差诊断公式

## 12.1 正则化误差

精确 GLOF 投影为

$$
P_Mf(x)
=
\sum_{\ell=0}^{M}
c_\ell
\mathcal S_\ell^{(\alpha,\beta,\lambda)}(x).
$$

正则化误差定义为

$$
RE_M
=
\|f-P_Mf\|.
$$

对模型函数

$$
f_0(x)=x^r(-\log x)^k
$$

可用理论比率

$$
R_{\alpha,\beta,\lambda}
=
\left|
\frac{2r+\lambda-\beta}
{2r+2+\lambda+\beta}
\right|.
$$

若

$$
R_{\alpha,\beta,\lambda}<1,
$$

则有典型估计

$$
\|f_0-P_Mf_0\|_{\chi^{(\alpha,\lambda)}}
\lesssim
C(k+1)!
M^{\frac{\alpha+1}{2}+k}
R_{\alpha,\beta,\lambda}^{M}.
$$

在匹配参数

$$
\beta=\lambda+2r
$$

下，

$$
R_{\alpha,\beta,\lambda}=0.
$$

这说明主奇异项 $x^r(-\log x)^k$ 被 GLOF 空间精确捕捉，实际误差主要来自解析因子 $g(x)$ 的非平凡变化。

---

## 12.2 截断误差

Fourier 数据重构误差分解为

$$
f-R_{N_F,M}f
=
(f-P_Mf)
+
(P_Mf-R_{N_F,M}f).
$$

定义

$$
TE_{N_F,M}
=
\|P_Mf-R_{N_F,M}f\|.
$$

系数误差为

$$
\delta c_{\ell,N_F}
=
c_\ell-\tilde c_{\ell,N_F}.
$$

由定义，

$$
\delta c_{\ell,N_F}
=
\frac1{\gamma_\ell}
\int_0^1
(f(x)-f_{N_F}(x))
\mathcal S_\ell^{(\alpha,\beta,\lambda)}(x)
\chi^{(\alpha,\lambda)}(x)\,dx.
$$

利用 Fourier 尾项

$$
f(x)-f_{N_F}(x)
=
\sum_{|n|>N_F}\hat f_n e^{2\pi i n x},
$$

得到

$$
\delta c_{\ell,N_F}
=
\sum_{|n|>N_F}\hat f_n I_{\ell n},
$$

其中

$$
I_{\ell n}
=
\frac1{\gamma_\ell}
\int_0^1
e^{2\pi i n x}
K_\ell(x)\,dx,
$$

核函数定义为

$$
K_\ell(x)
=
\mathcal S_\ell^{(\alpha,\beta,\lambda)}(x)
\chi^{(\alpha,\lambda)}(x).
$$

展开后

$$
K_\ell(x)
=
x^\nu
L_\ell^{(\alpha)}\left(-(\beta+1)\log x\right)
(-\log x)^\alpha x^\lambda.
$$

即

$$
\boxed{
K_\ell(x)
=
x^{\lambda+\nu}
(-\log x)^\alpha
L_\ell^{(\alpha)}
\left(-(\beta+1)\log x\right).
}
$$

若采用匹配参数 $\nu=r$，则

$$
K_\ell(x)
=
x^{\lambda+r}
(-\log x)^\alpha
L_\ell^{(\alpha)}
\left(-(\lambda+2r+1)\log x\right).
$$

当 $\lambda\sim N_F$ 时，$x=0$ 端点具有高阶消失。

---

## 12.3 分部积分型截断误差机制

若 $K_\ell$ 及其前 $t-1$ 阶导数在端点处边界项消失，则

$$
I_{\ell n}
=
\frac{1}{(2\pi i n)^t}
\frac1{\gamma_\ell}
\int_0^1
e^{2\pi i n x}
K_\ell^{(t)}(x)\,dx.
$$

因此

$$
|I_{\ell n}|
\le
\frac1{\gamma_\ell}
(2\pi |n|)^{-t}
\|K_\ell^{(t)}\|_{L^1(0,1)}.
$$

于是

$$
|\delta c_{\ell,N_F}|
\le
\sum_{|n|>N_F}
|\hat f_n|
\frac1{\gamma_\ell}
(2\pi |n|)^{-t}
\|K_\ell^{(t)}\|_{L^1}.
$$

若仅使用粗略有界性

$$
|\hat f_n|\le C,
$$

则

$$
\sum_{|n|>N_F}|n|^{-t}
\lesssim
N_F^{-t+1}.
$$

于是

$$
|\delta c_{\ell,N_F}|
\lesssim
\frac{N_F^{-t+1}}{\gamma_\ell}
\|K_\ell^{(t)}\|_{L^1}.
$$

若可取

$$
t\sim cN_F
$$

并证明

$$
\|K_\ell^{(t)}\|_{L^1}
\le
C^t(\lambda+\ell)^t
\times
\text{moderate factors},
$$

则可能得到

$$
TE_{N_F,M}
\le
C N_F^p q_T^{N_F},
\qquad 0<q_T<1.
$$

⚠️ 该 Fourier-GLOF 截断误差界不是直接由标准 GLOF 投影理论自动推出的。关键困难是 $K_\ell$ 在 $x=1$ 处是否也有足够高阶的边界消失。

---

# 13. 端点处理与稳定性注意事项

## 13.1 不要直接在 $x=0$ 处评价

因为

$$
y(x)=-(\beta+1)\log x\to+\infty,
\qquad x\to0^+.
$$

程序中不应直接令 $x=0$ 代入。

若 $\nu=r>0$，则

$$
x^r L_\ell^{(\alpha)}(-(\beta+1)\log x)\to0,
\qquad x\to0^+.
$$

因此可定义

$$
R_{N_F,M}f(0)=0.
$$

⚠️ 若 $\nu=0$ 或 $r=0$，则端点极限可能不存在或无界，应避免在 $x=0$ 处评价。

---

## 13.2 右端点 $x=1$ 的处理

在 $x=1$ 附近，

$$
-\log x\sim 1-x.
$$

若 $\alpha$ 固定，则

$$
(-\log x)^\alpha\sim(1-x)^\alpha.
$$

这只提供固定阶右端点消失，不足以支撑 $O(N_F)$ 次分部积分。

因此若需要严格的 Fourier 后处理指数误差机制，可考虑：

### 方案 A：令 $\alpha$ 随 $N_F$ 增长

$$
\alpha=\alpha_{N_F}=\eta N_F.
$$

此时

$$
(-\log x)^\alpha\sim(1-x)^{\eta N_F}.
$$

⚠️ 这会改变 Laguerre 大参数性质，需重新检查数值稳定性和误差证明。

### 方案 B：右端点窗口

引入窗口函数

$$
w_R(x)=(1-x)^\mu,
\qquad \mu\sim N_F.
$$

改用核

$$
\widetilde K_\ell(x)
=
w_R(x)K_\ell(x).
$$

⚠️ 该方案会破坏标准 GLOF 正交性，需要重新定义投影或采用加权最小二乘。

### 方案 C：周期匹配或端点修正

若 $f(1)$、$f'(1)$ 等可获得，可减去光滑修正函数 $p(x)$，使残差满足

$$
\partial_x^j(f-p)(1)=0,
\qquad j=0,\ldots,J.
$$

然后对残差重构。

---

## 13.3 数值溢出与下溢

当 $\lambda$、$\beta$ 或 $Q$ 较大时，节点 $x_j$ 可能极小：

$$
x_j=\exp\left(-\frac{y_j}{\beta+1}\right).
$$

建议：

1. 使用 log-gamma 计算 $\gamma_\ell$；
2. 使用递推计算 Laguerre 多项式；
3. 对 $x_j^\nu$ 使用

$$
x_j^\nu=\exp(\nu\log x_j);
$$

4. 若 $\nu\log x_j$ 很小导致下溢，可使用缩放基函数；
5. 对高阶 Laguerre 多项式可考虑正交归一化版本。

---

# 14. 推荐实现函数接口

## 14.1 主函数

函数名：

$$
\texttt{fourier\_glof\_reconstruct}
$$

输入：

$$
\{\hat f_n\}_{n=-N_F}^{N_F},
\quad
r,
\quad
k,
\quad
\alpha,
\quad
\theta,
\quad
\gamma,
\quad
Q,
\quad
\{z_i\}.
$$

内部参数：

$$
\lambda=\lfloor \theta N_F\rfloor,
$$

$$
\beta=\lambda+2r,
$$

$$
M=\lfloor \gamma N_F\rfloor,
$$

$$
\nu=r.
$$

输出：

$$
\{R_{N_F,M}f(z_i)\}_{i=0}^{N_{\mathrm{eval}}}.
$$

---

## 14.2 辅助函数 1：Laguerre-Gauss 节点

输入：

$$
Q,\alpha.
$$

输出：

$$
\{y_j,\omega_j\}_{j=0}^{Q}.
$$

要求：

$$
\int_0^\infty p(y)y^\alpha e^{-y}\,dy
=
\sum_{j=0}^{Q}p(y_j)\omega_j
$$

对 $p\in\mathbb P_{2Q+1}$ 精确。

---

## 14.3 辅助函数 2：GLOF 基函数矩阵

输入：

$$
\{x_j\},
\quad
\{y_j\},
\quad
\alpha,
\quad
\beta,
\quad
\lambda,
\quad
M.
$$

输出：

$$
\Phi_{j\ell}
=
x_j^{(\beta-\lambda)/2}L_\ell^{(\alpha)}(y_j).
$$

---

## 14.4 辅助函数 3：Fourier 部分和

输入：

$$
\{\hat f_n\}_{n=-N_F}^{N_F},
\quad
\{x_j\}.
$$

输出：

$$
f_j
=
\sum_{n=-N_F}^{N_F}\hat f_n e^{2\pi i n x_j}.
$$

---

## 14.5 辅助函数 4：GLOF 系数

输入：

$$
\{f_j\},
\quad
\Phi,
\quad
\{\chi_j\},
\quad
\{\gamma_\ell\}.
$$

输出：

$$
\tilde c_\ell
=
\frac1{\gamma_\ell}
\sum_{j=0}^{Q}
f_j\Phi_{j\ell}\chi_j.
$$

---

## 14.6 辅助函数 5：评价重构

输入：

$$
\{\tilde c_\ell\}_{\ell=0}^{M},
\quad
\{z_i\},
\quad
\alpha,
\quad
\beta,
\quad
\lambda.
$$

输出：

$$
R_i
=
\sum_{\ell=0}^{M}
\tilde c_\ell
z_i^{(\beta-\lambda)/2}
L_\ell^{(\alpha)}
\left(-(\beta+1)\log z_i\right).
$$

---

# 15. 最小可运行数学流程汇总

给定

$$
\hat f_n,\quad |n|\le N_F.
$$

选择

$$
M=\lfloor\gamma N_F\rfloor,
\qquad
\lambda=\lfloor\theta N_F\rfloor,
\qquad
\beta=\lambda+2r,
\qquad
\nu=r.
$$

生成 Laguerre-Gauss 节点和权重：

$$
\{y_j,\omega_j\}_{j=0}^{Q}.
$$

映射节点：

$$
x_j=e^{-y_j/(\beta+1)}.
$$

计算权重：

$$
\chi_j=x_j^{\lambda-\beta}(\beta+1)^{-\alpha-1}\omega_j.
$$

计算 Fourier 部分和：

$$
f_j=\sum_{n=-N_F}^{N_F}\hat f_n e^{2\pi i n x_j}.
$$

计算基函数矩阵：

$$
\Phi_{j\ell}=x_j^r L_\ell^{(\alpha)}(y_j).
$$

计算归一化：

$$
\gamma_\ell=
\frac{\Gamma(\ell+\alpha+1)}
{(\beta+1)^{\alpha+1}\Gamma(\ell+1)}.
$$

计算 GLOF 系数：

$$
\tilde c_\ell=
\frac1{\gamma_\ell}
\sum_{j=0}^{Q}f_j\Phi_{j\ell}\chi_j.
$$

在评价点 $z$ 上计算：

$$
Y(z)=-(\beta+1)\log z,
$$

$$
\Psi_\ell(z)=z^rL_\ell^{(\alpha)}(Y(z)).
$$

输出重构：

$$
\boxed{
R_{N_F,M}f(z)
=
\sum_{\ell=0}^{M}
\tilde c_\ell z^r
L_\ell^{(\alpha)}
\left(-(\beta+1)\log z\right).
}
$$

---

# 16. 实验建议

## 16.1 测试函数

可先测试模型函数

$$
f(x)=x^r(-\log x)^k.
$$

再测试

$$
f(x)=x^r(-\log x)^k e^x.
$$

或

$$
f(x)=x^r(-\log x)^k\sin x.
$$

---

## 16.2 误差计算

若精确函数已知，在测试网格 $\{z_i\}$ 上计算

$$
E_\infty
=
\max_i
|f(z_i)-R_{N_F,M}f(z_i)|.
$$

也可计算加权离散 $L^2$ 误差

$$
E_{2,\chi}
=
\left(
\sum_i
|f(z_i)-R_{N_F,M}f(z_i)|^2
w_i
\right)^{1/2}.
$$

---

## 16.3 无精确解时的收敛诊断

使用双分辨率比较：

$$
E_{\mathrm{proxy}}
=
\max_i
\left|
R_{N_F,M}f(z_i)
-
R_{2N_F,2M}f(z_i)
\right|.
$$

或比较相邻 GLOF 阶数：

$$
E_{\mathrm{tail}}
=
\left(
\sum_{\ell=M-p}^{M}
|\tilde c_\ell|^2\gamma_\ell
\right)^{1/2}.
$$

其中 $p$ 可取

$$
p=5,\ldots,10.
$$

---

# 17. 关键注意事项总结

1. GLOF 的自然变量是

$$
y=-(\beta+1)\log x.
$$

2. GLOF 的自然正交基是 Laguerre 多项式，而不是 Gegenbauer 或 Legendre 多项式。

3. 为匹配

$$
x^r(-\log x)^k
$$

的幂奇性，应选择

$$
\beta=\lambda+2r.
$$

4. 从 Fourier 系数到 GLOF 系数的核心矩阵为

$$
B=G^{-1}\Phi^TWE.
$$

5. 重构公式为

$$
R_{N_F,M}f(x)
=
\sum_{\ell=0}^{M}
\tilde c_{\ell,N_F}
x^r
L_\ell^{(\alpha)}
\left(-(\lambda+2r+1)\log x\right).
$$

6. 直接在 $x=0$ 处评价是数值不稳定的；若 $r>0$，可用极限值

$$
R_{N_F,M}f(0)=0.
$$

7. 若 Fourier 周期延拓在 $x=1$ 处有跳跃，必须先处理右端点，否则重构会被远端 Gibbs 现象污染。

8. 标准 GLOF 理论保证对一端点弱奇性有良好逼近；但从有限 Fourier 系数恢复指数精度的完整截断误差证明，还需要额外处理 $x=1$ 的边界项。

# 改进算法：通过直接计算 $B_{\ell n}$ 的 Fourier-to-GLOF 系数空间重构

# 3. Fourier-to-GLOF 矩阵 $B_{\ell n}$

## 3.1 原始矩阵元定义

定义

$$
\boxed{
B_{\ell n}
=
\frac1{\gamma_\ell^{(\alpha,\beta)}}
\int_0^1
e^{2\pi i n x}
\mathcal S_\ell^{(\alpha,\beta,\lambda)}(x)
\chi^{(\alpha,\lambda)}(x)\,dx.
}
$$

展开 GLOF 基函数和权函数：

$$
B_{\ell n}
=
\frac1{\gamma_\ell^{(\alpha,\beta)}}
\int_0^1
e^{2\pi i n x}
x^\nu
L_\ell^{(\alpha)}\left(-(\beta+1)\log x\right)
(-\log x)^\alpha x^\lambda\,dx.
$$

即

$$
\boxed{
B_{\ell n}
=
\frac1{\gamma_\ell^{(\alpha,\beta)}}
\int_0^1
e^{2\pi i n x}
x^{\lambda+\nu}
(-\log x)^\alpha
L_\ell^{(\alpha)}\left(-(\beta+1)\log x\right)
\,dx.
}
$$

匹配参数下，$\nu=r$，所以

$$
\boxed{
B_{\ell n}
=
\frac1{\gamma_\ell^{(\alpha,\beta)}}
\int_0^1
e^{2\pi i n x}
x^{\lambda+r}
(-\log x)^\alpha
L_\ell^{(\alpha)}\left(-(\lambda+2r+1)\log x\right)
\,dx.
}
$$

该公式不需要在极小 GLOF 节点上计算 $f_{N_F}(x_j)$。

---

## 3.2 原始矩阵元与未归一化矩阵元

为了后续边界适配基函数方便，定义未归一化矩阵元

$$
D_{\ell n}
=
\int_0^1
e^{2\pi i n x}
\mathcal S_\ell^{(\alpha,\beta,\lambda)}(x)
\chi^{(\alpha,\lambda)}(x)\,dx.
$$

则

$$
\boxed{
B_{\ell n}
=
\frac{D_{\ell n}}
{\gamma_\ell^{(\alpha,\beta)}}.
}
$$

即

$$
D_{\ell n}
=
\gamma_\ell^{(\alpha,\beta)}B_{\ell n}.
$$

---

# 4. $B_{\ell n}$ 的稳定计算公式

## 4.1 $y$ 域表达式

使用

$$
y=-(\beta+1)\log x,
\qquad
x=e^{-y/(\beta+1)}.
$$

有

$$
dx
=
\frac1{\beta+1}e^{-y/(\beta+1)}\,dy.
$$

又因为

$$
x^{\lambda+\nu}
=
e^{-(\lambda+\nu)y/(\beta+1)},
$$

以及

$$
(-\log x)^\alpha
=
\left(\frac{y}{\beta+1}\right)^\alpha,
$$

所以

$$
D_{\ell n}
=
\frac1{(\beta+1)^{\alpha+1}}
\int_0^\infty
e^{2\pi i n e^{-y/(\beta+1)}}
L_\ell^{(\alpha)}(y)
y^\alpha
e^{-\sigma y}
\,dy,
$$

其中

$$
\boxed{
\sigma
=
\frac{\lambda+\nu+1}{\beta+1}.
}
$$

由于

$$
\nu=\frac{\beta-\lambda}{2},
$$

可写成

$$
\boxed{
\sigma
=
\frac{\beta+\lambda+2}{2(\beta+1)}.
}
$$

匹配参数 $\beta=\lambda+2r$ 下，

$$
\boxed{
\sigma
=
\frac{\lambda+r+1}{\lambda+2r+1}.
}
$$

因此

$$
\boxed{
B_{\ell n}
=
\frac1{\gamma_\ell^{(\alpha,\beta)}(\beta+1)^{\alpha+1}}
\int_0^\infty
e^{2\pi i n e^{-y/(\beta+1)}}
L_\ell^{(\alpha)}(y)
y^\alpha
e^{-\sigma y}
\,dy.
}
$$

---

## 4.2 缩放到标准 Laguerre 权重

若直接用标准 Laguerre-Gauss 权重 $y^\alpha e^{-y}$，当 $\sigma<1$ 时会出现额外因子 $e^{(1-\sigma)y}$，数值上可能不稳定。

因此推荐使用缩放变量

$$
z=\sigma y.
$$

于是

$$
y=\frac z\sigma,
\qquad
dy=\frac{dz}{\sigma}.
$$

并且

$$
y^\alpha e^{-\sigma y}dy
=
\sigma^{-\alpha-1}
z^\alpha e^{-z}dz.
$$

所以

$$
D_{\ell n}
=
\frac1{(\beta+1)^{\alpha+1}\sigma^{\alpha+1}}
\int_0^\infty
e^{2\pi i n \exp\left(-\frac{z}{\sigma(\beta+1)}\right)}
L_\ell^{(\alpha)}\left(\frac z\sigma\right)
z^\alpha e^{-z}
\,dz.
$$

因此

$$
\boxed{
B_{\ell n}
=
\frac1{\gamma_\ell^{(\alpha,\beta)}(\beta+1)^{\alpha+1}\sigma^{\alpha+1}}
\int_0^\infty
e^{2\pi i n \exp\left(-\frac{z}{\sigma(\beta+1)}\right)}
L_\ell^{(\alpha)}\left(\frac z\sigma\right)
z^\alpha e^{-z}
\,dz.
}
$$

这是推荐用于程序实现的公式。

---

## 4.3 Laguerre-Gauss 求积计算 $B_{\ell n}$

取 $Q_B+1$ 个 Laguerre-Gauss 节点和权重

$$
\{z_j,\omega_j^{(\alpha)}\}_{j=0}^{Q_B},
$$

满足

$$
\int_0^\infty p(z)z^\alpha e^{-z}\,dz
=
\sum_{j=0}^{Q_B}p(z_j)\omega_j^{(\alpha)}
$$

对 $p\in\mathbb P_{2Q_B+1}$ 精确。

定义

$$
x_j^{(B)}
=
\exp\left(
-\frac{z_j}{\sigma(\beta+1)}
\right).
$$

则

$$
D_{\ell n}^{(Q_B)}
=
\frac1{(\beta+1)^{\alpha+1}\sigma^{\alpha+1}}
\sum_{j=0}^{Q_B}
\omega_j^{(\alpha)}
e^{2\pi i n x_j^{(B)}}
L_\ell^{(\alpha)}\left(\frac{z_j}{\sigma}\right).
$$

因此

$$
\boxed{
B_{\ell n}^{(Q_B)}
=
\frac1{\gamma_\ell^{(\alpha,\beta)}(\beta+1)^{\alpha+1}\sigma^{\alpha+1}}
\sum_{j=0}^{Q_B}
\omega_j^{(\alpha)}
e^{2\pi i n x_j^{(B)}}
L_\ell^{(\alpha)}\left(\frac{z_j}{\sigma}\right).
}
$$

注意：这里虽然仍然出现很小的 $x_j^{(B)}$，但只需要计算有界相位

$$
e^{2\pi i n x_j^{(B)}}.
$$

不需要计算

$$
f_{N_F}(x_j^{(B)})
$$

或

$$
x_j^{-r}f_{N_F}(x_j).
$$

因此避免了极小节点处的 Fourier 点值污染与放大。

---

# 5. Laguerre 多项式递推

在计算

$$
L_\ell^{(\alpha)}\left(\frac{z_j}{\sigma}\right)
$$

时，可使用三项递推。

令

$$
Y_j=\frac{z_j}{\sigma}.
$$

则

$$
L_0^{(\alpha)}(Y_j)=1,
$$

$$
L_1^{(\alpha)}(Y_j)=-Y_j+\alpha+1,
$$

对 $\ell\ge1$，

$$
L_{\ell+1}^{(\alpha)}(Y_j)
=
\frac{2\ell+\alpha+1-Y_j}{\ell+1}
L_\ell^{(\alpha)}(Y_j)
-
\frac{\ell+\alpha}{\ell+1}
L_{\ell-1}^{(\alpha)}(Y_j).
$$

