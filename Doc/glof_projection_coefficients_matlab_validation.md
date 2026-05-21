# MATLAB 验证 GLOF 投影系数计算正确性的测试指南

本文档用于指导 MATLAB 编写和验证 **GLOF projection coefficients** 的计算程序。目标不是直接测试重构算法本身，而是首先确保如下核心模块正确：

\[
\widehat f_n^{(\alpha,\beta,\lambda)} = \frac{1}{\gamma_n^{(\alpha,\beta)}}
\int_0^1
f(t)S_n^{(\alpha,\beta,\lambda)}(t)
(-\log t)^\alpha t^\lambda\,dt,
\qquad n=0,\ldots,m.
\]

其中 GLOF 系数积分应使用 Chen--Shen 文章中的 **Gauss--LOF / Gauss--GLOF quadrature**，而不是普通的 `integral`、复合梯形或复合 Simpson 积分。普通积分方法容易把端点弱奇异性、节点聚集误差和 GLOF 投影误差混在一起，导致无法判断系数代码是否正确。

---

## 1. 数学定义与 MATLAB 记号

### 1.1 Log mapping

在 Chen--Shen 文章中，基本映射为

\[
y=-(\beta+1)\log t,
\qquad t=e^{-y/(\beta+1)},
\qquad t\in(0,1),\quad y\in(0,\infty).
\]

MATLAB 中建议使用

```matlab
y = ...;                         % Laguerre-Gauss nodes on [0, infinity)
logt = -y/(beta+1);              % log(t)
t = exp(logt);                   % mapped nodes in (0,1)
neglogt = y/(beta+1);            % -log(t)
```

不要用 `-log(t)` 重新计算 `neglogt`，因为当 `t` 很小时会有额外舍入误差。

---

### 1.2 LOF and GLOF

LOF 定义为

\[
S_n^{(\alpha,\beta)}(t) = L_n^{(\alpha)}\bigl(-(\beta+1)\log t\bigr) = L_n^{(\alpha)}(y).
\]

GLOF 定义为

\[
S_n^{(\alpha,\beta,\lambda)}(t) = t^{(\beta-\lambda)/2}S_n^{(\alpha,\beta)}(t) = t^{(\beta-\lambda)/2}L_n^{(\alpha)}(y).
\]

对应权函数为

\[
\chi_{\alpha,\lambda}(t)=(-\log t)^\alpha t^\lambda.
\]

GLOF 正交性为

\[
\int_0^1
S_n^{(\alpha,\beta,\lambda)}(t)
S_m^{(\alpha,\beta,\lambda)}(t)
(-\log t)^\alpha t^\lambda\,dt = \gamma_n^{(\alpha,\beta)}\delta_{nm},
\]

其中

\[
\gamma_n^{(\alpha,\beta)} = \frac{\Gamma(n+\alpha+1)}{(\beta+1)^{\alpha+1}\Gamma(n+1)}.
\]

MATLAB 中应使用 `gammaln` 计算：

```matlab
n = (0:m).';
gamma_n = exp(gammaln(n + alpha + 1) - gammaln(n + 1)) / (beta + 1)^(alpha + 1);
```

---

## 2. Gauss--LOF / Gauss--GLOF 求积

### 2.1 Laguerre--Gauss nodes and weights

令 \(\{y_j^{(\alpha)},\omega_j^{(\alpha)}\}_{j=0}^Q\) 为 \(Q+1\) 点 generalized Laguerre--Gauss 求积节点和权重：

\[
\int_0^\infty p(y)y^\alpha e^{-y}\,dy = \sum_{j=0}^Q p(y_j^{(\alpha)})\omega_j^{(\alpha)},
\qquad p\in \mathbb P_{2Q+1}.
\]

建议用 Golub--Welsch 方法计算节点与权重。设节点数为 `nq = Q+1`，Jacobi 矩阵为

\[
J_{ii}=2i+\alpha+1,
\qquad i=0,\ldots,Q,
\]

\[
J_{i,i-1}=J_{i-1,i}=\sqrt{i(i+\alpha)},
\qquad i=1,\ldots,Q.
\]

MATLAB 函数：

```matlab
function [y,w] = laguerre_gauss_gw(Q, alpha)
% Return Q+1 generalized Laguerre-Gauss nodes and weights for
% integral_0^inf p(y) y^alpha exp(-y) dy.
%
% alpha > -1.

    nq = Q + 1;
    j = (0:Q).';
    diagJ = 2*j + alpha + 1;

    if nq == 1
        J = diagJ;
    else
        i = (1:Q).';
        subJ = sqrt(i .* (i + alpha));
        J = diag(diagJ) + diag(subJ,1) + diag(subJ,-1);
    end

    [V,D] = eig(J);
    y = diag(D);
    [y,idx] = sort(y);
    V = V(:,idx);

    w = gamma(alpha + 1) * (V(1,:).^2).';
end
```

---

### 2.2 Mapped Gauss--LOF quadrature

映射节点为

\[
t_j^{(\alpha,\beta)} = \exp\left(-\frac{y_j^{(\alpha)}}{\beta+1}\right).
\]

LOF 求积权重为

\[
\chi_j^{(\alpha,\beta)} = (\beta+1)^{-\alpha-1}\omega_j^{(\alpha)}.
\]

于是

\[
\int_0^1 p(t)(-\log t)^\alpha t^\beta\,dt = \sum_{j=0}^Q p(t_j)\chi_j^{(\alpha,\beta)},
\qquad p\in\operatorname{span}\{1,\log t,\ldots,(\log t)^{2Q+1}\}.
\]

---

### 2.3 Mapped Gauss--GLOF quadrature

GLOF 求积权重为

\[
\chi_j^{(\alpha,\beta,\lambda)} = (t_j^{(\alpha,\beta)})^{\lambda-\beta}
\chi_j^{(\alpha,\beta)}.
\]

于是

\[
\int_0^1 F(t)(-\log t)^\alpha t^\lambda\,dt
\approx
\sum_{j=0}^Q F(t_j)\chi_j^{(\alpha,\beta,\lambda)}.
\]

在计算投影系数时，若直接使用

\[
F(t)=f(t)S_n^{(\alpha,\beta,\lambda)}(t),
\]

则

\[
F(t_j)\chi_j^{(\alpha,\beta,\lambda)} = f(t_j)t_j^{(\beta-\lambda)/2}L_n^{(\alpha)}(y_j)
\cdot t_j^{\lambda-\beta}\chi_j^{(\alpha,\beta)}.
\]

合并幂次可得更稳定的形式：

\[
F(t_j)\chi_j^{(\alpha,\beta,\lambda)} = f(t_j)t_j^{(\lambda-\beta)/2}
L_n^{(\alpha)}(y_j)
\chi_j^{(\alpha,\beta)}.
\]

因此投影系数应按下面的矩阵公式计算：

\[
\boxed{
\widehat f_n^{(Q)} = \frac{1}{\gamma_n^{(\alpha,\beta)}}
\sum_{j=0}^Q
f(t_j)t_j^{(\lambda-\beta)/2}
L_n^{(\alpha)}(y_j)
\chi_j^{(\alpha,\beta)}.
}
\]

这比先构造 GLOF 基函数和 GLOF 权重再相乘更稳定。

---

## 3. MATLAB 核心函数

### 3.1 计算 generalized Laguerre polynomial table

使用三项递推：

\[
L_0^{(\alpha)}(y)=1,
\qquad
L_1^{(\alpha)}(y)=-y+\alpha+1,
\]

\[
L_{n+1}^{(\alpha)}(y) = \frac{(2n+\alpha+1-y)L_n^{(\alpha)}(y)-(n+\alpha)L_{n-1}^{(\alpha)}(y)}{n+1}.
\]

MATLAB 函数：

```matlab
function L = genlaguerre_table(m, alpha, y)
% L(j,n+1) = L_n^{alpha}(y_j), n=0,...,m.

    y = y(:);
    nq = length(y);
    L = zeros(nq, m+1);

    L(:,1) = 1;
    if m >= 1
        L(:,2) = -y + alpha + 1;
    end

    for n = 1:m-1
        L(:,n+2) = ((2*n + alpha + 1 - y).*L(:,n+1) ...
                   - (n + alpha).*L(:,n)) / (n + 1);
    end
end
```

---

### 3.2 计算 GLOF 投影系数

```matlab
function c = glof_project_coeffs(fhandle, m, Q, alpha, beta, lambda)
% Compute GLOF projection coefficients c_n, n=0,...,m.
%
% Projection:
% c_n = 1/gamma_n * int_0^1 f(t) S_n^{a,b,lambda}(t)
%       (-log t)^alpha t^lambda dt.
%
% The coefficient integral is computed by Gauss-GLOF quadrature.

    [y, wLag] = laguerre_gauss_gw(Q, alpha);

    logt = -y/(beta + 1);
    t = exp(logt);

    L = genlaguerre_table(m, alpha, y);       % size (Q+1) x (m+1)

    chiLOF = (beta + 1)^(-alpha - 1) * wLag;

    % Stable combined factor for projection:
    % f(t_j) * t_j^((lambda-beta)/2) * chiLOF_j
    scale = exp(((lambda - beta)/2) * logt);
    fvals = fhandle(t);
    rhs_weighted = fvals(:) .* scale(:) .* chiLOF(:);

    raw = L.' * rhs_weighted;                 % raw(n+1) = integral part

    n = (0:m).';
    gamma_n = exp(gammaln(n + alpha + 1) - gammaln(n + 1)) ...
              / (beta + 1)^(alpha + 1);

    c = raw ./ gamma_n;
end
```

---

### 3.3 用 GLOF 系数重构函数

```matlab
function u = glof_eval_series(t, c, alpha, beta, lambda)
% Evaluate sum_{n=0}^m c_n S_n^{alpha,beta,lambda}(t).

    t = t(:);
    m = length(c) - 1;

    % Avoid evaluating exactly at t=0 in this routine.
    y = -(beta + 1) * log(t);
    L = genlaguerre_table(m, alpha, y);

    theta = (beta - lambda)/2;
    u = (t.^theta) .* (L * c(:));
end
```

如果需要在 `t=0` 处取值，且 `theta>0`，可单独设置 GLOF 部分为 0；不要直接调用 `log(0)`。

---

## 4. Test 0：Laguerre--Gauss 求积精确性测试

这是所有后续测试的基础。对 \(p=0,1,\ldots,2Q+1\)，应有

\[
\sum_{j=0}^Q \omega_j y_j^p
\approx
\int_0^\infty y^{p+\alpha}e^{-y}\,dy = \Gamma(p+\alpha+1).
\]

MATLAB 测试：

```matlab
function test_laguerre_gauss_exactness()
    alpha_list = [0, 0.5, 2];
    Q = 30;

    for alpha = alpha_list
        [y,w] = laguerre_gauss_gw(Q, alpha);
        maxerr = 0;
        for p = 0:(2*Q+1)
            approx = sum(w .* y.^p);
            exact = gamma(p + alpha + 1);
            relerr = abs(approx - exact) / max(1, abs(exact));
            maxerr = max(maxerr, relerr);
        end
        fprintf('Laguerre exactness: alpha=%g, Q=%d, max relerr=%e\n', ...
                alpha, Q, maxerr);
    end
end
```

期望结果：

```text
max relerr approximately 1e-13 to 1e-11
```

当 `Q` 很大或 `p` 很大时，`gamma(p+alpha+1)` 会很大，此测试可能受浮点溢出影响。建议先用 `Q=20` 到 `Q=40`。

---

## 5. Test 1：GLOF Gram matrix 正交性测试

该测试直接验证 GLOF 节点、权重、基函数和归一化常数是否一致。

令

\[
\Phi_{j,n}=S_n^{(\alpha,\beta,\lambda)}(t_j),
\qquad j=0,
\ldots,Q,
\quad n=0,\ldots,m.
\]

用 Gauss--GLOF 求积计算 Gram 矩阵：

\[
G_{n\ell}^{(Q)} = \sum_{j=0}^Q
S_n^{(\alpha,\beta,\lambda)}(t_j)
S_\ell^{(\alpha,\beta,\lambda)}(t_j)
\chi_j^{(\alpha,\beta,\lambda)}.
\]

理论上，当 \(Q\ge m\) 时，

\[
G_{n\ell}^{(Q)}=\gamma_n^{(\alpha,\beta)}\delta_{n\ell},
\qquad 0\le n,\ell\le m.
\]

MATLAB 测试：

```matlab
function test_glof_orthogonality()
    alpha = 0;
    beta = 5;
    lambda = 0;
    m = 40;
    Q = m;

    [y,wLag] = laguerre_gauss_gw(Q, alpha);
    logt = -y/(beta+1);
    t = exp(logt);
    L = genlaguerre_table(m, alpha, y);

    theta = (beta - lambda)/2;
    Phi = (t.^theta) .* L;

    chiLOF = (beta+1)^(-alpha-1) * wLag;
    chiGLOF = exp((lambda-beta)*logt) .* chiLOF;

    G = Phi.' * (chiGLOF .* Phi);

    n = (0:m).';
    gamma_n = exp(gammaln(n+alpha+1)-gammaln(n+1)) / (beta+1)^(alpha+1);
    Gexact = diag(gamma_n);

    relFrob = norm(G - Gexact, 'fro') / norm(Gexact, 'fro');
    maxOffDiag = max(max(abs(G - diag(diag(G)))));

    fprintf('GLOF Gram relFrob = %.3e, maxOffDiag = %.3e\n', ...
            relFrob, maxOffDiag);
end
```

期望结果：

```text
relFrob approximately 1e-13 to 1e-11
maxOffDiag approximately 1e-14 to 1e-12
```

如果该测试失败，通常说明以下某一项有误：

1. `t_j = exp(-y_j/(beta+1))` 写错；
2. LOF 权重少乘了 `(beta+1)^(-alpha-1)`；
3. GLOF 权重少乘了 `t_j^(lambda-beta)`；
4. GLOF 基函数的幂次应为 `(beta-lambda)/2`，不是 `(lambda-beta)/2`；
5. `gamma_n` 归一化常数写错。

---

## 6. Test 2：投影一个 GLOF 基函数，检查 Kronecker delta 系数

取

\[
f(t)=S_{n_0}^{(\alpha,\beta,\lambda)}(t).
\]

则理论投影系数为

\[
\widehat f_n=\delta_{n,n_0}.
\]

MATLAB 测试：

```matlab
function test_project_one_basis()
    alpha = 0;
    beta = 5;
    lambda = 0;
    m = 50;
    Q = m;
    n0_list = [0, 1, 5, 20, 40];

    for n0 = n0_list
        fhandle = @(t) glof_single_basis(t, n0, alpha, beta, lambda);
        c = glof_project_coeffs(fhandle, m, Q, alpha, beta, lambda);

        cExact = zeros(m+1,1);
        cExact(n0+1) = 1;

        err = norm(c - cExact, inf);
        fprintf('Basis projection: n0=%d, coeff inf err=%e\n', n0, err);
    end
end

function val = glof_single_basis(t, n0, alpha, beta, lambda)
    t = t(:);
    y = -(beta + 1) * log(t);
    L = genlaguerre_table(n0, alpha, y);
    theta = (beta - lambda)/2;
    val = t.^theta .* L(:,n0+1);
end
```

期望结果：

```text
coeff inf err approximately 1e-13 to 1e-11
```

这个测试是最直接的投影系数代码正确性测试。它不依赖解析公式，也不依赖函数重构图像。

---

## 7. Test 3：对 \(f(t)=t^r(-\log t)^k\) 的解析系数验证

这是与研究问题最相关的测试。

### 7.1 解析系数公式

考虑

\[
f(t)=t^r(-\log t)^k,
\qquad r\ge0,
\quad k\in\mathbb N_0.
\]

GLOF 系数为

\[
\widehat f_n = \frac{1}{\gamma_n^{(\alpha,\beta)}}
\int_0^1
f(t)S_n^{(\alpha,\beta,\lambda)}(t)
(-\log t)^\alpha t^\lambda\,dt.
\]

令

\[
s=
\frac{\beta+\lambda+2r+2}{2(\beta+1)}.
\]

则

\[
\widehat f_n = (\beta+1)^{-k}
\frac{\Gamma(n+1)}{\Gamma(n+\alpha+1)}
I_{n,k},
\]

其中

\[
I_{n,k} = \int_0^\infty e^{-sy}y^{\alpha+k}L_n^{(\alpha)}(y)\,dy.
\]

为了在 MATLAB 中对任意小阶数稳定验证，可用 Laguerre 多项式展开：

\[
L_n^{(\alpha)}(y) = \sum_{p=0}^n
\frac{(-1)^p}{p!}
\binom{n+\alpha}{n-p}y^p.
\]

于是

\[
I_{n,k} = \sum_{p=0}^n
\frac{(-1)^p}{p!}
\binom{n+\alpha}{n-p}
\frac{\Gamma(\alpha+k+p+1)}{s^{\alpha+k+p+1}}.
\]

MATLAB 参考实现：

```matlab
function cExact = glof_exact_coeff_tlog(m, alpha, beta, lambda, r, k)
% Exact GLOF coefficients for f(t)=t^r*(-log t)^k.
% Suitable for validation with moderate m.

    s = (beta + lambda + 2*r + 2) / (2*(beta + 1));
    cExact = zeros(m+1,1);

    for n = 0:m
        I = 0;
        for p = 0:n
            % binom(n+alpha, n-p) = Gamma(n+alpha+1) /
            %                         [Gamma(n-p+1) Gamma(alpha+p+1)]
            logBinom = gammaln(n + alpha + 1) ...
                     - gammaln(n - p + 1) ...
                     - gammaln(alpha + p + 1);

            term = (-1)^p * exp(logBinom - gammaln(p + 1) ...
                   + gammaln(alpha + k + p + 1) ...
                   - (alpha + k + p + 1)*log(s));

            I = I + term;
        end

        cExact(n+1) = (beta + 1)^(-k) ...
                    * exp(gammaln(n + 1) - gammaln(n + alpha + 1)) ...
                    * I;
    end
end
```

> 注意：该展开对较大 `m` 可能出现交错求和的精度损失。因此它主要用于 `m <= 60` 左右的单元测试。更大规模时可使用 `vpa` 或只做 `Q`-doubling 检查。

---

### 7.2 调谐参数 sanity check

若选择

\[
\frac{\beta-\lambda}{2}=r,
\]

则

\[
S_n^{(\alpha,\beta,\lambda)}(t) = t^rL_n^{(\alpha)}(y).
\]

因此

\[
f(t)=t^r(-\log t)^k
=t^r\left(\frac{y}{\beta+1}\right)^k
\in
\operatorname{span}\{S_0^{(\alpha,\beta,\lambda)},\ldots,S_k^{(\alpha,\beta,\lambda)}\}.
\]

所以当 \(m\ge k\) 时，理论上只需要前 \(k+1\) 个 GLOF 模式即可精确表示。

MATLAB 测试：

```matlab
function test_tuned_tlog_exact_recovery()
    alpha = 0;
    r = 0.3;
    k = 3;

    lambda = 0;
    beta = 2*r;             % ensures (beta-lambda)/2 = r

    m = 12;
    Q = 30;

    fhandle = @(t) t.^r .* (-log(t)).^k;

    cNum = glof_project_coeffs(fhandle, m, Q, alpha, beta, lambda);
    cExact = glof_exact_coeff_tlog(m, alpha, beta, lambda, r, k);

    coeffErr = norm(cNum - cExact, inf);
    tailErr = norm(cNum((k+2):end), inf);

    fprintf('Tuned t^r log^k: coeffErr=%e, tailErr n>k=%e\n', ...
            coeffErr, tailErr);

    % Reconstruction test on a dense grid excluding t=0.
    tt = linspace(1e-12, 1, 2000).';
    ftrue = fhandle(tt);
    frec = glof_eval_series(tt, cNum, alpha, beta, lambda);
    recErr = norm(frec - ftrue, inf);
    fprintf('Tuned reconstruction Linf err=%e\n', recErr);
end
```

期望结果：

```text
coeffErr approximately 1e-12 to 1e-10
tailErr n>k approximately 1e-12 to 1e-10
reconstruction Linf err approximately 1e-11 to 1e-9
```

若该测试失败，优先检查 `t^((lambda-beta)/2)` 的合并权重是否写错。

---

## 8. Test 4：复现 Chen--Shen 的 Gauss--LOF quadrature 测试

Chen--Shen 文章的第一个数值验证是 Gauss--LOF quadrature 对如下函数的积分误差呈指数衰减：

\[
f(t)=\sin t,
\qquad f(t)=e^t,
\qquad f(t)=t^{-1/3},
\qquad f(t)=t^{1/10}.
\]

先取最简单的权重

\[
\alpha=0,
\qquad \beta=0.
\]

此时积分为普通积分

\[
\int_0^1 f(t)\,dt.
\]

精确值为

\[
\int_0^1 \sin t\,dt=1-\cos1,
\]

\[
\int_0^1 e^t\,dt=e-1,
\]

\[
\int_0^1 t^{-1/3}\,dt=\frac32,
\]

\[
\int_0^1 t^{1/10}\,dt=\frac{10}{11}.
\]

MATLAB 测试：

```matlab
function test_gauss_lof_quadrature_fig3_style()
    alpha = 0;
    beta = 0;

    funcs = {
        @(t) sin(t),      1 - cos(1),       'sin(t)';
        @(t) exp(t),      exp(1) - 1,       'exp(t)';
        @(t) t.^(-1/3),   3/2,              't^{-1/3}';
        @(t) t.^(1/10),   10/11,            't^{1/10}'
    };

    Qlist = 2:2:50;
    figure; hold on;

    for id = 1:size(funcs,1)
        fhandle = funcs{id,1};
        exact = funcs{id,2};
        label = funcs{id,3};
        errs = zeros(size(Qlist));

        for iq = 1:length(Qlist)
            Q = Qlist(iq);
            [y,wLag] = laguerre_gauss_gw(Q, alpha);
            t = exp(-y/(beta+1));
            chiLOF = (beta+1)^(-alpha-1) * wLag;
            approx = sum(fhandle(t) .* chiLOF);
            errs(iq) = abs(approx - exact);
        end

        semilogy(Qlist, errs, '-o', 'DisplayName', label);
    end

    xlabel('Q'); ylabel('absolute quadrature error');
    legend('Location','southwest'); grid on;
    title('Gauss--LOF quadrature validation');
end
```

期望现象：

1. 四条曲线均应快速下降；
2. 即使 \(t^{-1/3}\) 在 \(t=0\) 奇异，也应观察到快速收敛；
3. 误差最终到达双精度平台，大约 \(10^{-13}\) 到 \(10^{-15}\)。

该测试通过后，说明 mapped Laguerre quadrature 的节点和权重基本正确。

---

## 9. Test 5：复现 Chen--Shen 的 GLOF approximation 测试

Chen--Shen 文章中使用 GLOF 验证了对弱奇异函数的指数逼近。典型参数为

\[
\alpha=0,
\qquad \beta=5,
\qquad \lambda=0.
\]

此时 GLOF 为

\[
S_n^{(0,5,0)}(t)=t^{5/2}L_n^{(0)}(-6\log t).
\]

### 9.1 测试 \(f(t)=t^r\)

取

\[
r=0.1,0.3,0.5,0.7,0.9
\]

以及整数指数

\[
r=1,2,3,4,5.
\]

误差使用加权 \(L^2_{\chi_{\alpha,\lambda}}\) 范数；当 \(\alpha=\lambda=0\) 时，就是普通 \(L^2(0,1)\) 范数。

精确误差可由 Parseval 关系计算：

\[
\|f-P_m f\|_{\chi_{\alpha,\lambda}}^2 = \|f\|_{\chi_{\alpha,\lambda}}^2
-
\sum_{n=0}^m \gamma_n^{(\alpha,\beta)}|\widehat f_n|^2.
\]

对于

\[
f(t)=t^r(-\log t)^k,
\]

有

\[
\|f\|_{\chi_{\alpha,\lambda}}^2 = \int_0^1 t^{2r+\lambda}(-\log t)^{\alpha+2k}\,dt = \frac{\Gamma(\alpha+2k+1)}{(2r+\lambda+1)^{\alpha+2k+1}}.
\]

MATLAB 函数：

```matlab
function err = glof_projection_error_from_coeff(c, alpha, beta, lambda, r, k)
% Weighted L2 projection error for f=t^r*(-log t)^k,
% computed from coefficients.

    m = length(c) - 1;
    n = (0:m).';
    gamma_n = exp(gammaln(n+alpha+1)-gammaln(n+1)) / (beta+1)^(alpha+1);

    normsq = gamma(alpha + 2*k + 1) / (2*r + lambda + 1)^(alpha + 2*k + 1);
    projnormsq = sum(gamma_n .* abs(c(:)).^2);

    errsq = max(normsq - projnormsq, 0);
    err = sqrt(errsq);
end
```

用数值投影系数复现实验：

```matlab
function test_glof_approx_tr_fig5_style()
    alpha = 0;
    beta = 5;
    lambda = 0;
    k = 0;

    rlist = [0.1, 0.3, 0.5, 0.7, 0.9];
    mlist = 0:5:70;

    figure; hold on;

    for r = rlist
        errs = zeros(size(mlist));
        for ii = 1:length(mlist)
            m = mlist(ii);
            Q = max(80, 3*m + 20);
            fhandle = @(t) t.^r;
            c = glof_project_coeffs(fhandle, m, Q, alpha, beta, lambda);
            errs(ii) = glof_projection_error_from_coeff(c, alpha, beta, lambda, r, k);
        end
        semilogy(mlist, errs, '-o', 'DisplayName', sprintf('r=%g', r));
    end

    xlabel('m'); ylabel('weighted L2 projection error');
    legend('Location','southwest'); grid on;
    title('GLOF projection: f(t)=t^r, alpha=lambda=0, beta=5');
end
```

期望现象：

1. 每条曲线应近似直线下降，即指数收敛；
2. \(r\) 越靠近与 \(\beta\) 的最佳适配范围，收敛越快；
3. 最后会进入双精度平台。

---

### 9.2 测试 \(f(t)=t(-\log t)^k\) 和 \(f(t)=t^2(-\log t)^k\)

取

\[
k=1,2,3,4,5.
\]

MATLAB 测试：

```matlab
function test_glof_approx_tlog_fig6_style()
    alpha = 0;
    beta = 5;
    lambda = 0;

    rlist = [1, 2];
    klist = 1:5;
    mlist = 0:5:60;

    for r = rlist
        figure; hold on;
        for k = klist
            errs = zeros(size(mlist));
            for ii = 1:length(mlist)
                m = mlist(ii);
                Q = max(100, 3*m + 30);
                fhandle = @(t) t.^r .* (-log(t)).^k;
                c = glof_project_coeffs(fhandle, m, Q, alpha, beta, lambda);
                errs(ii) = glof_projection_error_from_coeff(c, alpha, beta, lambda, r, k);
            end
            semilogy(mlist, errs, '-o', 'DisplayName', sprintf('k=%d', k));
        end

        xlabel('m'); ylabel('weighted L2 projection error');
        legend('Location','southwest'); grid on;
        title(sprintf('GLOF projection: f(t)=t^{%g}(-log t)^k', r));
    end
end
```

期望现象：

1. 对所有 \(k\)，误差应指数下降；
2. \(k\) 越大，前期误差可能越大，但仍应保持指数趋势；
3. 这组测试最接近你的研究模型 \(t^r(-\log t)^k\)。

---

## 10. Test 6：数值系数与解析系数逐项比较

该测试直接检查 `glof_project_coeffs` 是否正确。

```matlab
function test_coeff_against_exact_formula()
    alpha = 0;
    beta = 5;
    lambda = 0;

    cases = [
        0.1, 0;
        0.3, 0;
        1.0, 1;
        1.0, 3;
        2.0, 2
    ];

    m = 40;
    Q = 160;

    for icase = 1:size(cases,1)
        r = cases(icase,1);
        k = cases(icase,2);

        fhandle = @(t) t.^r .* (-log(t)).^k;
        cNum = glof_project_coeffs(fhandle, m, Q, alpha, beta, lambda);
        cExact = glof_exact_coeff_tlog(m, alpha, beta, lambda, r, k);

        absErr = norm(cNum - cExact, inf);
        relErr = absErr / max(1, norm(cExact, inf));

        fprintf('r=%g, k=%d: coeff absErr=%e, relErr=%e\n', ...
                r, k, absErr, relErr);
    end
end
```

期望结果：

```text
coeff relErr roughly 1e-11 to 1e-9 for moderate m
```

若误差明显更大，先做以下检查：

1. 将 `Q` 增大到 `2Q`，检查是否是求积阶数不足；
2. 用 Test 1 检查 Gram matrix；
3. 用 Test 2 检查单个基函数投影；
4. 检查 `glof_exact_coeff_tlog` 是否因高阶交错求和失去精度。必要时用 `vpa` 或只测试较小 `m`。

---

## 11. Q-doubling 检查：避免把求积误差误认为投影误差

在实际重构实验中，投影系数来自 Fourier 近似或 Chebyshev 配点近似，不一定有解析系数。此时必须做 `Q`-doubling 检查。

对于同一 \(m,\alpha,\beta,\lambda\)，计算

\[
\widehat f_n^{(Q)},
\qquad
\widehat f_n^{(2Q)}.
\]

定义

\[
E_Q^{\rm coeff} = \frac{
\left(\sum_{n=0}^m \gamma_n|\widehat f_n^{(2Q)}-\widehat f_n^{(Q)}|^2\right)^{1/2}
}{
\left(\sum_{n=0}^m \gamma_n|\widehat f_n^{(2Q)}|^2\right)^{1/2}
}.
\]

MATLAB 函数：

```matlab
function errQ = coeff_q_doubling_error(fhandle, m, Q, alpha, beta, lambda)
    cQ = glof_project_coeffs(fhandle, m, Q, alpha, beta, lambda);
    c2Q = glof_project_coeffs(fhandle, m, 2*Q, alpha, beta, lambda);

    n = (0:m).';
    gamma_n = exp(gammaln(n+alpha+1)-gammaln(n+1)) / (beta+1)^(alpha+1);

    num = sqrt(sum(gamma_n .* abs(c2Q - cQ).^2));
    den = sqrt(sum(gamma_n .* abs(c2Q).^2));
    errQ = num / max(den, eps);
end
```

建议规则：

```matlab
if errQ < 1e-12
    % coefficient quadrature is sufficiently resolved
else
    % increase Q, or inspect beta/node distribution/stability
end
```

对于后处理实验，通常取

\[
Q=2\max\{N_F,m\}
\]

作为起点，再用 `Q`-doubling 检查。

---

## 12. 节点过度靠近零时的诊断

GLOF / LOF 节点为

\[
t_j=e^{-y_j/(\beta+1)}.
\]

当 \(Q\) 很大而 \(\beta\) 固定时，最大 Laguerre 节点 \(y_{\max}\) 较大，导致

\[
t_{\min}=e^{-y_{\max}/(\beta+1)}
\]

可能极小。此时若使用调谐参数 \((\beta-\lambda)/2=r\)，投影系数公式中会出现

\[
t_j^{(\lambda-\beta)/2}=t_j^{-r}.
\]

对精确函数 \(f(t)=t^r(-\log t)^k\) 来说，乘积稳定；但对 Fourier partial sum 或 Chebyshev interpolant 来说，端点附近的近似函数未必严格满足 \(t^r\) 行为，可能放大截断误差。

建议在每次实验中打印：

```matlab
fprintf('Q=%d, min(t)=%e, max(t)=%e\n', Q, min(t), max(t));
```

若 `min(t)` 小于 `1e-300`，可能已经接近下溢；此时不应盲目增大 `Q`，而应考虑：

1. 增大 \(\beta\)，使节点远离零；
2. 使用随 \(N\) 变化的参数，例如 \(\beta_N=2r+bN\)，\(\lambda_N=bN\)；
3. 使用 `Q`-doubling 检查，而不是无限制增加 `Q`。

---

## 13. 推荐的完整测试顺序

建议按以下顺序验证 MATLAB 程序：

### Step 1：验证 Laguerre--Gauss 节点和权重

运行：

```matlab
test_laguerre_gauss_exactness();
```

若失败，不要继续测试 GLOF。

---

### Step 2：验证 GLOF 正交性

运行：

```matlab
test_glof_orthogonality();
```

若失败，重点检查映射、权重和归一化常数。

---

### Step 3：验证单个基函数投影

运行：

```matlab
test_project_one_basis();
```

若失败，重点检查投影系数公式中的合并因子

\[
t^{(\lambda-\beta)/2}.
\]

---

### Step 4：验证调谐情形精确恢复

运行：

```matlab
test_tuned_tlog_exact_recovery();
```

这是最接近你的模型 \(t^r(-\log t)^k\) 的精确性测试。

---

### Step 5：复现 Chen--Shen Gauss--LOF quadrature 图像

运行：

```matlab
test_gauss_lof_quadrature_fig3_style();
```

应看到所有曲线快速下降。

---

### Step 6：复现 Chen--Shen GLOF approximation 图像

运行：

```matlab
test_glof_approx_tr_fig5_style();
test_glof_approx_tlog_fig6_style();
```

应看到加权 \(L^2\) 投影误差随 \(m\) 指数下降。

---

### Step 7：逐项比较解析系数

运行：

```matlab
test_coeff_against_exact_formula();
```

这是确认投影系数代码正确性的最终单元测试。

---

## 14. 常见错误清单

### 错误 1：把 GLOF 权重写成 LOF 权重

错误：

\[
\chi_j=\chi_j^{(\alpha,\beta)}.
\]

正确：

\[
\chi_j^{(\alpha,\beta,\lambda)}
=t_j^{\lambda-\beta}\chi_j^{(\alpha,\beta)}.
\]

或者在投影系数中使用合并形式：

\[
f(t_j)t_j^{(\lambda-\beta)/2}L_n^{(\alpha)}(y_j)
\chi_j^{(\alpha,\beta)}.
\]

---

### 错误 2：GLOF 基函数幂次符号写反

正确：

\[
S_n^{(\alpha,\beta,\lambda)}(t)
=t^{(\beta-\lambda)/2}L_n^{(\alpha)}(y).
\]

不是

\[
t^{(\lambda-\beta)/2}L_n^{(\alpha)}(y).
\]

注意：投影系数的合并形式里出现的是 \(t^{(\lambda-\beta)/2}\)，这是因为基函数和 GLOF 权重相乘后合并得到的，不是基函数本身。

---

### 错误 3：忘记 \((\beta+1)^{-\alpha-1}\)

LOF 权重必须包含

\[
\chi_j^{(\alpha,\beta)}=(\beta+1)^{-\alpha-1}\omega_j^{(\alpha)}.
\]

---

### 错误 4：归一化常数 \(\gamma_n\) 写错

正确：

\[
\gamma_n^{(\alpha,\beta)} = \frac{\Gamma(n+\alpha+1)}{(\beta+1)^{\alpha+1}\Gamma(n+1)}.
\]

不是 Laguerre 原始范数

\[
\frac{\Gamma(n+\alpha+1)}{\Gamma(n+1)}.
\]

二者相差 \((\beta+1)^{-\alpha-1}\)。

---

### 错误 5：在 `t=0` 处直接计算 GLOF

GLOF 节点不会等于 0，但重构图像的测试网格如果包含 `t=0`，直接 `log(0)` 会出错。应单独处理端点值。

---

## 15. 最小主程序模板

下面是一个建议的主程序结构：

```matlab
clear; clc; close all;

fprintf('Step 1: Laguerre-Gauss exactness\n');
test_laguerre_gauss_exactness();

fprintf('\nStep 2: GLOF orthogonality\n');
test_glof_orthogonality();

fprintf('\nStep 3: Project one GLOF basis\n');
test_project_one_basis();

fprintf('\nStep 4: Tuned t^r(-log t)^k exact recovery\n');
test_tuned_tlog_exact_recovery();

fprintf('\nStep 5: Chen-Shen Gauss-LOF quadrature style test\n');
test_gauss_lof_quadrature_fig3_style();

fprintf('\nStep 6: Chen-Shen GLOF approximation style tests\n');
test_glof_approx_tr_fig5_style();
test_glof_approx_tlog_fig6_style();

fprintf('\nStep 7: Coefficients against exact formula\n');
test_coeff_against_exact_formula();
```

---

## 16. 通过标准

在双精度 MATLAB 中，合理的通过标准为：

| 测试 | 期望误差 |
|---|---:|
| Laguerre--Gauss monomial exactness | \(10^{-13}\sim10^{-11}\) |
| GLOF Gram matrix | \(10^{-13}\sim10^{-11}\) |
| Projection of one basis function | \(10^{-13}\sim10^{-11}\) |
| Tuned \(t^r(-\log t)^k\) tail coefficients | \(10^{-12}\sim10^{-10}\) |
| Coefficients vs analytic formula | \(10^{-11}\sim10^{-9}\) |
| Chen--Shen quadrature curves | exponential decay to machine precision |
| Chen--Shen GLOF approximation curves | exponential decay to machine precision |

若以上测试全部通过，则可以认为 GLOF 投影系数计算模块基本正确，可以进入下一步：

1. 用 Fourier partial sum \(F_N\) 代替精确 \(f\)，测试从 Fourier 系数计算 GLOF 系数；
2. 用 Chebyshev interpolant \(I_N f\) 代替精确 \(f\)，测试从配点值计算 GLOF 系数；
3. 比较重构误差、传统 Fourier/Chebyshev 误差和 `Q`-doubling 误差。

