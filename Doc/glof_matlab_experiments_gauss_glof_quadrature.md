# MATLAB 数值实验指南：用 Chen--Shen 的 Gauss--GLOF 求积计算 GLOF 重构系数

本文件是对前一版数值实验方案的关键修正：**从 Fourier 系数或配点值计算近似 GLOF 系数时，不建议逐项调用 MATLAB `integral`，也不建议直接使用普通 Gauss--Legendre 求积。应优先使用 Chen--Shen 文章中的 Gauss--LOF / Gauss--GLOF 求积节点与权重。**

目标函数先取

$$
f(x)=x^r(-\log x)^k,\qquad x\in(0,1],\qquad r>0,\qquad k\in\mathbb N.
$$

约定

$$
f(0)=0.
$$

---

## 1. 实验总目标

需要比较三类近似。

### 1.1 传统 Fourier 表示

在 $[0,1]$ 上使用周期 Fourier 部分和

$$
f_N^{\rm F}(x)=\sum_{|\ell|\le N}c_\ell e^{2\pi i\ell x},
\qquad
c_\ell=\int_0^1 f(x)e^{-2\pi i\ell x}\,dx.
$$

由于纯模型函数满足 $f(0)=f(1)=0$，周期延拓在函数值层面没有跳跃，所以未必出现典型的跳跃 Gibbs overshoot；更应观察的是端点奇异导致的 Gibbs-type 振荡、端点误差聚集以及谱精度退化。

### 1.2 传统 Chebyshev 配点插值

使用 $[0,1]$ 上的 Chebyshev--Lobatto 点

$$
x_j=\frac{1+\cos(j\pi/N)}2,\qquad j=0,1,\ldots,N.
$$

由于 $x=0$ 处导数奇异，Chebyshev 插值误差通常不会呈指数衰减。

### 1.3 GLOF 后处理重构

从以下两类标准谱数据出发：

1. Fourier 系数 $\{c_\ell\}_{|\ell|\le N}$；
2. Chebyshev 配点值 $\{f(x_j)\}_{j=0}^N$；

先构造标准谱近似 $F_N$，然后计算近似 GLOF 系数

$$
g_n^{(N)}
\approx
\frac{1}{\gamma_n^{(\alpha,\beta)}}
\int_0^1
F_N(x)S_n^{(\alpha,\beta,\lambda)}(x)
(-\log x)^\alpha x^\lambda\,dx,
\qquad n=0,\ldots,m.
$$

最后重构

$$
R_{N,m}^{\rm GLOF}f(x)=\sum_{n=0}^m g_n^{(N)}S_n^{(\alpha,\beta,\lambda)}(x).
$$

本文件重点说明上述积分如何用 Gauss--GLOF 求积高效计算。

---

## 2. GLOF 定义与正交性

令

$$
y=-(\beta+1)\log x,
\qquad
x=e^{-y/(\beta+1)},
\qquad y\in[0,\infty).
$$

LOF 为

$$
S_n^{(\alpha,\beta)}(x)=L_n^{(\alpha)}(y),
$$

GLOF 为

$$
S_n^{(\alpha,\beta,\lambda)}(x)= x^{(\beta-\lambda)/2}L_n^{(\alpha)}(-(\beta+1)\log x).
$$

记

$$
\theta=\frac{\beta-\lambda}{2}.
$$

则

$$
S_n^{(\alpha,\beta,\lambda)}(x)=x^\theta L_n^{(\alpha)}(y).
$$

GLOF 的正交性为

$$
\int_0^1
S_n^{(\alpha,\beta,\lambda)}(x)
S_m^{(\alpha,\beta,\lambda)}(x)
(-\log x)^\alpha x^\lambda\,dx=\gamma_n^{(\alpha,\beta)}\delta_{mn},
$$

其中

$$
\gamma_n^{(\alpha,\beta)}=
\frac{\Gamma(n+\alpha+1)}{(\beta+1)^{\alpha+1}\Gamma(n+1)}.
$$

---

## 3. 参数选择

对

$$
f(x)=x^r(-\log x)^k
$$

最自然的调谐条件是

$$
\theta=\frac{\beta-\lambda}{2}=r.
$$

也就是

$$
\lambda=\beta-2r.
$$

这样

$$
S_n^{(\alpha,\beta,\lambda)}(x) = x^rL_n^{(\alpha)}(y),
$$

而 $(-
\log x)^k=(y/(\beta+1))^k$ 是 $y$ 的 $k$ 次多项式。因此精确 GLOF 投影在 $m\ge k$ 时可以表示纯模型函数。

### 3.1 固定参数 sanity check

先取

$$
\alpha=0,
\qquad
\beta=2r,
\qquad
\lambda=0.
$$

这是最简单的调谐选择。

### 3.2 更接近 Chen--Shu 重构思想的 $N$-dependent 参数

为了从截断 Fourier / Chebyshev 数据中恢复精度，更推荐测试

$$
\alpha_N=aN,
\qquad
\beta_N=2r+bN,
\qquad
\lambda_N=bN,
\qquad
m=\max\{k,\lfloor \gamma N\rfloor\}.
$$

此时仍有

$$
\frac{\beta_N-\lambda_N}{2}=r.
$$

建议初值：

$$
a=0.03\sim0.08,
\qquad
b=0.03\sim0.10,
\qquad
\gamma=0.03\sim0.10.
$$

---

## 4. Chen--Shen Gauss--LOF / Gauss--GLOF 求积

### 4.1 Laguerre--Gauss 节点

令 $\{y_q^{(\alpha)},\omega_q^{(\alpha)}\}_{q=0}^Q$ 为 $(Q+1)$ 点 Laguerre--Gauss 求积节点和权重，使

$$
\int_0^\infty p(y)y^\alpha e^{-y}\,dy = \sum_{q=0}^Q \omega_q^{(\alpha)}p(y_q^{(\alpha)})
$$

对 $\deg p\le 2Q+1$ 精确。

### 4.2 LOF 映射节点和权重

由 Chen--Shen 的 log mapping，定义

$$
x_q=e^{-y_q^{(\alpha)}/(\beta+1)},
$$

以及 LOF 权重

$$
\chi_q^{(\alpha,\beta)} = (\beta+1)^{-\alpha-1}\omega_q^{(\alpha)}.
$$

于是

$$
\int_0^1 p(x)(-\log x)^\alpha x^\beta\,dx
\approx
\sum_{q=0}^Q p(x_q)\chi_q^{(\alpha,\beta)}.
$$

当 $p$ 是 $\log x$ 的次数不超过 $2Q+1$ 的多项式时，上式精确。

### 4.3 GLOF 权重

GLOF 求积权重为

$$
\chi_q^{(\alpha,\beta,\lambda)} = x_q^{\lambda-\beta}\chi_q^{(\alpha,\beta)} = x_q^{\lambda-\beta}(\beta+1)^{-\alpha-1}\omega_q^{(\alpha)}.
$$

于是

$$
\int_0^1 h(x)(-\log x)^\alpha x^\lambda\,dx
\approx
\sum_{q=0}^Q h(x_q)\chi_q^{(\alpha,\beta,\lambda)}.
$$

在理论上，该公式对

$$
h(x)\in x^{\beta-\lambda}\operatorname{span}\{1,\log x,\ldots,(\log x)^{2Q+1}\}
$$

精确。

---

## 5. GLOF 系数积分的离散公式

对任意由标准谱数据构造出来的近似函数 $F_N$，定义近似 GLOF 系数

$$
g_n^{(N,Q)} = \frac{1}{\gamma_n^{(\alpha,\beta)}} \sum_{q=0}^Q F_N(x_q) S_n^{(\alpha,\beta,\lambda)}(x_q) \chi_q^{(\alpha,\beta,\lambda)},
\qquad n=0,\ldots,m.
$$

这就是程序中的核心公式。不要在这里逐项调用 `integral`。

由于

$$
S_n^{(\alpha,\beta,\lambda)}(x_q)=x_q^\theta L_n^{(\alpha)}(y_q),
\qquad
\theta=\frac{\beta-\lambda}{2},
$$

也可以写成

$$
g_n^{(N,Q)} = \frac{1}{\gamma_n^{(\alpha,\beta)}} \sum_{q=0}^Q F_N(x_q) x_q^\theta L_n^{(\alpha)}(y_q) x_q^{\lambda-\beta}(\beta+1)^{-\alpha-1}\omega_q^{(\alpha)}.
$$

即

$$
g_n^{(N,Q)} = \frac{(\beta+1)^{-\alpha-1}}{\gamma_n^{(\alpha,\beta)}} \sum_{q=0}^Q F_N(x_q) L_n^{(\alpha)}(y_q) x_q^{-\theta} \omega_q^{(\alpha)}.
$$

最后一个形式在数学上等价，但在数值上可能出现 $x_q^{-\theta}$ 较大的情况。程序里建议使用矩阵形式，并检查求积阶数 $Q$ 的收敛性。

---

## 6. 矩阵形式

令

$$
F_q=F_N(x_q),
$$

$$
S_{nq}=S_n^{(\alpha,\beta,\lambda)}(x_q),
\qquad n=0,\ldots,m,
\qquad q=0,\ldots,Q,
$$

$$
W_q=\chi_q^{(\alpha,\beta,\lambda)},
$$

$$
D_\gamma=\operatorname{diag}(\gamma_0^{(\alpha,\beta)},\ldots,\gamma_m^{(\alpha,\beta)}).
$$

则

$$
\mathbf g =
D_\gamma^{-1}S\operatorname{diag}(W)\mathbf F.
$$

MATLAB 中直接写成

```matlab = ```

---

## 7. MATLAB 基础函数

### 7.1 模型函数

```matlab
function val = f_logsing(x,r,k)
    val = zeros(size(x));
    idx = (x > 0);
    val(idx) = x(idx).^r .* (-log(x(idx))).^k;
    val(~idx) = 0;
end
```

### 7.2 Laguerre--Gauss 节点和权重

使用 Golub--Welsch 算法。

```matlab
function [y,w] = gauss_laguerre(Q,alpha)
    % Q means degree Q: returns Q+1 nodes.
    n = Q + 1;
    j = (0:n-1)';
    a = 2*j + alpha + 1;
    b = sqrt((1:n-1)'.*((1:n-1)' + alpha));
    J = diag(a) + diag(b,1) + diag(b,-1);
    [V,D] = eig(J);
    y = diag(D);
    [y,idx] = sort(y);
    V = V(:,idx);
    w = gamma(alpha+1)*(V(1,:)').^2;
end
```

### 7.3 Laguerre 多项式矩阵

返回

$$
L_{nq}=L_n^{(\alpha)}(y_q),
\qquad n=0,\ldots,m.
$$

```matlab
function L = laguerre_matrix(m,alpha,y)
    y = y(:).';
    nq = length(y);
    L = zeros(m+1,nq);
    L(1,:) = 1;
    if m >= 1
        L(2,:) = -y + alpha + 1;
    end
    for n = 1:m-1
        L(n+2,:) = ((2*n + alpha + 1 - y)/(n+1)).*L(n+1,:) ...
               - ((n + alpha)/(n+1))*L(n,:);
    end
end
```

### 7.4 GLOF 求积节点、权重、基函数矩阵

```matlab
function [xq,Wg,S,L,y,wLag,gamma_vec] = glof_quadrature_matrices(Q,m,alpha,beta,lambda)
    theta = (beta-lambda)/2;

    % Laguerre--Gauss nodes and weights for y^alpha exp(-y).
    [y,wLag] = gauss_laguerre(Q,alpha);

    % Chen--Shen mapped nodes.
    xq = exp(-y/(beta+1));

    % GLOF quadrature weights for (-log x)^alpha x^lambda dx.
    Wg = (xq.^(lambda-beta)) .* (beta+1)^(-alpha-1) .* wLag;

    % Laguerre values and GLOF basis values.
    L = laguerre_matrix(m,alpha,y);
    S = (xq(:)'.^theta) .* L;

    % Norm constants gamma_n.
    gamma_vec = zeros(m+1,1);
    for n = 0:m
        gamma_vec(n+1) = exp(gammaln(n+alpha+1) - gammaln(n+1) ...
                         - (alpha+1)*log(beta+1));
    end
end
```

---

## 8. 从 Fourier 系数计算 GLOF 系数

### 8.1 已知 Fourier 系数

假设已知

$$
\{c_\ell\}_{|\ell|\le N}.
$$

构造

$$
F_N(x)=f_N^{\rm F}(x)=\sum_{|\ell|\le N}c_\ell e^{2\pi i\ell x}.
$$

在 GLOF 求积节点上评价：

$$
F_q=\sum_{|\ell|\le N}c_\ell e^{2\pi i\ell x_q}.
$$

### 8.2 MATLAB 函数

```matlab
function Fq = eval_fourier_on_nodes(xq,ell,c)
    % xq: column vector, ell and c are row or column vectors.
    xq = xq(:);
    ell = ell(:).';
    c = c(:).';
    E = exp(2*pi*1i*(xq*ell));
    Fq = real(E*c.');
end
```

```matlab
function g = glof_coeffs_from_fourier_gaussglof(ell,c,m,alpha,beta,lambda,Q)
    [xq,Wg,S,~,~,~,gamma_vec] = glof_quadrature_matrices(Q,m,alpha,beta,lambda);
    Fq = eval_fourier_on_nodes(xq,ell,c);
    g = (S * (Wg(:).*Fq(:))) ./ gamma_vec(:);
    g = real(g);
end
```

### 8.3 Fourier 系数计算

为了把误差来源集中在重构上，第一阶段可以用高精度数值积分计算 Fourier 系数：

```matlab
function [ell,c] = fourier_coeffs_quad(N,r,k)
    ell = -N:N;
    c = zeros(size(ell));
    for j = 1:length(ell)
        l = ell(j);
        integrand = @(x) f_logsing(x,r,k).*exp(-2*pi*1i*l*x);
        c(j) = integral(integrand,0,1, ...
                        'RelTol',1e-12,'AbsTol',1e-14, ...
                        'ArrayValued',true);
    end
end
```

注意：这里的 `integral` 只用于生成参考 Fourier 系数。**GLOF 系数积分本身不要用 `integral`，而用 Gauss--GLOF 求积。**

---

## 9. 从 Chebyshev 配点值计算 GLOF 系数

### 9.1 Chebyshev--Lobatto 节点

```matlab
function x = cheb_lobatto_nodes(N)
    j = 0:N;
    x = (1 + cos(pi*j/N))/2;
    x = fliplr(x);  % increasing order from 0 to 1
end
```

### 9.2 Barycentric 权重和评价

```matlab
function w = bary_weights_cheb_lobatto(N)
    w = ones(1,N+1);
    w(1) = 0.5;
    w(end) = 0.5;
    w = w .* (-1).^(0:N);
    w = fliplr(w);
end
```

```matlab
function p = bary_eval(xeval,xnodes,fnodes,w)
    xeval = xeval(:);
    p = zeros(size(xeval));
    for q = 1:numel(xeval)
        x = xeval(q);
        diff = x - xnodes(:).';
        hit = find(abs(diff) < 1e-14,1);
        if ~isempty(hit)
            p(q) = fnodes(hit);
        else
            tmp = w(:).' ./ diff;
            p(q) = sum(tmp .* fnodes(:).') / sum(tmp);
        end
    end
end
```

### 9.3 GLOF 系数计算

设 $I_Nf$ 是由 Chebyshev 配点值构造的插值函数。令

$$
F_N(x)=I_Nf(x).
$$

在 GLOF 求积节点上评价

$$
F_q=(I_Nf)(x_q).
$$

然后使用同一个 Gauss--GLOF 系数公式。

```matlab
function g = glof_coeffs_from_cheb_gaussglof(xnodes,fnodes,wbar,m,alpha,beta,lambda,Q)
    [xq,Wg,S,~,~,~,gamma_vec] = glof_quadrature_matrices(Q,m,alpha,beta,lambda);
    Fq = bary_eval(xq,xnodes,fnodes,wbar);
    g = (S * (Wg(:).*Fq(:))) ./ gamma_vec(:);
    g = real(g);
end
```

---

## 10. GLOF 重构评价

```matlab
function R = eval_glof_reconstruction(x,g,alpha,beta,lambda)
    m = length(g)-1;
    theta = (beta-lambda)/2;
    xrow = x(:).';
    idx = (xrow > 0);
    y = zeros(size(xrow));
    y(idx) = -(beta+1)*log(xrow(idx));
    y(~idx) = 0;
    L = laguerre_matrix(m,alpha,y);
    S = (xrow.^theta).*L;
    S(:,~idx) = 0;  % theta>0, so S_n(0)=0
    R = real(g(:).'*S);
    R = reshape(R,size(x));
end
```

---

## 11. 求积阶数 $Q$ 的选择与检查

由于 $F_N$ 来自 Fourier 或 Chebyshev 标准数据，一般不属于 GLOF 求积的精确空间，所以必须检查 $Q$-收敛。

建议初始取

$$
Q=\max\{2N+20,\,2m+2k+20\}.
$$

更稳妥的选择是

$$
Q=3N\quad\text{or}\quad Q=4N.
$$

但要注意：如果 $\beta$ 固定且 $Q$ 很大，则

$$
x_{\min}\approx e^{-O(Q)/(\beta+1)}
$$

可能非常小，导致 $x_q^{\lambda-\beta}$ 或等价的 $x_q^{-\theta}$ 产生数值放大。因此：

1. 固定参数实验中，不要盲目把 $Q$ 取到几百以上；
2. $N$-dependent 参数 $\beta_N=2r+bN$ 会显著缓解节点过度靠近零的问题；
3. 每次实验都计算

$$
\frac{\|g^{(Q)}-g^{(2Q)}\|_2}{\|g^{(2Q)}\|_2+10^{-300}}
$$

作为求积误差诊断。

MATLAB 检查：

```matlab
function rel = check_quadrature_convergence_fourier(ell,c,m,alpha,beta,lambda,Q)
    g1 = glof_coeffs_from_fourier_gaussglof(ell,c,m,alpha,beta,lambda,Q);
    g2 = glof_coeffs_from_fourier_gaussglof(ell,c,m,alpha,beta,lambda,2*Q);
    rel = norm(g1-g2) / (norm(g2) + 1e-300);
end
```

---

## 12. 纯模型的精确性 sanity check

当

$$
\theta=r
$$

时，

$$
f(x)S_n^{(\alpha,\beta,\lambda)}(x) = x^{2r}(-\log x)^k L_n^{(\alpha)}(y).
$$

由于

$$
\beta-\lambda=2r,
$$

上式属于

$$
x^{\beta-\lambda}\operatorname{span}\{1,\log x,\ldots,(\log x)^{n+k}\}.
$$

因此 Gauss--GLOF 求积对精确函数计算 GLOF 系数时，只要

$$
2Q+1\ge m+k,
$$

就应该达到机器精度。这是检验代码是否正确的第一步。

实现：

```matlab
function g = glof_coeffs_from_exact_function_gaussglof(r,k,m,alpha,beta,lambda,Q)
    [xq,Wg,S,~,~,~,gamma_vec] = glof_quadrature_matrices(Q,m,alpha,beta,lambda);
    Fq = f_logsing(xq,r,k);
    g = (S * (Wg(:).*Fq(:))) ./ gamma_vec(:);
    g = real(g);
end
```

检查：

```matlab
r = 0.5; k = 1;
alpha = 0; beta = 2*r; lambda = 0;
m = 8; Q = 20;
g_exact = glof_coeffs_from_exact_function_gaussglof(r,k,m,alpha,beta,lambda,Q);

xref = linspace(0,1,5001);
fref = f_logsing(xref,r,k);
R = eval_glof_reconstruction(xref,g_exact,alpha,beta,lambda);
err = max(abs(fref - R));
fprintf('direct GLOF projection error = %.3e\n',err);
```

预期：若 $m\ge k$，误差应接近机器精度。

---

## 13. 主实验流程

### 13.1 准备

```matlab
r = 0.5;
k = 1;
Ns = [16 32 64 128 256];
xref = linspace(0,1,5001);
fref = f_logsing(xref,r,k);
```

### 13.2 主循环

```matlab
EF  = zeros(size(Ns));   % Fourier partial sum error
EC  = zeros(size(Ns));   % Chebyshev interpolation error
EGF = zeros(size(Ns));   % GLOF from Fourier data
EGC = zeros(size(Ns));   % GLOF from Chebyshev data
QrelF = zeros(size(Ns));
QrelC = zeros(size(Ns));

for ii = 1:length(Ns)
    N = Ns(ii);

    % ---------------------------
    % parameter choice
    % ---------------------------
    mode = 'Ndependent';  % 'fixed' or 'Ndependent'
    if strcmp(mode,'fixed')
        alpha = 0;
        beta = 2*r;
        lambda = 0;
        m = max(k,ceil(0.05*N));
        Q = max(40,2*N+20);
    else
        a = 0.05;
        b = 0.05;
        gamma_m = 0.05;
        alpha = a*N;
        beta = 2*r + b*N;
        lambda = b*N;
        m = max(k,ceil(gamma_m*N));
        Q = max(60,3*N);
    end

    % ---------------------------
    % Fourier data
    % ---------------------------
    [ell,c] = fourier_coeffs_quad(N,r,k);
    fF = eval_fourier_on_nodes(xref(:),ell,c);
    EF(ii) = max(abs(fref(:) - fF(:)));

    % GLOF reconstruction from Fourier coefficients
    gF = glof_coeffs_from_fourier_gaussglof(ell,c,m,alpha,beta,lambda,Q);
    RF = eval_glof_reconstruction(xref,gF,alpha,beta,lambda);
    EGF(ii) = max(abs(fref(:) - RF(:)));

    % quadrature convergence check
    gF2 = glof_coeffs_from_fourier_gaussglof(ell,c,m,alpha,beta,lambda,2*Q);
    QrelF(ii) = norm(gF-gF2)/(norm(gF2)+1e-300);

    % ---------------------------
    % Chebyshev collocation data
    % ---------------------------
    xnodes = cheb_lobatto_nodes(N);
    fnodes = f_logsing(xnodes,r,k);
    wbar = bary_weights_cheb_lobatto(N);

    fC = bary_eval(xref(:),xnodes,fnodes,wbar);
    EC(ii) = max(abs(fref(:) - fC(:)));

    % GLOF reconstruction from Chebyshev values
    gC = glof_coeffs_from_cheb_gaussglof(xnodes,fnodes,wbar,m,alpha,beta,lambda,Q);
    RC = eval_glof_reconstruction(xref,gC,alpha,beta,lambda);
    EGC(ii) = max(abs(fref(:) - RC(:)));

    % quadrature convergence check
    gC2 = glof_coeffs_from_cheb_gaussglof(xnodes,fnodes,wbar,m,alpha,beta,lambda,2*Q);
    QrelC(ii) = norm(gC-gC2)/(norm(gC2)+1e-300);
end
```

### 13.3 画图

```matlab
figure;
semilogy(Ns,EF,'-o'); hold on;
semilogy(Ns,EC,'-s');
semilogy(Ns,EGF,'-^');
semilogy(Ns,EGC,'-d');
grid on;
legend('Fourier partial sum','Chebyshev interp', ...
       'GLOF from Fourier','GLOF from Chebyshev', ...
       'Location','southwest');
xlabel('N'); ylabel('max error');
title(sprintf('f(x)=x^{%.2f}(-log x)^{%d}',r,k));
```

求积误差诊断：

```matlab
figure;
semilogy(Ns,QrelF,'-o'); hold on;
semilogy(Ns,QrelC,'-s');
grid on;
legend('Fourier data: ||g_Q-g_{2Q}||','Cheb data: ||g_Q-g_{2Q}||');
xlabel('N'); ylabel('relative coefficient difference');
title('Gauss-GLOF quadrature convergence check');
```

---

## 14. 必须报告的误差量

建议每组 $(r,k)$ 报告：

1. $E_N^{\rm F}=\|f-f_N^{\rm F}\|_\infty$；
2. $E_N^{\rm C}=\|f-I_Nf\|_\infty$；
3. $E_N^{\rm GF}=\|f-R_{N,m}^{\rm GLOF,F}\|_\infty$；
4. $E_N^{\rm GC}=\|f-R_{N,m}^{\rm GLOF,C}\|_\infty$；
5. $\|g_Q-g_{2Q}\|_2/\|g_{2Q}\|_2$，用于说明 GLOF 系数积分已经充分收敛。

如果第 5 项不小，则不能把误差归因于重构方法；首先要增加或调整 $Q$、$\alpha$、$\beta$。

---

## 15. 数值稳定性注意事项

### 15.1 不要让 GLOF 求积节点过度靠近零

固定 $\beta$ 且大 $Q$ 时，最大 Laguerre 节点大约为 $O(Q)$，于是

$$
x_{\min}\approx e^{-O(Q)/(\beta+1)}.
$$

如果 $x_{\min}$ 下溢为零，GLOF 权重中的 $x_q^{\lambda-\beta}$ 会出问题。

解决方法：

1. 先使用中等 $Q$，例如 $Q=2N+20$；
2. 使用 $N$-dependent $\beta_N=2r+bN$，让节点分布不过度压向零；
3. 若仍需很大 $Q$，改用 log-weight 形式实现加权乘积。

### 15.2 检查 $F_N(x_q)/x_q^r$ 的行为

在调谐情形 $\theta=r$ 下，GLOF 系数公式等价于对

$$
F_N(x_q)x_q^{-r}
$$

进行加权求和。精确函数满足 $f(x)x^{-r}=(-\log x)^k$，但 Fourier 或 Chebyshev 近似 $F_N$ 在极小 $x_q$ 处可能不完全满足这种端点行为。

因此，如果观察到重构误差异常大，应检查：

```matlab
ratio = Fq ./ (xq.^r);
semilogy(xq,abs(ratio));
```

若 `ratio` 在极小节点处爆炸，说明标准谱近似在端点附近不稳定，或者 $Q$ 过大导致采样到过小的节点。

### 15.3 推荐先使用 $N$-dependent 参数

在后处理实验中，建议优先使用

$$
\alpha_N=aN,
\qquad
\beta_N=2r+bN,
\qquad
\lambda_N=bN.
$$

这既保持 $\theta=r$，又可以改善求积节点分布和截断误差行为。

---

## 16. 进一步测试函数

纯模型

$$
f(x)=x^r(-\log x)^k
$$

是必要的第一步，但它在调谐 GLOF 空间中有限维可表示。完成 sanity check 后，应继续测试：

### 16.1 解析乘子

$$
f(x)=e^x x^r(-\log x)^k.
$$

### 16.2 多项式解析乘子

$$
f(x)=(1+x+x^2)x^r(-\log x)^k.
$$

### 16.3 解析部分加奇异部分

$$
f(x)=1+x+x^r(-\log x)^k.
$$

此时单一 GLOF 块无法自然表示常数项，建议使用

$$
\operatorname{span}\{1\}
\oplus \mathcal V_m^{(1)}
\oplus \mathcal V_m^{(r)}.
$$

---

## 17. 最小实现清单

建议按以下顺序实现：

1. `f_logsing.m`
2. `gauss_laguerre.m`
3. `laguerre_matrix.m`
4. `glof_quadrature_matrices.m`
5. `eval_glof_reconstruction.m`
6. `glof_coeffs_from_exact_function_gaussglof.m`
7. `fourier_coeffs_quad.m`
8. `eval_fourier_on_nodes.m`
9. `glof_coeffs_from_fourier_gaussglof.m`
10. `cheb_lobatto_nodes.m`
11. `bary_weights_cheb_lobatto.m`
12. `bary_eval.m`
13. `glof_coeffs_from_cheb_gaussglof.m`
14. `main_experiment.m`

第一阶段只做 $r=1/2,k=1$，并确认：

1. 精确 GLOF 投影达到机器精度；
2. Fourier / Chebyshev 传统近似误差不是指数衰减；
3. 从 Fourier / Chebyshev 数据得到的 GLOF 后处理误差明显改善；
4. $Q$-收敛诊断足够小。

---

## 18. 核心结论

本实验的关键不是单纯实现 GLOF 展开，而是实现如下后处理链条：

$$
\boxed{
\text{Fourier coefficients or Chebyshev values}
\longrightarrow
F_N(x_q)
\longrightarrow
\text{Gauss--GLOF quadrature for }g_n
\longrightarrow
\sum_{n=0}^m g_nS_n^{(\alpha,\beta,\lambda)}(x).
}
$$

其中中间的积分计算必须使用与 log mapping 匹配的 Gauss--LOF / Gauss--GLOF 求积。这样才能既保持效率，又准确处理 $x=0$ 附近的对数端点奇异性。
