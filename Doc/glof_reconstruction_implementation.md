# Implementation plan for GLOF reconstruction of logarithmic endpoint singularities

This file is intended as a programming-oriented companion to the LaTeX note.  The target model is

\[
f(x)=a(x)+b(x)x^r(-\log x)^k,\qquad x\in[0,1],\quad r>0,
\]

with analytic `a` and `b`.  The singular endpoint is `x=0`.

---

## 1. Baseline experiment: show the loss of accuracy of standard spectral data

### 1.1 Choose manufactured test functions

Use functions with known exact values, for example

```text
f1(x) = cos(x) + sin(x) * x^r * (-log(x))^k
f2(x) = 1 + x + exp(x) * x^r * (-log(x))^k
```

At `x=0`, define the singular term as zero because `x^r (-log x)^k -> 0` for `r>0`.

Recommended starting parameters:

```text
r = 0.35, 0.5, 0.8, 1.3
k = 1 or 2
N = 32, 64, 128, 256, 512
```

Use a dense evaluation grid, for example 5000 to 20000 points clustered near `x=0`.

---

## 2. Traditional Chebyshev collocation baseline

### 2.1 Nodes

Use Chebyshev-Lobatto points mapped from `[-1,1]` to `[0,1]`:

```text
xi_j = cos(pi*j/N),       j = 0,...,N
x_j  = (1 + xi_j)/2
```

This includes `x=0` and `x=1`.  If evaluating `(-log x)^k`, handle `x=0` by assigning the limiting value of the full singular term `x^r(-log x)^k`, namely zero.

### 2.2 Interpolant

Use barycentric interpolation or a DCT-based Chebyshev coefficient transform.

Programming choices:

```python
# robust and simple
scipy.interpolate.BarycentricInterpolator

# faster for many N
scipy.fft.dct for Chebyshev coefficients
```

### 2.3 Diagnostics

For each `N`:

1. Build the Chebyshev interpolant from `f(x_j)`.
2. Evaluate it on the dense grid.
3. Plot `f_N(x)` and `f(x)` near `x=0`.
4. Plot the pointwise error on a logarithmic scale.
5. Record

```text
max_error = max(abs(f_N(x_dense) - f(x_dense)))
endpoint_error_window = max error on x in [0, 0.05]
```

Expected behavior: the convergence is generally algebraic rather than exponential because the endpoint derivatives are singular.

---

## 3. Traditional Fourier baseline

### 3.1 Uniform grid

Use `M = 2N+1` or an FFT-friendly even number `M` on `[0,1)`:

```text
x_j = j/M,     j = 0,...,M-1
```

### 3.2 Fourier coefficients and partial sum

Compute coefficients by FFT:

```python
fvals = f(x)
fhat = fft(fvals) / M
```

To evaluate a truncated partial sum on a dense grid, either use direct summation

```python
sum_{|ell|<=N} fhat_ell * exp(2*pi*i*ell*x)
```

or zero-pad Fourier coefficients and use an inverse FFT.

### 3.3 Diagnostics

Because the periodic extension of a non-periodic function usually has a mismatch between `x=0` and `x=1`, Fourier partial sums may show Gibbs-type oscillations near the endpoints.  Plot:

```text
f_N(x) vs f(x)
pointwise error in log scale
zoom near x=0 and x=1
```

If the goal is to isolate only the logarithmic singularity at `x=0`, also test a periodic-compatible manufactured function, for example by subtracting endpoint mismatches or multiplying by a smooth cutoff.

---

## 4. GLOF basis construction

The GLOF basis is

\[
\Phi_n^{(\theta;A,B)}(x)
= x^\theta L_n^{(A)}((B+1)(-\log x)),
\]

where `L_n^(A)` is the generalized Laguerre polynomial.

For the model

\[
f(x)=a(0)+x a_1(x)+x^r(-\log x)^k b(x),
\]

use the enriched dictionary

```text
constant mode: 1
regular analytic remainder modes: theta = 1
log-singular modes: theta = r
```

That is,

```text
V_N = span{1}
      + span{x * L_n^A((B1+1)*(-log x)), n=0,...,m1}
      + span{x^r * L_n^A((Br+1)*(-log x)), n=0,...,mr}
```

At `x=0`, set all `theta>0` GLOF modes to zero.

---

## 5. Parameter choices

Start with conservative values:

```text
B = b*N,        b in [0.05, 0.5]
m = gamma*N,    gamma in [0.02, 0.25]
```

Usually require

```text
gamma < b
```

because the regularization mechanism improves when the scaling parameter `B` is larger than the retained degree `m`.

For Chebyshev experiments, start with

```text
A = 0
```

because `(-log x)^k` is analytic at `x=1` for integer `k` and the singular endpoint is only `x=0`.

For Fourier-periodic data, consider

```text
A = a*N,        a in [0.02, 0.2]
```

so that the factor `(-log x)^A ~ (1-x)^A` suppresses the endpoint `x=1`.  This is more delicate numerically and theoretically.

---

## 6. Practical reconstruction algorithm: data-only least squares

This is the recommended implementation when only values or Fourier coefficients of `f` are available and the decomposition into `a_1` and `b` is not known.

### 6.1 Build sample values

For Chebyshev data, use the original Chebyshev values or oversample by evaluating the Chebyshev interpolant on a larger grid.

For Fourier data, evaluate the Fourier partial sum on an oversampled grid.

Recommended oversampling:

```text
J >= 2 * number_of_basis_functions
J = 4 * number_of_basis_functions is safer
```

Use sample points clustered near `x=0`, for example Chebyshev-Lobatto points or mapped GLOF/Laguerre points.

### 6.2 Build the design matrix

For each sample point `x_j`, form columns:

```text
A[:, 0] = 1
A[:, 1+n] = x_j * L_n^(A1)((B1+1)*(-log x_j)),          n=0,...,m1
A[:, ...] = x_j^r * L_n^(Ar)((Br+1)*(-log x_j)),         n=0,...,mr
```

Use `scipy.special.eval_genlaguerre(n, A, z)` or recurrence relations for Laguerre polynomials.

### 6.3 Scaling and stabilization

Column scaling is essential:

```python
col_norms = np.linalg.norm(A, axis=0)
A_scaled = A / col_norms
```

Solve by QR or SVD, not by normal equations:

```python
coef_scaled, *_ = scipy.linalg.lstsq(W @ A_scaled, W @ f_samples)
coef = coef_scaled / col_norms
```

If the matrix is ill-conditioned, use truncated SVD or Tikhonov regularization:

```python
min ||W(Ac-d)||_2^2 + lambda_reg^2 ||c||_2^2
```

Record the condition number for every run.

---

## 7. Projection-based reconstruction for manufactured tests

When `a`, `b`, and hence `u_1`, `u_r` are known, test the theory componentwise.

For a block

\[
u_\theta(x)=x^\theta(-\log x)^\kappa c(x),
\]

use

\[
\Phi_n^{(\theta;A,B)}(x)=x^\theta L_n^{(A)}((B+1)(-\log x)),
\]

and compute exact or quadrature-based coefficients

\[
\widehat u_n
=\frac{1}{h_n^{(A,B)}}\int_0^1u_\theta(x)\Phi_n^{(\theta;A,B)}(x)(-\log x)^A x^{B-2\theta}\,dx.
\]

Use Gauss-Laguerre quadrature after the map

```text
z = (B+1)*(-log x)
x = exp(-z/(B+1))
```

The quadrature weight is `z^A exp(-z)`.  This is often more stable than integrating directly in `x`.

---

## 8. Pseudocode for the enriched least-squares GLOF postprocessor

```python
def glof_reconstruct_from_samples(x_samples, f_samples, r, N,
                                  A1=0, Ar=0,
                                  b1=0.2, br=0.2,
                                  gamma1=0.08, gammar=0.08,
                                  weights=None,
                                  lambda_reg=0.0):
    B1 = b1 * N
    Br = br * N
    m1 = int(round(gamma1 * N))
    mr = int(round(gammar * N))

    # avoid log(0); handle x=0 by limits
    y = np.zeros_like(x_samples)
    mask = x_samples > 0
    y[mask] = -np.log(x_samples[mask])

    columns = []
    columns.append(np.ones_like(x_samples))

    # theta = 1 block
    z1 = (B1 + 1.0) * y
    for n in range(m1 + 1):
        col = np.zeros_like(x_samples)
        col[mask] = x_samples[mask] * eval_genlaguerre(n, A1, z1[mask])
        columns.append(col)

    # theta = r block
    zr = (Br + 1.0) * y
    for n in range(mr + 1):
        col = np.zeros_like(x_samples)
        col[mask] = (x_samples[mask] ** r) * eval_genlaguerre(n, Ar, zr[mask])
        columns.append(col)

    A = np.column_stack(columns)

    # optional row weights
    if weights is None:
        Wsqrt = np.ones_like(f_samples)
    else:
        Wsqrt = np.sqrt(weights)

    Aw = A * Wsqrt[:, None]
    dw = f_samples * Wsqrt

    # column scaling
    col_norms = np.linalg.norm(Aw, axis=0)
    col_norms[col_norms == 0] = 1.0
    Aws = Aw / col_norms

    if lambda_reg == 0.0:
        coef_scaled, *_ = scipy.linalg.lstsq(Aws, dw)
    else:
        top = Aws
        bottom = lambda_reg * np.eye(Aws.shape[1])
        rhs = np.concatenate([dw, np.zeros(Aws.shape[1])])
        coef_scaled, *_ = scipy.linalg.lstsq(np.vstack([top, bottom]), rhs)

    coef = coef_scaled / col_norms
    return coef, dict(B1=B1, Br=Br, m1=m1, mr=mr, A1=A1, Ar=Ar,
                      b1=b1, br=br, gamma1=gamma1, gammar=gammar)
```

Evaluation uses the same dictionary columns at new points and multiplies them by `coef`.

---

## 9. What to report in numerical experiments

For each method and each `N`, report:

```text
1. max-norm error on a dense grid
2. max error near x=0
3. weighted L2 error if quadrature weights are available
4. condition number of the reconstruction matrix
5. chosen parameters A, B, m
6. sensitivity to noise in data
```

Suggested plots:

```text
1. baseline Chebyshev/Fourier approximation vs exact f
2. baseline pointwise error in log scale
3. GLOF reconstruction pointwise error in log scale
4. max error vs N in semilog scale
5. condition number vs N
6. parameter sweep over b and gamma
```

---

## 10. Numerical cautions

1. Do not evaluate `log(0)`.  Use the endpoint limits of the full basis functions.
2. Laguerre polynomials can overflow for large `n` or large `z`.  Use recurrence with scaling if needed.
3. Always scale columns before least squares.
4. Oversampling is strongly recommended.
5. The optimal parameters are problem-dependent.  A small grid search over `(b, gamma)` is useful.
6. Fourier data may include artificial endpoint jumps caused by periodic extension.  This can mask the logarithmic singularity unless handled carefully.
7. A pure LOF basis is not recommended for uniform endpoint reconstruction because its modes grow like powers of `-log(x)` near zero.  Use GLOFs plus a separate constant mode.
