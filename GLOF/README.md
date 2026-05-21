# MATLAB experiments: GLOF reconstruction of $f(x)=x^r(-\log x)^k$

This folder implements the numerical experiments described in
`../glof_matlab_experiments_gauss_glof_quadrature.md`.

The model function is

$$
  f(x) = x^r(-\log x)^k, \qquad x \in [0,1], \quad r > 0, \quad k \in \mathbb{N},
$$

with the convention $f(0) = 0$.

There are two goals:

1. Show that the standard truncated Fourier series and Chebyshev
   collocation interpolant lose spectral accuracy because $f$ has a
   logarithmic singularity at the left endpoint.
2. Demonstrate that a Generalized Log-Orthogonal-Function (GLOF) post
   processing of the same standard spectral data recovers high accuracy.

## Folder layout

```
matlab/
  README.md                          this file
  main_experiment.m                  master driver, all cases together
  f_logsing.m                        f(x) = x^r * (-log x)^k
  fourier_coeffs_quad.m              accurate Fourier coefficients via integral()
  fourier_coeffs_fft.m               FFT-sampled Fourier coefficients
  eval_fourier_sum.m                 evaluate a truncated Fourier series
  cheb_lobatto_nodes.m               nodes on [0,1] in ascending order
  bary_weights_cheb_lobatto.m        barycentric weights for those nodes
  bary_eval.m                        barycentric interpolation evaluator
  laguerre_matrix.m                  generalized Laguerre polynomials, recurrence
  gauss_laguerre.m                   Golub-Welsch Gauss-Laguerre nodes/weights
  glof_quadrature_matrices.m         Chen-Shen Gauss-GLOF nodes/weights/matrices
  glof_norm_const.m                  gamma_n^{(alpha,beta)} (log-space stable)
  glof_matrix.m                      evaluate S_n^{(alpha,beta,lambda)}(x)
  glof_coeffs_exact.m                sanity-check projection from exact f
  glof_coeffs_from_fourier.m         g_n from truncated Fourier coefficients
  glof_coeffs_from_cheb_interp.m     g_n from Chebyshev values via Gauss-GLOF
  glof_coeffs_from_cheb_weighted.m   backward-compatible wrapper
  eval_glof_reconstruction.m         evaluate sum_n g_n S_n(x)

  tests/
    add_paths.m                      add ../ to MATLAB path
    run_all_tests.m                  driver, runs all test_*.m in order
    test_01_f_logsing.m              f(x) sanity check
    test_02_laguerre_and_quadrature  Laguerre orthogonality and quadrature
    test_03_glof_orthogonality       S_n L^2 orthogonality at numeric level
    test_04_glof_exact_projection    exact projection: roundoff for m >= k
    test_05_fourier_baseline         traditional Fourier failure
    test_06_chebyshev_baseline       traditional Chebyshev failure
    test_07_glof_from_fourier        GLOF recovery from Fourier data
    test_08_glof_from_cheb_A         GLOF Variant A from Chebyshev data
    test_09_glof_from_cheb_B         compatibility wrapper for old Variant B name
    test_10_glof_analytic_multiplier single-block limitation on e^x * x^r * (-log x)^k
```

## How to run

### Run the full experiment

From the `matlab/` folder in MATLAB:

```matlab
main_experiment
```

This loops over four test cases and several truncation sizes $N$, and
produces a summary figure plus a text table of maximum errors for

* the bare Fourier partial sum,
* the bare Chebyshev interpolant,
* GLOF reconstruction from Fourier coefficients (called `GLOF-F`),
* GLOF reconstruction from Chebyshev values (`GLOF-C`),
* the relative Gauss--GLOF quadrature diagnostics
  `||g_Q-g_{2Q}||_2/||g_{2Q}||_2`.

The top of the script has one switch that you may toggle:

* `mode = 'Ndependent'` (Chen--Shen-style parameters, $N$-dependent)
  or `'fixed'` (fixed diagnostic parameters).

### Step-by-step diagnostics

Each numbered script under `tests/` validates one stage of the pipeline.
Running them in order is the easiest way to verify the implementation:

| Script | What it checks |
| --- | --- |
| `test_01_f_logsing` | model function values, $f(0)=f(1)=0$ |
| `test_02_laguerre_and_quadrature` | Laguerre orthogonality and Gauss--Laguerre exactness |
| `test_03_glof_orthogonality` | $\int_0^1 S_m S_n w_{\alpha,\lambda} dx = \gamma_n\delta_{mn}$ numerically |
| `test_04_glof_exact_projection` | direct projection of $f$ is exact for $m\ge k$ in the tuned basis |
| `test_05_fourier_baseline` | algebraic convergence of the Fourier partial sum |
| `test_06_chebyshev_baseline` | algebraic convergence of Chebyshev interpolation |
| `test_07_glof_from_fourier` | GLOF recovery from $(2N+1)$ Fourier coefficients |
| `test_08_glof_from_cheb_A` | GLOF recovery from Chebyshev values using Gauss--GLOF |
| `test_09_glof_from_cheb_B` | compatibility wrapper for the old Variant-B function name |
| `test_10_glof_analytic_multiplier` | failure mode: single-block GLOF for $e^x \cdot x^r(-\log x)^k$ |

Or run all of them in one go:

```matlab
cd tests
run_all_tests
```

## Parameter conventions

The GLOF basis used here is

$$
  S_n^{(\alpha,\beta,\lambda)}(x) = x^{(\beta-\lambda)/2}
    L_n^{(\alpha)}\bigl(-(\beta+1)\log x\bigr),
$$

orthogonal on $[0,1]$ with respect to $(-\log x)^\alpha x^\lambda$.

For the singular function $x^r(-\log x)^k$ we always tune

$$
  \theta := \frac{\beta-\lambda}{2} = r,
$$

so that the prefactor $x^\theta$ matches the singularity exponent
exactly.  The Laguerre family then spans every polynomial in $-\log x$
and the *direct* GLOF projection is exact once $m\ge k$.

The two parameter families used by the scripts are:

* **fixed**: $\alpha = 0$, $\beta = 2r$, $\lambda = 0$, $m = \max(k,\lceil 0.05N\rceil)$.
* **N-dependent** (Chen--Shen style): $\alpha_N = aN$, $\beta_N = 2r+bN$,
  $\lambda_N = bN$, $m = \max(k,\lceil \gamma_m N\rceil)$ with
  $a = b = \gamma_m = 0.05$ as the default starting values.

## Performance notes

* `glof_coeffs_from_fourier.m`, `glof_coeffs_from_cheb_interp.m`, and the
  compatibility wrapper `glof_coeffs_from_cheb_weighted.m` all use the same
  Chen--Shen Gauss--GLOF quadrature matrices.  They do not call
  MATLAB's adaptive `integral` for GLOF coefficients.
* `fourier_coeffs_quad.m` does one `integral` call per requested Fourier
  mode and accepts a function handle plus a mode vector such as `(-N:N).'`.
  For larger $N$ switch to `fourier_coeffs_fft.m` with a large $M$.
* All GLOF post-processing routines use only matrix products and the
  Gauss--Laguerre Jacobi-matrix eigendecomposition and remain fast.

## Expected behaviour

For the pure model $f = x^r(-\log x)^k$ with tuned GLOF parameters the
direct projection is exact, so the nontrivial question is whether the
GLOF coefficients computed *from truncated Fourier or Chebyshev data*
converge to those exact coefficients fast enough.

The reference run (`main_experiment` with `mode='Ndependent'`) reports the
plain Fourier and Chebyshev errors, the two Gauss--GLOF post-processed
errors, and the required quadrature diagnostics
`||g_Q-g_{2Q}||_2/||g_{2Q}||_2`.  If these diagnostics are not small, first
increase or retune $Q,\alpha,\beta$ before interpreting reconstruction
errors.

### Parameter pitfall

For the pure model the exact GLOF expansion has only $k+1$ nonzero
coefficients, so `main_experiment` and the diagnostic scripts default to
$m = k$.  Larger $m$ is useful for genuinely multi-mode functions, but it
amplifies tiny coefficient
errors because the basis functions $S_n(x) = x^r L_n^{(\alpha)}(y)$ grow
like $[(\beta+1)/r]^n$ near $x = 0$; for $\beta = 7.4$ and $r = 0.5$
that is already $16.8^7 \sim 10^8$ at $n = 7$.  This is documented in
the comments of the test scripts.

The same endpoint conditioning affects the quadrature order.  Very large
$Q$ can place Gauss--GLOF nodes far below the resolution scale of the
Fourier partial sum, making $F_N(x_q)/x_q^r$ enormous and polluting high
GLOF modes.  The default scripts therefore use moderate $Q$ for the pure
model and report the `Qrel` diagnostics so this failure mode is visible.

For functions with more genuine multi-mode content (Taylor multiplier,
mixed analytic/singular parts) the multi-block decomposition of §13.2
of the spec is needed; `test_10` demonstrates the failure mode that
motivates it.
