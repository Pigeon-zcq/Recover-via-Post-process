# GLOF Python

High-precision Python implementation of the Fourier-to-GLOF reconstruction
workflow in `code.md`.

The core implementation lives in `glof_reconstruction.py` and uses `mpmath`
for quadrature nodes, weights, Laguerre recurrences, and complex arithmetic.

The default Fourier-to-GLOF coefficient path uses the stable direct
`B_{ell n}` formula from `code.md`; it does not first evaluate the Fourier
partial sum at tiny GLOF quadrature nodes. The older node-value path is still
available with `coefficient_method="nodes"` for comparisons.

`Q` controls the ordinary GLOF quadrature used by the node-value path.
`QB` controls the Laguerre-Gauss quadrature used for the stable
`B_{ell n}` integrals.  Use `adaptive_b=True` to compare successive `QB`
values and report convergence diagnostics.

## Quick Demo

```bash
python GLOF_Python/run_glof_reconstruction.py --NF 12 --r 0.5 --k 1 --alpha 0.5 --dps 80
python GLOF_Python/run_glof_reconstruction.py --NF 20 --QB 160 --adaptive-B --QB-max 320
```

## Tests

```bash
python -m unittest discover -s GLOF_Python -p 'test_*.py'
```

## Minimal API

```python
from GLOF_Python.glof_reconstruction import (
    fourier_coefficients,
    fourier_glof_reconstruct,
    linspace,
    model_log_singularity,
)

f = model_log_singularity(r="0.5", k=1)
f_hat = fourier_coefficients(f, nf=12, dps=80)
z = linspace("0.005", "0.995", 80)
values = fourier_glof_reconstruct(
    f_hat,
    r="0.5",
    z_eval=z,
    k=1,
    alpha="0.5",
    theta="0.05",
    gamma="0.1",
    QB=96,
    dps=80,
    coefficient_method="stable_b",
    adaptive_b=True,
    q_b_max=192,
    b_tolerance="1e-8",
)
```
