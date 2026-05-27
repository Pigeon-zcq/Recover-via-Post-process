# Gegenbauer Python conversion

This folder now contains a high-precision Python translation of the supplied
Fortran Gegenbauer reconstruction programs.

## Files

- `gegenbauer_reconstruction.py`: core mpmath implementation.
- `run_reconstruction.py`: command-line runners for paper examples, the
  periodic `fs` path, parameter search, and the Fourier-Galerkin transport
  path.
- `validate_examples.py`: compares Chen-Shu 2015 Examples 6.1 and 6.2 against
  the published linear-choice maximum-error table.
- `reproduce_chen_shu_2015_examples.py`: reproduces Chen-Shu 2015 Examples 6.1
  and 6.2 tables plus Fig. 1 pointwise-error plots.
- `test_reconstruction.py`: lightweight regression tests.

## Fortran mapping

- `fs`, `fv` -> `fs_periodic_singularity`, `fv_transport_exact`
- `transf` -> `transform_y_to_x`, `transform_x_to_y`
- `getPointValuefromFC` -> `composite_gauss_fourier_coefficients`
- `gegen`, `gegenP` -> `gegenbauer_grid_values`, `gegenbauer_value`
- `gcoe` -> `gegenbauer_coefficients_from_fourier`
- `ccoe` -> `collocation_spectral_coefficients`
- `FourierGalerkinTransport` -> `fourier_galerkin_transport_coefficients`
- `VariableCoeffsLinearTransportOptimize.f90` parameter scan ->
  `find_best_parameters`

## Usage

Run the unit tests:

```bash
python -m unittest discover -s Gegenbauer_Python -p 'test_*.py'
```

Validate the published linear-choice rows for Examples 6.1 and 6.2:

```bash
python Gegenbauer_Python/validate_examples.py --example 6.1 6.2 --N 40 80 --dps 60
```

Generate reproduction tables and Fig. 1 images for Examples 6.1 and 6.2:

```bash
python Gegenbauer_Python/reproduce_chen_shu_2015_examples.py --quick --dps 60
python Gegenbauer_Python/reproduce_chen_shu_2015_examples.py --full --dps 80
```

Run a single paper example:

```bash
python Gegenbauer_Python/run_reconstruction.py paper --example 6.1 --N 40 --dps 60
```

Search for an optimal pair:

```bash
python Gegenbauer_Python/run_reconstruction.py scan --example 6.1 --N 40 --max-lambda 20 --max-m 20
```

Run the variable-coefficient transport reconstruction:

```bash
python Gegenbauer_Python/run_reconstruction.py transport --N 20 --lambda 1 --m 5
```
