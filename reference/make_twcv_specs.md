# Create a TWCV specification

Defines one TWCV weighting configuration, including balancing variables,
discretization level, and shrinkage strength.

## Usage

``` r
make_twcv_specs(
  predictor_vars,
  include_distance = TRUE,
  balance_by = 0.2,
  shrink_lambda = 0.2,
  name = "twcv_extended"
)
```

## Arguments

- predictor_vars:

  Character vector of predictor variables.

- include_distance:

  Logical; include prediction distance `d`.

- balance_by:

  Quantile spacing used for discretization.

- shrink_lambda:

  Shrinkage parameter for calibration weights.

- name:

  Name of the specification.

## Value

Named list with one TWCV specification.
