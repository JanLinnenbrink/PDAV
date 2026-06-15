# Generate buffered leave-one-out tasks

Constructs buffered LOO validation tasks by excluding observations
within a radius around each test point.

## Usage

``` r
generate_buffered_loo_tasks(
  sample_dat,
  buffer_radii = NULL,
  target_d = NULL,
  n_candidates = 30,
  n_tasks = NULL,
  radius_correction = 0.6,
  max_holdout_frac = 0.2,
  min_train_n = floor(nrow(sample_dat)/2),
  max_dist = NULL,
  include_zero = TRUE,
  verbose = 0,
  seed = NULL
)
```

## Arguments

- sample_dat:

  Data with coordinates and IDs.

- buffer_radii:

  Optional vector of radii.

- target_d:

  Optional target distance distribution.

- n_candidates:

  Number of candidate radii per point.

- n_tasks:

  Optional number of tasks to retain.

- radius_correction:

  Scaling factor for derived radii.

- max_holdout_frac:

  Maximum fraction of excluded data.

- min_train_n:

  Minimum training size.

- max_dist:

  Maximum allowed prediction distance.

- include_zero:

  Logical; include zero radius.

- verbose:

  Verbosity level.

- seed:

  Optional random seed.

## Value

Object of class `"twcv_buffered_tasks"`.
