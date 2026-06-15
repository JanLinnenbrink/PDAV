# Estimate deployment-oriented predictive performance from buffered tasks

Analogous to
[`compute_cv_estimators()`](https://janlinnenbrink.github.io/PDAV/reference/compute_cv_estimators.md),
but based on externally generated buffered validation tasks, typically
from buffered leave-one-out resampling. Validation losses are computed
for the selected buffered tasks, augmented with realized task
descriptors, and summarized using unweighted, DWCV, TWCV, and optionally
IWCV estimators.

## Usage

``` r
compute_buffered_estimators(
  sample_dat,
  grid_dat,
  task_obj,
  model = c("rf", "lm", "ked_het_x1", "ked", "ked_het_pop"),
  response = NULL,
  fit_fun = fit_model,
  predict_fun = predict_model,
  verbose = 0,
  twcv_specs = NULL,
  predictor_vars = NULL,
  env_vars = NULL,
  iwcv_vars = NULL,
  run_iwcv = FALSE,
  iwcv_shrink_lambda = 0,
  ...
)
```

## Arguments

- sample_dat:

  Data frame of sampled observations used for validation.

- grid_dat:

  Data frame representing deployment or prediction locations.

- task_obj:

  Buffered task object, typically created by
  [`generate_buffered_loo_tasks()`](https://janlinnenbrink.github.io/PDAV/reference/generate_buffered_loo_tasks.md).

- model:

  Character string identifying the prediction model.

- response:

  Optional response variable name. If `NULL`, the function tries `z` and
  then `outcome`.

- fit_fun:

  Model-fitting function passed to `compute_buffered_task_losses()`. It
  must accept at least `train_dat`, `model`, and `response`.

- predict_fun:

  Prediction function passed to `compute_buffered_task_losses()`. It
  must accept a fitted model object and `newdata`.

- verbose:

  Verbosity level.

- twcv_specs:

  Optional named list of TWCV specifications.

- predictor_vars:

  Optional character vector of predictor variables used by the
  predictive model.

- env_vars:

  Optional character vector of environmental variables used as task
  descriptors. Defaults to `predictor_vars`.

- iwcv_vars:

  Optional character vector of variables for IWCV density-ratio
  estimation. Prediction distance `d` is appended internally when IWCV
  is enabled.

- run_iwcv:

  Logical; if `TRUE`, also compute IWCV.

- iwcv_shrink_lambda:

  Shrinkage parameter for IWCV weights.

- ...:

  Additional arguments passed to `fit_fun`.

## Value

A list with elements:

- losses:

  Buffered validation-loss data frame augmented with task descriptors.

- sample_tasks_bal:

  Currently `NULL`; retained for backward compatibility.

- grid_tasks_bal:

  Currently `NULL`; retained for backward compatibility.

- estimators:

  Named list of performance summaries.

- weights:

  Named list of weight objects for DWCV, TWCV, and optionally IWCV.

- twcv_specs:

  The TWCV specification set actually used.

## Details

Model fitting and prediction are supplied explicitly through `fit_fun`
and `predict_fun`, allowing the shared buffered-task estimator engine to
be reused across case studies without redefining global adapter
functions.

## See also

[`compute_cv_estimators()`](https://janlinnenbrink.github.io/PDAV/reference/compute_cv_estimators.md),
[`generate_buffered_loo_tasks()`](https://janlinnenbrink.github.io/PDAV/reference/generate_buffered_loo_tasks.md)
