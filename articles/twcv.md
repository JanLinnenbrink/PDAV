# TWCV

## Introduction

Target-Weighted Cross-Validation (TWCV) was developed by Brenning and
Suesse (2026). It uses iterative calibration weighting to align the
validation task (predictor values of the training points and
nearest-neighbour distances (NNDs) between folds) to the prediction task
(predictor values of the prediction points and NNDs between prediction
and training points).

The more detailed workflow is:

1.  Compute CV losses by training a model on the training folds fold and
    predicting to the held-out fold
2.  Compute the CV validation task and the prediction task (NNDs and
    predictor values)
3.  Calculate quantiles of the balancing variables for training and
    prediction points separately (NNDs and predictors)
4.  Calculate the proportion of training points and prediction points
    falling in each quantile
5.  Applies iterative proportional fitting (“raking”) to align the
    proportions of the training points to those in the prediction
    points:
    1.  For each balancing variable, calculate the current number of
        training points in each quantile.
    2.  Calculate the target number of training points in each quantile
        based on the target proportions.
    3.  Calculate multiplicative adjustment factors that transform the
        current counts to the target counts.
    4.  Because the weighting is done for each variable seperately
        (“marginal”), correcting one variable can disturb variables
        adjusted earlier.
    5.  Repeat this process until the maximum relative change in weights
        between successive iterations falls below a threshold.
6.  Normalize weights and shrink them towards 1 to mitigate extreme
    values.

## Setup

``` r

library(PDAV)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(caret)
#> Loading required package: ggplot2
#> Loading required package: lattice
library(terra)
#> terra 1.9.27
library(sf)
#> Linking to GEOS 3.12.1, GDAL 3.8.4, PROJ 9.4.0; sf_use_s2() is TRUE
library(simsam)
library(ggplot2)
library(cowplot)
library(tidyterra)
#> 
#> Attaching package: 'tidyterra'
#> The following object is masked from 'package:stats':
#> 
#>     filter
library(ggnewscale)

set.seed(100)
```

## Simulate predictors and response

``` r

r <- PDAV:::generate_rast()
predictor_stack <- r[[setdiff(names(r), "outcome")]]
cate_rasters <- which(names(r) %in% c("forest", "grass"))

n_sample <- 100

sampling_r <- r
sampling_r[sampling_r$elev > 1] <- NA

samples <- sam_field(
    x = sampling_r,
    size = n_sample,
    method = sample_clustered(nclusters = 10, radius = 30, na.rm = TRUE)
)

predictor_vars <- names(predictor_stack)
env_vars <- predictor_vars
task_vars <- env_vars
response <- "outcome"
pred_points <- spatSample(predictor_stack, 1000, method = "regular", na.rm = TRUE, as.points = TRUE) |>
    st_as_sf()

sample_coords <- st_coordinates(samples) |> as.data.frame() |> rename("x" = X, "y" = Y)
grid_coords <- st_coordinates(pred_points) |> as.data.frame() |> rename("x" = X, "y" = Y)
sample_dat <- st_drop_geometry(samples) |>
    mutate(id = row_number()) |>
    cbind(sample_coords)
grid_dat <- st_drop_geometry(pred_points) |>
    mutate(id = row_number()) |>
    cbind(grid_coords)

twcv_specs <- PDAV:::make_twcv_specs(
    predictor_vars = predictor_vars,
    include_distance = TRUE,
    balance_by = 0.2,
    shrink_lambda = 0.2,
    name = "twcv_extended"
)

fit_fun <- function(train_dat, model, response, predictor_vars = NULL, ...) {
    if (is.null(predictor_vars)) {
        predictor_vars <- setdiff(names(train_dat), response)
    }

    dat <- train_dat[, c(response, predictor_vars), drop = FALSE]

    ranger::ranger(
        formula = stats::reformulate(predictor_vars, response = response),
        data = dat,
        ...
    )
}

predict_fun <- function(fit, newdata, ...) {
    predict(fit, data = newdata)$predictions
}

model = "rf"

folds <- sample(rep(1:5, length.out = n_sample))

rel_columns <- c("fold", "obs", "pred", "error", "elev")
```

![Figure 1: The predictor
stack](twcv_files/figure-html/unnamed-chunk-3-1.png)

Figure 1: The predictor stack

## Run TWCV

### 1. Compute CV losses

Iterate over each fold and calculate the predictions and the resulting
pointwise CV error:

``` r

pred <- PDAV:::compute_cv_predictions(
    sample_dat = sample_dat,
    folds = folds,
    model = model,
    response = response,
    fit_fun = fit_fun,
    predict_fun = predict_fun
)

pe <- PDAV:::compute_pointwise_errors(
    obs = sample_dat[[response]],
    pred = pred
)

cv_losses <- data.frame(
    id = sample_dat$id,
    fold = folds,
    obs = pe$obs,
    pred = pe$pred,
    error = pe$error,
    se = pe$se,
    ae = pe$ae
)

knitr::kable(head(cv_losses))
```

|  id | fold |        obs |       pred |      error |        se |        ae |
|----:|-----:|-----------:|-----------:|-----------:|----------:|----------:|
|   1 |    4 |  0.0250339 |  0.0477339 | -0.0226999 | 0.0005153 | 0.0226999 |
|   2 |    4 |  0.0194739 |  0.0017141 |  0.0177599 | 0.0003154 | 0.0177599 |
|   3 |    2 |  0.0295682 |  0.0369364 | -0.0073683 | 0.0000543 | 0.0073683 |
|   4 |    3 | -0.1403769 | -0.0757283 | -0.0646486 | 0.0041794 | 0.0646486 |
|   5 |    1 | -0.1468337 | -0.0656559 | -0.0811778 | 0.0065898 | 0.0811778 |
|   6 |    2 | -0.0060736 |  0.0310038 | -0.0370774 | 0.0013747 | 0.0370774 |

### 2. Compute the CV validation task and the prediction task

``` r

d_sample <- PDAV:::nearest_neighbor_distance(
    query_coords = sample_dat[, c("x", "y")],
    ref_coords = sample_dat[, c("x", "y")],
    exclude_self = TRUE
)

d_grid <- PDAV:::nearest_neighbor_distance(
    query_coords = grid_dat[, c("x", "y")],
    ref_coords = sample_dat[, c("x", "y")],
    exclude_self = FALSE
)

sample_tasks <- data.frame(
    id = sample_dat$id,
    d = d_sample,
    stringsAsFactors = FALSE
)

grid_tasks <- data.frame(
    id = grid_dat$id,
    d = d_grid,
    stringsAsFactors = FALSE
)

for (v in env_vars) {
    sample_tasks[[v]] <- sample_dat[[v]]
    grid_tasks[[v]] <- grid_dat[[v]]
}

tasks <- list(
    sample_tasks = sample_tasks,
    grid_tasks = grid_tasks
)

sample_desc <- tasks$sample_tasks
idx <- match(cv_losses$id, sample_desc$id)

# combine CV losses with task descriptors for training points
out <- cv_losses
for (v in task_vars) {
    if (!(v %in% names(sample_desc))) {
        stop("Variable '", v, "' missing in sample_tasks.", call. = FALSE)
    }
    out[[v]] <- sample_desc[[v]][idx]
}

# Calculates the NND between folds based on distance matrix / matrices
d_realized <- PDAV:::compute_cv_prediction_distances(
    sample_dat = sample_dat,
    folds = cv_losses$fold
)

# Update the NND attached to the training points from NND between points to NND between folds
idx_d <- match(out$id, sample_dat$id)
out$d <- d_realized[idx_d]

aug <- list(
    losses = out,
    grid_tasks = tasks$grid_tasks
)

knitr::kable(head(aug$losses))
```

| id | fold | obs | pred | error | se | ae | temp | moisture | ph | slope | solar | dist_road | prod | elev | forest | grass | d |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 4 | 0.0250339 | 0.0477339 | -0.0226999 | 0.0005153 | 0.0226999 | 0.8769668 | -1.5636978 | -1.4725255 | 0.5039787 | -1.5788051 | -0.0977249 | 1.0431762 | -0.2515965 | 1 | 1 | 26.683328 |
| 2 | 4 | 0.0194739 | 0.0017141 | 0.0177599 | 0.0003154 | 0.0177599 | -0.1627803 | -2.4732208 | -1.7179238 | -0.6153057 | -1.2322361 | 0.8162922 | 0.0095712 | -0.9679169 | 1 | 1 | 9.848858 |
| 3 | 2 | 0.0295682 | 0.0369364 | -0.0073683 | 0.0000543 | 0.0073683 | 1.1839176 | 0.0311049 | -0.5123982 | -1.1399273 | -1.2682281 | 0.5378851 | -0.2595845 | 0.4287385 | 1 | 1 | 11.661904 |
| 4 | 3 | -0.1403769 | -0.0757283 | -0.0646486 | 0.0041794 | 0.0646486 | -0.1536180 | -1.4972850 | -0.9786745 | -0.2723358 | -0.9876763 | 0.9412718 | 0.7548076 | -0.3165647 | 0 | 1 | 16.552945 |
| 5 | 1 | -0.1468337 | -0.0656559 | -0.0811778 | 0.0065898 | 0.0811778 | -0.1915660 | -1.4061724 | -1.3192238 | 0.6799157 | -0.7441159 | 0.5422058 | -0.1562645 | -0.9152544 | 0 | 1 | 8.602325 |
| 6 | 2 | -0.0060736 | 0.0310038 | -0.0370774 | 0.0013747 | 0.0370774 | -0.5612379 | -1.7635113 | -1.2023100 | -0.8075206 | -1.2668175 | -0.1518995 | -0.2905415 | 0.1555381 | 1 | 1 | 10.198039 |

``` r

knitr::kable(head(aug$grid_tasks))
```

| id | d | temp | moisture | ph | slope | solar | dist_road | prod | elev | forest | grass |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 7.615773 | -0.7063682 | 1.6181636 | 0.3247341 | -0.1572279 | -1.5364642 | -0.656340 | -0.0769430 | -1.9746046 | 1 | 1 |
| 2 | 4.472136 | -0.2595710 | 0.5600811 | 0.1500221 | -0.4940667 | -1.0587877 | -1.305025 | 0.3643323 | -1.5782595 | 1 | 1 |
| 3 | 2.828427 | -0.2050061 | 0.0364014 | 0.0129812 | -0.3565046 | -0.3649710 | -1.736716 | 0.8307905 | -1.1459060 | 1 | 1 |
| 4 | 8.246211 | 0.1409018 | -1.1082273 | 0.0015563 | -0.3931227 | 0.0028835 | -1.469470 | 0.8723531 | -0.8795467 | 1 | 1 |
| 5 | 10.049876 | 0.3179103 | -0.5781679 | -0.6160732 | -0.2139655 | 0.6954836 | -1.924319 | 0.3423136 | -1.5066579 | 1 | 1 |
| 6 | 12.206556 | 0.5716950 | -0.6276057 | -0.2341143 | 0.0015311 | 0.1396628 | -1.451022 | 0.5832729 | -1.1956650 | 1 | 1 |

### 3. Calculate quantiles of the balancing variables

``` r

cv_losses <- aug$losses
grid_tasks <- aug$grid_tasks

unweighted_losses <- list(
    unweighted = PDAV:::summarize_losses(cv_losses)
)

sample_tasks_bal <- PDAV:::prepare_for_balancing(
    df = cv_losses,
    vars = twcv_specs$twcv_extended$balancing_vars,
    ref_df = grid_tasks,
    by = twcv_specs$twcv_extended$balance_by
)

grid_tasks_bal <- PDAV:::prepare_for_balancing(
    df = grid_tasks,
    vars = twcv_specs$twcv_extended$balancing_vars,
    ref_df = grid_tasks,
    by = twcv_specs$twcv_extended$balance_by
)

bal <- list(
    sample_tasks_bal = sample_tasks_bal,
    grid_tasks_bal = grid_tasks_bal
)
```

![](twcv_files/figure-html/plot-1.png)

### 4. Calculate frequencies for the quantiles

``` r

balance_df <- as.data.frame(
    lapply(twcv_specs$twcv_extended$balancing_vars, function(v) sample_tasks_bal[[paste0(v, "_cat")]])
)
names(balance_df) <- twcv_specs$twcv_extended$balancing_vars

# Calculates proportion of predpoints in each quantile of each predictor used for weighting
target_margins <- PDAV:::compute_target_margins_generic(
    grid_tasks_bal = grid_tasks_bal,
    balancing_vars = twcv_specs$twcv_extended$balancing_vars
)
target_margins$elev
#> [1] 0.2001953 0.2001953 0.1992188 0.2001953 0.2001953

# For visualization only: margins of samples
levs <- seq_along(target_margins[["elev"]])
freq_samples <- table(factor(as.integer(balance_df[["elev"]]), levels = levs))
as.numeric(freq_samples) / sum(freq_samples)
#> [1] 0.40 0.16 0.24 0.19 0.01
```

### 5. Apply Raking

In this example, one elevation quintile of the prediction locations is
not covered by the training points. However, this does not return an
error. Instead, if a quintile is not covered by the training points, the
weights shrink towards 0 and the algorithm “converges” A better approach
would likely be to stop the algorithm and return an error message
hinting towards avoiding this extreme extrapolation.

``` r

margin_names <- names(target_margins)
max_iter <- 500
tol = 1e-6
n <- nrow(balance_df)
w <- rep(1, n)

for (iter in seq_len(max_iter)) {
    w_old <- w

    for (m in margin_names) {
        x <- as.integer(balance_df[[m]])
        levs <- seq_along(target_margins[[m]])
        target_prop <- target_margins[[m]]

        # calculate the number of training points currently in each quantile
        # uses the weights from the previous balancing variable
        # -> already weighted, but likely not ideally for this variable
        current_totals <- tapply(w, factor(x, levels = levs), sum)
        current_totals[is.na(current_totals)] <- 0

        # calculate the desired number of training points in each quantile that matches the target margins:
        # Number of training points * target proportion vector
        target_totals <- sum(w) * target_prop

        adj <- rep(1, length(levs))
        ok <- current_totals > 0

        # calculates the weight needed to adjust the training point distribution to the target margins
        # (weights > 1 are used to up-weigh underrepresented classes, weights < 1 to down-weight over-represented ones)
        adj[ok] <- target_totals[ok] / current_totals[ok]
        adj[!ok] <- NA_real_

        # Checks for all training points if they fall in the distribution of the grid
        # However, if one quintile of the prediction domain is not covered by the training data, this does not return an error
        # Instead, if a quintile is not covered by the training points, the weights shrink towards 0 and the algorithm "converges"
        valid <- !is.na(adj[x])
        w[valid] <- w[valid] * adj[x[valid]]
    }

    # Compute the relative strength of the absolute change of weights from base (or previous) to new weights
    rel_change <- max(abs(w - w_old) / pmax(abs(w_old), 1e-12))

    # When the changes converge (i.e., when the weight difference from previous iteration to new one is small), stop and return weights
    # (when the weights for predictor B are changed, weights for predictor A might be off again, and another iteration starts,
    # until the diff between them is small)
    if (rel_change < tol) break
}

converged <- iter < max_iter || rel_change < tol

tw <- list(
    weights = w,
    converged = converged,
    iterations = iter
)
tw$weights
```

### 6. Normalize and shrink weights

``` r

# Weights are normalized by their mean and shrinked towards 1 to mitigate extreme values
shrink_lambda <- 0
tw$weights_raw <- PDAV:::normalize_weights(tw$weights)
tw$weights <- PDAV:::shrink_weights(tw$weights_raw, lambda = shrink_lambda)
tw$shrink_lambda <- twcv_specs$twcv_extended$shrink_lambda
tw$balancing_vars <- twcv_specs$twcv_extended$balancing_vars

tw$weights
#>   [1] 2.503053e+01 6.715845e-47 1.670722e-23 3.589109e-19 8.420921e-28
#>   [6] 1.456392e+00 2.555495e-27 1.613776e-19 1.270596e-47 2.552497e-10
#>  [11] 1.014327e+00 1.607507e-46 4.583316e-66 1.492707e+00 7.007203e-39
#>  [16] 4.460577e+00 4.407537e-11 1.737345e-42 7.184545e-20 3.314342e-47
#>  [21] 4.231939e+00 3.475718e-31 3.636416e-01 8.928592e-01 1.699966e-48
#>  [26] 4.028438e-26 1.127379e-63 7.457782e-16 3.151407e-52 1.292677e+01
#>  [31] 2.191856e+00 1.398897e-24 2.725420e-33 3.068857e-12 2.220361e+00
#>  [36] 3.695944e-07 1.193136e-25 2.804604e+00 5.095991e-11 4.686885e-26
#>  [41] 5.729552e-19 4.286948e+00 3.858182e-01 1.580222e-46 1.278465e-30
#>  [46] 9.483699e-72 2.115043e-18 1.117528e-02 1.475736e-24 7.675915e-15
#>  [51] 1.795980e-23 3.787657e-16 8.930879e-02 3.509206e-05 6.501129e-25
#>  [56] 7.897270e-02 1.873253e-36 1.370882e+00 2.521778e-27 1.009272e-28
#>  [61] 7.405940e-18 7.134251e+00 2.408758e-11 1.474823e+00 7.202142e-12
#>  [66] 2.293770e-40 6.370046e-31 3.358196e-68 3.221041e-45 7.018954e-40
#>  [71] 4.272597e-35 2.921522e-02 2.686407e-31 4.280209e-38 3.949011e+00
#>  [76] 8.289632e-10 2.035687e-25 2.570220e-03 3.210328e+00 1.176227e-18
#>  [81] 4.236346e-13 2.270758e-38 1.718837e-12 2.562944e+00 8.168795e-23
#>  [86] 6.185815e-01 5.422943e-11 3.058058e-16 4.881522e-12 3.459297e-37
#>  [91] 6.569227e-27 8.606627e-06 1.359560e-18 2.277638e-05 3.458782e+00
#>  [96] 6.261614e-10 2.596027e+00 7.705747e-28 3.425543e+00 6.228199e+00
```

### 7. Return weighted error estimates

``` r

est_list <- PDAV:::summarize_losses(cv_losses, tw$weights)
weight_objects <- tw

result <- list(
    losses = cv_losses,
    estimators = est_list,
    weights = weight_objects,
    twcv_specs = twcv_specs
)

# Unweighted error:
unweighted_losses
#> $unweighted
#>         bias          mse         rmse          mae 
#> 0.0009066178 0.0032602974 0.0570990138 0.0432061638

# Weighted error:
result$estimators
#>         bias          mse         rmse          mae 
#> -0.005962343  0.003049602  0.055223204  0.042842620
```

Brenning, Alexander, and Thomas Suesse. 2026. *Aligning Validation with
Deployment: Target-Weighted Cross-Validation for Spatial Prediction*.
arXiv. <https://doi.org/10.48550/ARXIV.2603.29981>.
