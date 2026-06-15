# Dissimilarity-Adaptive Cross-Validation (DA-CV)

This function implements dissimilarity-adaptive cross-validation
(DA-CV), which combines random CV (RDM-CV) and spatial CV (SP-CV) based
on adversarial validation (AV). The method follows the framework
described in [Wang et al.,
2025](https://doi.org/10.1016/j.ecoinf.2025.103287).

## Usage

``` r
da_cv(
  samples,
  predictors,
  response,
  folds_k = 5,
  autoc_threshold = 0.2,
  seed = NULL,
  ...
)
```

## Arguments

- samples:

  An `sf` object with point samples. Must include the response variable
  and predictor variables as attributes.

- predictors:

  A
  [`terra::rast`](https://rspatial.github.io/terra/reference/rast.html)
  object containing predictor rasters aligned to the study area.

- response:

  Name of the response column in `samples`.

- folds_k:

  Integer. Number of folds for cross-validation.

- autoc_threshold:

  Numeric. Spatial autocorrelation threshold used in `spatial_plus_cv`.

- seed:

  Optional integer. Random seed for reproducibility (passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html)).

- ...:

  Parameters passed to
  [`spatial_plus_cv()`](https://janlinnenbrink.github.io/PDAV/reference/spatial_plus_cv.md).
  Most importantly `method`.

## Value

A list with components:

- dissimilarity:

  Dissimilarity index D; ranging from 0 to 1

- threshold:

  Threshold T(D) = D \* 0.5

- similarity_raster:

  Similarity raster

- category_raster:

  Binarized similarity raster. 1 - Different, 2 - Similar

- weights:

  Relative weights of "similar" and "different" areas

- folds_RDM:

  folds from RDM-CV

- folds_SP:

  folds from SP-CV

## Details

Steps:

1.  Generate adversarial validation dataset by combining sample points
    with randomly sampled prediction locations from `predictors`.

2.  Train a classifier (random forest) to distinguish samples from
    prediction locations.

3.  Evaluate classifier via AUC → dissimilarity index D.

4.  Apply classifier to all raster cells → similarity map.

5.  Threshold map at T(D) = D \* 0.5 to classify cells into {similar,
    different}.

6.  Calculate random and spatial+ cross-validation splits

7.  Combine RMSEs via weighted average by area proportion. (Needs to be
    done manually after running `da_cv()`)

## Examples

``` r
if (FALSE) { # \dontrun{
library(sf)
library(terra)

# Sample points
pts <- st_as_sf(data.frame(
  x = runif(50, 0, 100),
  y = runif(50, 0, 100),
  yvar = rnorm(50)
), coords = c("x", "y"), crs = 32633)

# Predictor rasters
r <- rast(ncols=100, nrows=100, xmin=0, xmax=100, ymin=0, ymax=100)
values(r) <- runif(ncell(r))

# Run DA-CV (RDM_CV and spatial_plus_cv (must be available)
result <- da_cv(
  samples = pts,
  predictors = r,
  response = "yvar",
  folds_k = 5,
  autoc_threshold = 0.2
)
print(result$weights)
} # }
```
