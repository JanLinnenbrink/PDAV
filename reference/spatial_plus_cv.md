# Spatially-Aware Cross-Validation (SP-CV) Split

Performs spatially-aware cross-validation (SP-CV) by splitting samples
into folds. Used by the da_cv function. Supports three methods:

- **"SP1"**: Stage 1 spatial CV (hierarchical clustering, clusters
  assigned randomly to folds).

- **"SP2"**: Each point treated as a cluster (finest granularity).

- **"SP"**: General SP-CV with hierarchical clustering, k-means/k-modes,
  and ensemble majority voting over coordinates, environment, and target
  values.

## Usage

``` r
spatial_plus_cv(
  samples,
  response_name,
  cate_col_start = 0,
  cate_col_end = 0,
  k = 5,
  sp_threshold = 1,
  method = c("SP", "SP1", "SP2")
)
```

## Arguments

- samples:

  An `sf` object containing point geometries and a unique ID column.

- response_name:

  A data frame or matrix with target variable(s) (first column must
  match point IDs, last column treated as the target).

- cate_col_start:

  Integer, index of first categorical column in `env` (1-based).

- cate_col_end:

  Integer, index of last categorical column in `env` (1-based).

- k:

  Integer, number of folds/clusters.

- sp_threshold:

  Numeric, spatial distance threshold for hierarchical clustering
  (default 1). If set to 0, each point is its own cluster (SP2).

- method:

  Character, one of `"SP"`, `"SP1"`, or `"SP2"`.

## Value

A data frame with columns:

- `ID` - Original point ID

- `fold` - Assigned fold label

- `cluster` - Assigned cluster ID

## Examples

``` r
if (FALSE) { # \dontrun{
library(sf)
pts <- st_as_sf(data.frame(ID = 1:10, x = runif(10), y = runif(10)),
                coords = c("x", "y"), crs = 4326)
env <- data.frame(ID = 1:10, var1 = rnorm(10), var2 = runif(10), target = runif(10))
pts <- cbind(pts, env)
response_name <- "target"

folds <- spatial_plus_cv(samples = pts, response_name = response_name,
                        cate_col_start = 0, cate_col_end = 0,
                        k = 3, sp_threshold = 1, method = "SP")
} # }
```
