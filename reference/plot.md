# Plot PDAV classes

Generic plot function for prediction-domain adaptive validation results
Classes

## Usage

``` r
# S3 method for class 'nndm'
plot(x, type = "strict", stat = "ecdf", ...)

# S3 method for class 'knndm'
plot(x, type = "strict", stat = "ecdf", ...)

# S3 method for class 'da_cv'
plot(x, ...)
```

## Arguments

- x:

  An object of type *da_cv*.

- type:

  String, defaults to "strict" to show the original nearest neighbour
  distance definitions in the legend. Alternatively, set to "simple" to
  have more intuitive labels.

- stat:

  String, defaults to "ecdf" but can be set to "density" to estimate
  density functions.

- ...:

  other arguments.

## Author

Carles Milà

Jan Linnenbrink
