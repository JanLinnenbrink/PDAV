# Conditional logging with verbosity levels

Prints a timestamped message when the requested verbosity threshold is
met.

## Usage

``` r
log_message(verbose = 0, level = 1, ...)
```

## Arguments

- verbose:

  Integer verbosity level supplied by the caller.

- level:

  Minimum verbosity level required for printing.

- ...:

  Objects to be concatenated into the message text.

## Value

Invisibly `NULL`.
