# Published values used by the tests and vignettes

Reads `inst/extdata/reference.csv`: one row per published number, with
its source, the R call that recomputes it with this package (`call`),
the expression that extracts the number from the result (`extract`, in
terms of `x`), the printed value and the digits to compare. The same
file drives the test of published values and the vignette that
reproduces the published examples. Each distinct `call` is evaluated
once.

## Usage

``` r
bam_reference(source = NULL, evaluate = FALSE, level = 1)
```

## Arguments

- source:

  optional regular expression; keep rows whose `source` matches.

- evaluate:

  logical; evaluate the rows and add columns `computed` and `agrees`
  (default `FALSE`).

- level:

  largest `level` evaluated (1 = fast rows, 2 = slow rows; default 1);
  rows above it get `NA`.

## Value

A `data.table` with columns `id`, `source`, `quantity`, `call`,
`extract`, `value`, `digits`, `level`, `note` (and `computed`,
`agrees`).

## See also

Other bayesammi utilities:
[`bam_notation()`](https://myaseen208.com/bayesammi/reference/bam_notation.md),
[`bam_save()`](https://myaseen208.com/bayesammi/reference/bam_save.md),
[`bam_theme()`](https://myaseen208.com/bayesammi/reference/bam_theme.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
ref <-
  bam_reference(
      source   = NULL
    , evaluate = FALSE
    , level    = 1
    )
ref
#>                       id                                         source
#>                   <char>                                         <char>
#> 1: crossa2011_ls_lambda1 Crossa et al. (2011), Crop Science 51, p. 1463
#> 2: crossa2011_ls_lambda2 Crossa et al. (2011), Crop Science 51, p. 1463
#> 3:       crossa2011_mean Crossa et al. (2011), Crop Science 51, Table 1
#>                                              quantity
#>                                                <char>
#> 1:           least-squares lambda_1 of the maize data
#> 2:           least-squares lambda_2 of the maize data
#> 3: grand mean of the maize data (posterior mean 4858)
#>                                                                    call
#>                                                                  <char>
#> 1: ge_ammi(.data = Maiz, .y = y, .gen = entry, .env = site, .rep = rep)
#> 2: ge_ammi(.data = Maiz, .y = y, .gen = entry, .env = site, .rep = rep)
#> 3:             ge_mean(.data = Maiz, .y = y, .gen = entry, .env = site)
#>         extract value digits level   note
#>          <char> <int>  <int> <int> <char>
#> 1: x$lambdas[1]  5923      0     1       
#> 2: x$lambdas[2]  3070      0     1       
#> 3: x$grand_mean  4858      0     1       
# \donttest{
bam_reference(
    source   = "."
  , evaluate = TRUE
  , level    = 1
  )
#>                       id                                         source
#>                   <char>                                         <char>
#> 1: crossa2011_ls_lambda1 Crossa et al. (2011), Crop Science 51, p. 1463
#> 2: crossa2011_ls_lambda2 Crossa et al. (2011), Crop Science 51, p. 1463
#> 3:       crossa2011_mean Crossa et al. (2011), Crop Science 51, Table 1
#>                                              quantity
#>                                                <char>
#> 1:           least-squares lambda_1 of the maize data
#> 2:           least-squares lambda_2 of the maize data
#> 3: grand mean of the maize data (posterior mean 4858)
#>                                                                    call
#>                                                                  <char>
#> 1: ge_ammi(.data = Maiz, .y = y, .gen = entry, .env = site, .rep = rep)
#> 2: ge_ammi(.data = Maiz, .y = y, .gen = entry, .env = site, .rep = rep)
#> 3:             ge_mean(.data = Maiz, .y = y, .gen = entry, .env = site)
#>         extract value digits level   note computed agrees
#>          <char> <int>  <int> <int> <char>    <num> <lgcl>
#> 1: x$lambdas[1]  5923      0     1        5922.891   TRUE
#> 2: x$lambdas[2]  3070      0     1        3070.454   TRUE
#> 3: x$grand_mean  4858      0     1        4858.138   TRUE
# }
```
