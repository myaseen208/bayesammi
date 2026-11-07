# Environment effects

`e_eff()` computes the environment main effects, the deviations of the
environment means from the grand mean.

## Usage

``` r
e_eff(.data, .y, .gen, .env)

# Default S3 method
e_eff(.data, .y, .gen, .env)
```

## Arguments

- .data:

  a data frame with one row per plot.

- .y, .gen, .env:

  the columns of the response, genotypes and environments: bare names,
  strings or vectors.

## Value

An \\E \times 1\\ matrix with column `EEffs` (row names: environments).

## Details

\$\$\hat{\delta}\_j = \bar{y}\_{\cdot j} - \bar{y}, \qquad
\bar{y}\_{\cdot j} = \frac{1}{G}\sum_i \bar{y}\_{ij}.\$\$

## References

Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K.,
Liu, G. and Cornelius, P. L. (2011). Bayesian estimation of the additive
main effects and multiplicative interaction model. \emphCrop Science,
\strong51(4), 1458–1469. \doi10.2135/cropsci2010.06.0343

## See also

[`g_eff`](https://myaseen208.com/bayesammi/reference/g_eff.md),
[`ge_eff`](https://myaseen208.com/bayesammi/reference/ge_eff.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(Maiz)
e_eff(
    .data = Maiz
  , .y    = y
  , .gen  = entry
  , .env  = site
  )
#>          EEffs
#> 1  -1245.60722
#> 2   -648.69333
#> 3    246.78167
#> 4    371.85944
#> 5     97.09278
#> 6   1447.16778
#> 7  -1628.79611
#> 8   -830.45444
#> 9    112.50667
#> 10 -1921.71278
#> 11   449.05389
#> 12  2658.36222
#> 13  1473.96222
#> 14  1194.15944
#> 15   190.89000
#> 16   548.02611
#> 17    18.55111
#> 18  -306.59056
#> 19 -2206.03500
#> 20   -20.52389
```
