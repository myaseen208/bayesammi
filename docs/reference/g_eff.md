# Genotype effects

`g_eff()` computes the genotype main effects, the deviations of the
genotype means from the grand mean.

## Usage

``` r
g_eff(.data, .y, .gen, .env)

# Default S3 method
g_eff(.data, .y, .gen, .env)
```

## Arguments

- .data:

  a data frame with one row per plot.

- .y, .gen, .env:

  the columns of the response, genotypes and environments: bare names,
  strings or vectors.

## Value

A \\G \times 1\\ matrix with column `GEffs` (row names: genotypes).

## Details

\$\$\hat{\tau}\_i = \bar{y}\_{i\cdot} - \bar{y}, \qquad
\bar{y}\_{i\cdot} = \frac{1}{E}\sum_j \bar{y}\_{ij},\$\$ the
least-squares estimates of the AMMI model (equal to the deviations of
the plot means when the data are balanced).

## References

Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K.,
Liu, G. and Cornelius, P. L. (2011). Bayesian estimation of the additive
main effects and multiplicative interaction model. \emphCrop Science,
\strong51(4), 1458–1469. \doi10.2135/cropsci2010.06.0343

## See also

[`e_eff`](https://myaseen208.com/bayesammi/reference/e_eff.md),
[`ge_eff`](https://myaseen208.com/bayesammi/reference/ge_eff.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(Maiz)
g_eff(
    .data = Maiz
  , .y    = y
  , .gen  = entry
  , .env  = site
  )
#>        GEffs
#> 1 -241.10903
#> 2 -254.87903
#> 3  -35.22403
#> 4  359.71097
#> 5  388.73597
#> 6  472.03347
#> 7 -185.82153
#> 8 -574.94278
#> 9   71.49597
```
