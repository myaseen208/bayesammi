# Genotype-by-environment interaction effects

`ge_eff()` computes the interaction effects (the cell means centred on
genotype and environment means) and their singular value decomposition,
the least-squares estimate of the multiplicative part of the AMMI model.

## Usage

``` r
ge_eff(.data, .y, .gen, .env)

# Default S3 method
ge_eff(.data, .y, .gen, .env)
```

## Arguments

- .data:

  a data frame with one row per plot.

- .y, .gen, .env:

  the columns of the response, genotypes and environments: bare names,
  strings or vectors.

## Value

A list (class `c("bam_ge_eff", "bam_result", "list")`) with `ge_effects`
(the \\G \times E\\ matrix \\\mathbf{Z}\\) and `ge_svd`
([`svd()`](https://rdrr.io/r/base/svd.html) of it: `d`, `u`, `v`).

## Details

\$\$\mathbf{Z} = (z\_{ij}), \quad z\_{ij} = \bar{y}\_{ij} -
\bar{y}\_{i\cdot} - \bar{y}\_{\cdot j} + \bar{y}, \qquad \mathbf{Z} =
\sum_k \lambda_k\boldsymbol{\alpha}\_k\boldsymbol{\gamma}\_k^\top,\$\$
with singular values \\\lambda_1 \ge \lambda_2 \ge \cdots\\ and
orthonormal singular vectors orthogonal to \\\mathbf{1}\\.

## References

Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K.,
Liu, G. and Cornelius, P. L. (2011). Bayesian estimation of the additive
main effects and multiplicative interaction model. \emphCrop Science,
\strong51(4), 1458–1469. \doi10.2135/cropsci2010.06.0343

Gauch, H. G. (1992). *Statistical Analysis of Regional Yield Trials:
AMMI Analysis of Factorial Designs*. Elsevier, Amsterdam.

## See also

[`ge_ammi`](https://myaseen208.com/bayesammi/reference/ge_ammi.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(Maiz)
fit <-
  ge_eff(
      .data = Maiz
    , .y    = y
    , .gen  = entry
    , .env  = site
    )
fit$ge_svd$d
#> [1] 5.922891e+03 3.070454e+03 2.552407e+03 2.320222e+03 1.758671e+03
#> [6] 1.313415e+03 7.926324e+02 7.570277e+02 3.268710e-12
```
