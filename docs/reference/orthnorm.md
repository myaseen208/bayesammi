# Gram-Schmidt orthonormalisation

`orthnorm()` orthogonalises the columns of a matrix by the Gram-Schmidt
process, optionally completes them to a basis of the whole space and
scales them to unit length.

## Usage

``` r
orthnorm(u = NULL, basis = TRUE, norm = TRUE)

# Default S3 method
orthnorm(u = NULL, basis = TRUE, norm = TRUE)
```

## Arguments

- u:

  a matrix (or vector) whose columns are orthogonalised.

- basis:

  complete the columns to a basis of the whole space (default `TRUE`).

- norm:

  scale the columns to unit length (default `TRUE`).

## Value

The matrix of orthogonal (orthonormal with `norm = TRUE`) columns;
`NULL` for `u = NULL`.

## Details

For columns \\\mathbf{u}\_1, \ldots, \mathbf{u}\_m\\, \$\$\mathbf{v}\_1
= \mathbf{u}\_1, \qquad \mathbf{v}\_i = \mathbf{u}\_i - \sum\_{h \<
i}\frac{\mathbf{u}\_i^\top
\mathbf{v}\_h}{\mathbf{v}\_h^\top\mathbf{v}\_h}\mathbf{v}\_h,\$\$ each
divided by its length when `norm = TRUE`. With `basis = TRUE` and fewer
columns than rows, unit vectors are added (those least explained by the
given columns) to complete a basis.

## References

Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K.,
Liu, G. and Cornelius, P. L. (2011). Bayesian estimation of the additive
main effects and multiplicative interaction model. \emphCrop Science,
\strong51(4), 1458–1469. \doi10.2135/cropsci2010.06.0343

## See also

[`matrix_k`](https://myaseen208.com/bayesammi/reference/matrix_k.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
V <-
  orthnorm(
      u     = cbind(1, 1:4)
    , basis = TRUE
    , norm  = TRUE
    )
round(x = crossprod(V), digits = 10)
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    0    0    0
#> [2,]    0    1    0    0
#> [3,]    0    0    1    0
#> [4,]    0    0    0    1
```
