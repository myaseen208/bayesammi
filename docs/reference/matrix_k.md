# Orthonormal contrast matrix

`matrix_k()` returns the \\(n - 1) \times n\\ matrix \\\mathbf{K}\_n\\
whose rows are orthonormal and orthogonal to \\\mathbf{1}\_n\\ (Helmert
contrasts).
[`bayes_ammi()`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md)
writes the main effects as \\\boldsymbol{\tau} =
\mathbf{K}\_G^\top\mathbf{t}\\, which keeps them summing to zero.

## Usage

``` r
matrix_k(n)

# Default S3 method
matrix_k(n)
```

## Arguments

- n:

  number of columns (at least 2).

## Value

The \\(n - 1) \times n\\ matrix \\\mathbf{K}\_n\\.

## Details

Row \\r\\ (\\r = 1, \ldots, n - 1\\) has \\0\\ in its first \\r - 1\\
positions, \\(n - r)/c_r\\ in position \\r\\ and \\-1/c_r\\ after it,
with \\c_r = \sqrt{(n - r)(n - r + 1)}\\, so that
\$\$\mathbf{K}\_n\mathbf{K}\_n^\top = \mathbf{I}\_{n-1}, \qquad
\mathbf{K}\_n\mathbf{1}\_n = \mathbf{0}, \qquad
\mathbf{K}\_n^\top\mathbf{K}\_n = \mathbf{I}\_n -
\frac{1}{n}\mathbf{J}\_n.\$\$

## References

Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K.,
Liu, G. and Cornelius, P. L. (2011). Bayesian estimation of the additive
main effects and multiplicative interaction model. \emphCrop Science,
\strong51(4), 1458–1469. \doi10.2135/cropsci2010.06.0343

## See also

[`orthnorm`](https://myaseen208.com/bayesammi/reference/orthnorm.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
K <-
  matrix_k(
    n = 4
  )
K
#>           [,1]       [,2]       [,3]       [,4]
#> [1,] 0.8660254 -0.2886751 -0.2886751 -0.2886751
#> [2,] 0.0000000  0.8164966 -0.4082483 -0.4082483
#> [3,] 0.0000000  0.0000000  0.7071068 -0.7071068
round(x = K %*% t(K), digits = 10)
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
```
