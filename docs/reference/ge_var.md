# Genotype-by-environment variances

`ge_var()` computes the variance of the replicates of every
genotype-by-environment cell.

## Usage

``` r
ge_var(.data, .y, .gen, .env)

# Default S3 method
ge_var(.data, .y, .gen, .env)
```

## Arguments

- .data:

  a data frame with one row per plot.

- .y, .gen, .env:

  the columns of the response, genotypes and environments: bare names,
  strings or vectors.

## Value

A list (class `c("bam_ge_var", "bam_result", "list")`) with
`ge_variances`, the \\G \times E\\ matrix of variances.

## Details

\$\$s^2\_{ij} = \frac{1}{n\_{ij} - 1}\sum_l (y\_{ijl} -
\bar{y}\_{ij})^2,\$\$ `NA` for a cell with one replicate. The
within-cell sum of squares of the Bayesian AMMI model is
\\\sum\_{ij}(n\_{ij} - 1)s^2\_{ij}\\.

**Change from the earlier versions.**
[`bayes_ammi()`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md)
used the sum of these variances in place of the within-cell sum of
squares.

## References

Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K.,
Liu, G. and Cornelius, P. L. (2011). Bayesian estimation of the additive
main effects and multiplicative interaction model. \emphCrop Science,
\strong51(4), 1458–1469. \doi10.2135/cropsci2010.06.0343

## See also

[`ge_mean`](https://myaseen208.com/bayesammi/reference/ge_mean.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(Maiz)
ge_var(
    .data = Maiz
  , .y    = y
  , .gen  = entry
  , .env  = site
  )
#> <bam_ge_var> Genotype-by-environment variances
#> 
#> $ge_variances
#>         1      2       3      4       5       6      7       8       9     10
#> 1  232591 239434   69226 295490  234259  233029 176896  751234 3847788  28327
#> 2 1515737 651925  289329 440025  733619  112847 550173  505106  910349  62165
#> 3 1546435 560798  763984 120665  261589  167639 903642 1394066  207997 124323
#> 4  211897 426378 2003423 160082 1523162   55110 155237  767942 3042002 176054
#> 5  488114 386601  815919 317691  241103 1412073 572317  348099  514745 159698
#> 6 1083169 126145 3318749 304564  515868  500458  19295  185805  260306  88058
#> 7 1413644  84034 1117405 146933  142661  233631 404165  748141 4224452  19173
#> 8   55005  51764  425938 422738   82656  331409 336668 1724361 2084463  21066
#> 9 1739201 184724 1606346  61565   44665  820531 504980   75568 1655346  47196
#>        11      12      13      14      15     16      17      18      19     20
#> 1  411126   75804   87124 1151139  559363 263864 1449535 1095435 1802424 576507
#> 2  350665 1022490  236509  606342  350565  89720  356624 1965285  788734  92785
#> 3 1274303  726702  107722 1401692  344613 168845  367372 3108304 2035153 335183
#> 4 1091659  611209 2781308 1293484  168213 428856   95103  746374 3128493 406134
#> 5  206193  296823   33108  304752  727480 659732  687992 2180437 2238274 123916
#> 6 1009733 1724462  710806 2365610  726838 295147  615533 1062673 3818259 349410
#> 7  389083  351349 1190591 1784078  340745  22361  284678  291728 1069468 470154
#> 8 1982295  182001 1044828   46565  150887 278821  278014  908955 1511362 343330
#> 9 1820457  372354  436979 1921407 1325558 418682  809679  217594 2354966 448989
```
