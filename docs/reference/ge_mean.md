# Genotype-by-environment means

`ge_mean()` computes the table of genotype-by-environment cell means and
the grand mean, the data of the AMMI model.

## Usage

``` r
ge_mean(.data, .y, .gen, .env)

# Default S3 method
ge_mean(.data, .y, .gen, .env)
```

## Arguments

- .data:

  a data frame with one row per plot.

- .y, .gen, .env:

  the columns of the response, genotypes and environments: bare names,
  strings or vectors.

## Value

A list (class `c("bam_ge_mean", "bam_result", "list")`) with `ge_means`,
the \\G \times E\\ matrix of means (row and column names: genotypes and
environments), and `grand_mean`.

## Details

With \\y\_{ijl}\\ the response of genotype \\i\\ in replicate \\l\\ of
environment \\j\\ and \\n\_{ij}\\ the number of replicates of the cell,
\$\$\bar{y}\_{ij} = \frac{1}{n\_{ij}}\sum_l y\_{ijl}, \qquad \bar{y} =
\frac{1}{GE}\sum\_{ij}\bar{y}\_{ij}.\$\$ Genotypes and environments are
in the order of their factor levels (sorted values when the columns are
not factors).

## References

Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K.,
Liu, G. and Cornelius, P. L. (2011). Bayesian estimation of the additive
main effects and multiplicative interaction model. \emphCrop Science,
\strong51(4), 1458–1469. \doi10.2135/cropsci2010.06.0343

## See also

[`ge_var`](https://myaseen208.com/bayesammi/reference/ge_var.md),
[`ge_eff`](https://myaseen208.com/bayesammi/reference/ge_eff.md),
[`ge_ammi`](https://myaseen208.com/bayesammi/reference/ge_ammi.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(Maiz)
ge_mean(
    .data = Maiz
  , .y    = y
  , .gen  = entry
  , .env  = site
  )
#> <bam_ge_mean> Genotype-by-environment means
#> 
#> $ge_means
#>      1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
#> 1 3622 3728 5554 4566 4380 6437 2832 6011 4647 3100 4433 6873 6721 5849 4601
#> 2 3426 3919 4937 4963 5201 6036 2515 5278 4714 2972 4349 7571 5627 5932 4126
#> 3 3446 4082 5117 5136 4178 6459 3529 4731 5448 2785 4526 7727 6294 5886 4537
#> 4 3720 4539 4542 6030 5672 6678 2998 2516 4864 2843 7117 8385 7332 6439 6331
#> 5 3165 4079 6173 5831 5414 6882 3556 2732 5588 2688 5995 8106 7174 6359 6328
#> 6 4116 4878 5205 5980 5591 6916 3949 2983 5603 3024 6150 7637 7262 6380 5961
#> 7 3354 4767 5389 4342 4277 6745 3537 4206 4318 2889 5052 7444 5544 5820 4346
#> 8 4529 3393 5248 4442 4476 4986 3088 4484 4001 3353 3713 5816 4117 5522 4321
#> 9 3136 4500 3780 5781 5407 5610 3061 3309 5553 2774 6430 8091 6920 6282 4889
#>     16   17   18   19   20
#> 1 5010 4415 3344 1632 4587
#> 2 5196 4211 4415 2282 4396
#> 3 5455 4749 4295 3059 5018
#> 4 6351 5161 5618 2233 4988
#> 5 6070 5454 4498 3073 5776
#> 6 5730 5807 5333 3011 5088
#> 7 5013 3862 5276 3211 4056
#> 8 4551 5243 2940 2634 4806
#> 9 5278 4989 5244 2735 4822
#> 
#> $grand_mean
#> [1] 4858
```
