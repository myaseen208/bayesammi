# Bayesian stability of the genotypes

`bam_stability()` gives the posterior distributions of two stability
measures of every genotype: the AMMI stability value (ASV) of Purchase
et al. (2000) computed draw by draw, and the variance of the interaction
effects over the environments (Shukla's stability variance of the
model), with posterior means, HPD intervals, the probability of being
the most stable genotype and the posterior expected ranks, combined with
the mean into a Bayesian yield stability index.

## Usage

``` r
bam_stability(fit, higher_better = TRUE)
```

## Arguments

- fit:

  a result of
  [`bayes_ammi()`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md)
  or
  [`bayes_stability()`](https://myaseen208.com/bayesammi/reference/bayes_stability.md)
  (`bam_stability()` needs
  [`bayes_ammi()`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md)).

- higher_better:

  `TRUE` (default) when larger values of the trait are better.

## Value

An object of class `c("bam_stability", "bam_result", "list")` with
`stability` (`data.table` sorted by the index: `Gen`, `mean`, `ASV`,
`ASV_lower`, `ASV_upper`, `variance`, `prob_most_stable` (smallest ASV),
`rank_mean`, `rank_ASV`, `YSI`). Methods:
[`print()`](https://rdrr.io/r/base/print.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`bam_interpret()`](https://myaseen208.com/bayesammi/reference/bam_interpret.md).

## Details

With the biplot scores \\u\_{ik} = \lambda_k^{1/2}\alpha\_{ik}\\ and the
sums of squares of the terms proportional to \\\lambda_k^2\\,
\$\$\mathrm{ASV}\_i =
\sqrt{\left(\frac{\bar{\lambda}\_1^2}{\bar{\lambda}\_2^2}u\_{i1}\right)^2 +
u\_{i2}^2}, \qquad V_i = \frac{1}{E - 1}\sum_j ge\_{ij}^2,\$\$
(\\\mathrm{ASV}\_i = \|u\_{i1}\|\\ with one term), evaluated in every
draw. The weight of the first axis uses the posterior means
\\\bar{\lambda}\_k\\ of the singular values: the ratio of the draws
themselves has no finite mean when the posterior of \\\lambda_2\\
reaches zero. The Bayesian yield stability index is the sum of the
posterior expected ranks of the mean and of the ASV.

## References

Purchase, J. L., Hatting, H. and van Deventer, C. S. (2000). Genotype x
environment interaction of winter wheat in South Africa: II. Stability
analysis of yield performance. *South African Journal of Plant and
Soil*, **17**(3), 101–107.

Shukla, G. K. (1972). Some statistical aspects of partitioning
genotype-environmental components of variability. *Heredity*, **29**(2),
237–245. [doi:10.1038/hdy.1972.87](https://doi.org/10.1038/hdy.1972.87)

## See also

[`bam_superiority`](https://myaseen208.com/bayesammi/reference/bam_superiority.md)

Other decisions:
[`bam_predict()`](https://myaseen208.com/bayesammi/reference/bam_predict.md),
[`bam_ranks()`](https://myaseen208.com/bayesammi/reference/bam_ranks.md),
[`bam_superiority()`](https://myaseen208.com/bayesammi/reference/bam_superiority.md),
[`bam_winners()`](https://myaseen208.com/bayesammi/reference/bam_winners.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(Maiz)
fit <-
  bayes_ammi(
      .data   = Maiz
    , .y      = y
    , .gen    = entry
    , .env    = site
    , .rep    = rep
    , .nIter  = 1000
    , .k      = 2
    , .burnin = 200
    , .thin   = 1
    , .chains = 1
    , .priors = bam_priors(
          mu_mean     = 0
        , mu_var      = 1e+8
        , tau_mean    = 0
        , tau_var     = 1e+8
        , delta_mean  = 0
        , delta_var   = 1e+8
        , lambda_mean = 0
        , lambda_var  = 1e+8
        , nu          = 1
        , s2          = NULL
        )
    , .level  = 0.95
    , .seed   = 1
    , .align  = TRUE
    , .model  = "AMMI"
    , .het    = FALSE
    , .waic   = TRUE
    )
st <-
  bam_stability(
      fit           = fit
    , higher_better = TRUE
    )
st$stability
#>       Gen     mean      ASV  ASV_lower ASV_upper  variance prob_most_stable
#>    <char>    <num>    <num>      <num>     <num>     <num>            <num>
#> 1:      6 5331.927 213.0762 107.982997  332.3294 128362.66            0.024
#> 2:      5 5248.385 221.0675  96.155067  335.6657 149841.43            0.021
#> 3:      3 4822.457 129.1675  12.772790  239.9464  57420.18            0.319
#> 4:      7 4665.152 101.0916   5.065334  205.0764  64520.18            0.554
#> 5:      9 4924.780 254.3651 138.100674  378.1377 177020.35            0.008
#> 6:      4 5219.423 405.2636 291.030365  517.2402 410786.63            0.000
#> 7:      2 4606.430 193.0657  75.960309  319.8921 116309.35            0.070
#> 8:      1 4613.857 301.9509 198.330401  434.0974 242789.42            0.004
#> 9:      8 4285.844 379.3936 256.537235  499.2179 404930.79            0.000
#>    rank_mean rank_ASV    YSI
#>        <num>    <num>  <num>
#> 1:     1.464    4.338  5.802
#> 2:     2.153    4.578  6.731
#> 3:     4.998    2.283  7.281
#> 4:     6.553    1.730  8.283
#> 5:     4.247    5.394  9.641
#> 6:     2.417    8.392 10.809
#> 7:     7.154    3.779 10.933
#> 8:     7.037    6.554 13.591
#> 9:     8.977    7.952 16.929
```
