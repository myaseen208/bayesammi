# Linear mixed model of a multi-environment trial

`ge_model()` fits the two-way model with interaction and random
replicates within environments by restricted maximum likelihood
([`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html)). Its residual
variance is the error variance of the least-squares AMMI fit
([`ge_ammi()`](https://myaseen208.com/bayesammi/reference/ge_ammi.md),
which computes it in closed form).

## Usage

``` r
ge_model(.data, .y, .gen, .env, .rep)

# Default S3 method
ge_model(.data, .y, .gen, .env, .rep)
```

## Arguments

- .data:

  a data frame with one row per plot.

- .y, .gen, .env, .rep:

  the columns of the response, genotypes, environments and replicates:
  bare names, strings or vectors.

## Value

The `lmerMod` object of
`lme4::lmer(Y ~ Env + Gen + Env:Gen + (1 | Env:Rep))`.

## Details

\$\$y\_{ijl} = \mu + \delta_j + \tau_i + (\tau\delta)\_{ij} +
b\_{l(j)} + e\_{ijl}, \qquad b\_{l(j)} \sim N(0, \sigma^2_b), \quad
e\_{ijl} \sim N(0, \sigma^2).\$\$ The package lme4 is needed (it is
suggested, not required, since version 1.0.0).

## References

Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K.,
Liu, G. and Cornelius, P. L. (2011). Bayesian estimation of the additive
main effects and multiplicative interaction model. \emphCrop Science,
\strong51(4), 1458–1469. \doi10.2135/cropsci2010.06.0343

Bates, D., Maechler, M., Bolker, B. and Walker, S. (2015). Fitting
linear mixed-effects models using lme4. *Journal of Statistical
Software*, **67**(1), 1–48.
[doi:10.18637/jss.v067.i01](https://doi.org/10.18637/jss.v067.i01)

## See also

[`ge_ammi`](https://myaseen208.com/bayesammi/reference/ge_ammi.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(Maiz)
if (requireNamespace("lme4", quietly = TRUE)) {
  fm <-
    ge_model(
        .data = Maiz
      , .y    = y
      , .gen  = entry
      , .env  = site
      , .rep  = rep
      )
  stats::sigma(fm)^2
}
#> [1] 602844.8
```
