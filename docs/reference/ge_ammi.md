# Least-squares AMMI model

`ge_ammi()` computes the least-squares estimates of the AMMI model with
all interaction axes: grand mean, genotype and environment effects,
singular values and vectors of the interaction, and the error variance.
They are the starting values of
[`bayes_ammi()`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md)
and the classical estimates it is compared with.

## Usage

``` r
ge_ammi(.data, .y, .gen, .env, .rep)

# Default S3 method
ge_ammi(.data, .y, .gen, .env, .rep)
```

## Arguments

- .data:

  a data frame with one row per plot.

- .y, .gen, .env, .rep:

  the columns of the response, genotypes, environments and replicates:
  bare names, strings or vectors.

## Value

A list (class `c("bam_ge_ammi", "bam_result", "list")`) with the
components of the earlier versions:

- `g`, `e`, `Rep`, `k`:

  numbers of genotypes, environments, replicates and axes.

- `ge_means`:

  the result of
  [`ge_mean()`](https://myaseen208.com/bayesammi/reference/ge_mean.md).

- `mu`, `sigma2`, `tau`:

  grand mean, error variance and its reciprocal (precision).

- `tao`, `delta`:

  genotype and environment effects
  ([`g_eff()`](https://myaseen208.com/bayesammi/reference/g_eff.md),
  [`e_eff()`](https://myaseen208.com/bayesammi/reference/e_eff.md)).

- `lambdas`, `alphas`, `gammas`:

  singular values and genotype and environment singular vectors of the
  interaction.

## Details

With the cell means \\\bar{y}\_{ij}\\ of \\B\\ replicates, \\\hat{\mu} =
\bar{y}\\, \\\hat{\tau}\_i = \bar{y}\_{i\cdot} - \bar{y}\\,
\\\hat{\delta}\_j = \bar{y}\_{\cdot j} - \bar{y}\\ and the singular
value decomposition of the interaction effects
([`ge_eff()`](https://myaseen208.com/bayesammi/reference/ge_eff.md)),
\\K = \min(G, E) - 1\\ axes. The error variance is the residual variance
of the model with random replicates within environments
([`ge_model()`](https://myaseen208.com/bayesammi/reference/ge_model.md)),
in closed form for balanced data: \$\$\hat{\sigma}^2 =
\frac{\mathrm{SSW} - \mathrm{SS}\_{B(E)}}{E(B - 1)(G - 1)}\$\$ when the
replicate mean square is not smaller, otherwise (the replicate variance
component estimated as zero) \\\mathrm{SSW}/(GE(B - 1))\\; SSW is the
within-cell and \\\mathrm{SS}\_{B(E)}\\ the replicate-within-environment
sum of squares.

**Change from the earlier versions.** The error variance came from
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html); the closed
form gives the same value for balanced data without needing lme4.

## References

Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K.,
Liu, G. and Cornelius, P. L. (2011). Bayesian estimation of the additive
main effects and multiplicative interaction model. \emphCrop Science,
\strong51(4), 1458–1469. \doi10.2135/cropsci2010.06.0343

Gauch, H. G. (1992). *Statistical Analysis of Regional Yield Trials:
AMMI Analysis of Factorial Designs*. Elsevier, Amsterdam.

## See also

[`bayes_ammi`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md),
[`ge_eff`](https://myaseen208.com/bayesammi/reference/ge_eff.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(Maiz)
fit <-
  ge_ammi(
      .data = Maiz
    , .y    = y
    , .gen  = entry
    , .env  = site
    , .rep  = rep
    )
fit$lambdas
#> [1] 5922.8913 3070.4543 2552.4065 2320.2216 1758.6707 1313.4148  792.6324
#> [8]  757.0277
fit$sigma2
#> [1] 602844.9
```
