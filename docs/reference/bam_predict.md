# Predictions of every genotype-environment cell

`bam_predict()` gives the posterior distribution of the expected
response of every genotype in every environment, including the cells
that were not observed, and of the genotype means over the environments.

## Usage

``` r
bam_predict(fit, .level = NULL)
```

## Arguments

- fit:

  a result of
  [`bayes_ammi()`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md)
  or
  [`bayes_stability()`](https://myaseen208.com/bayesammi/reference/bayes_stability.md)
  ([`bam_stability()`](https://myaseen208.com/bayesammi/reference/bam_stability.md)
  needs
  [`bayes_ammi()`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md)).

- .level:

  probability of the HPD intervals (default: that of the fit).

## Value

An object of class `c("bam_predict", "bam_result", "list")`:

- `cells`:

  `data.table`: `Gen`, `Env`, `n` (plots observed), `observed` (their
  mean, `NA` when missing), `mean`, `sd`, `lower`, `upper`.

- `genotypes`:

  `data.table`: posterior mean, sd and HPD interval of every genotype
  mean.

Methods: [`print()`](https://rdrr.io/r/base/print.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html)
(`which = c("cells", "genotypes")`),
[`bam_interpret()`](https://myaseen208.com/bayesammi/reference/bam_interpret.md).

## Details

For draw \\s\\, \$\$\eta^{(s)}\_{ij} = \mu^{(s)} + \tau^{(s)}\_i +
\delta^{(s)}\_j + \sum\_{k=1}^{K}\lambda^{(s)}\_k
\alpha^{(s)}\_{ik}\gamma^{(s)}\_{jk}\$\$ (without \\\tau_i\\ in the GGE
model); the genotype mean over the environments is
\\\bar{\eta}\_{i\cdot} = E^{-1}\sum_j\eta\_{ij}\\. A missing cell is
predicted from the main effects and from the interaction pattern its
genotype and environment show in the other cells; the data augmentation
of
[`bayes_ammi()`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md)
has already propagated its uncertainty into all parameters.

## References

Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K.,
Liu, G. and Cornelius, P. L. (2011). Bayesian estimation of the additive
main effects and multiplicative interaction model. *Crop Science*,
**51**(4), 1458–1469.
[doi:10.2135/cropsci2010.06.0343](https://doi.org/10.2135/cropsci2010.06.0343)

## See also

[`bam_winners`](https://myaseen208.com/bayesammi/reference/bam_winners.md),
[`bam_superiority`](https://myaseen208.com/bayesammi/reference/bam_superiority.md)

Other decisions:
[`bam_ranks()`](https://myaseen208.com/bayesammi/reference/bam_ranks.md),
[`bam_stability()`](https://myaseen208.com/bayesammi/reference/bam_stability.md),
[`bam_superiority()`](https://myaseen208.com/bayesammi/reference/bam_superiority.md),
[`bam_winners()`](https://myaseen208.com/bayesammi/reference/bam_winners.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(Maiz)
d <- Maiz[!(Maiz$entry == 8 & Maiz$site == 8), ]
fit <-
  bayes_ammi(
      .data   = d
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
pr <-
  bam_predict(
      fit    = fit
    , .level = 0.95
    )
pr$cells[n == 0]
#>       Gen    Env     n observed     mean       sd    lower    upper
#>    <fctr> <fctr> <int>    <num>    <num>    <num>    <num>    <num>
#> 1:      8      8     0       NA 4986.848 2063.612 815.0165 7808.311
```
