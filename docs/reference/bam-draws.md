# The draws of a fit for posterior, bayesplot and coda

The draws of
[`bayes_ammi()`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md)
are a [`coda::mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html)
(component `draws`). These methods hand them to the packages posterior
and coda; every plot of bayesplot (`mcmc_trace()`, `mcmc_areas()`,
`mcmc_rank_overlay()`, ...) then works on them, as does
[`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html)
(rank-normalised \\\hat{R}\\, bulk and tail effective sample sizes).

## Usage

``` r
# S3 method for class 'bam_bayes_ammi'
as_draws_df(x, ...)

# S3 method for class 'bam_bayes_ammi'
as_draws(x, ...)

# S3 method for class 'bam_bayes_ammi'
as.mcmc.list(x, ...)
```

## Arguments

- x:

  a result of
  [`bayes_ammi()`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md).

- ...:

  unused.

## Value

`as_draws_df()`: a
[`posterior::draws_df`](https://mc-stan.org/posterior/reference/draws_df.html);
`as_draws()`: a
[`posterior::draws_list`](https://mc-stan.org/posterior/reference/draws_list.html);
`as.mcmc.list()`: the
[`coda::mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html).

## References

Gabry, J., Simpson, D., Vehtari, A., Betancourt, M. and Gelman, A.
(2019). Visualization in Bayesian workflow. *Journal of the Royal
Statistical Society A*, **182**(2), 389–402.
[doi:10.1111/rssa.12378](https://doi.org/10.1111/rssa.12378)

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Burkner, P.-C.
(2021). Rank-normalization, folding, and localization: an improved R-hat
for assessing convergence of MCMC. *Bayesian Analysis*, **16**(2),
667–718. [doi:10.1214/20-BA1221](https://doi.org/10.1214/20-BA1221)

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
    , .nIter  = 500
    , .k      = 2
    , .burnin = 100
    , .thin   = 1
    , .chains = 2
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
if (requireNamespace("posterior", quietly = TRUE)) {
  dr <- posterior::as_draws_df(x = fit)
  posterior::summarise_draws(
      .x = posterior::subset_draws(
          x        = dr
        , variable = c("mu", "lambda[1]", "lambda[2]")
        )
    )
}
#> # A tibble: 3 × 10
#>   variable   mean median    sd   mad    q5   q95  rhat ess_bulk ess_tail
#>   <chr>     <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
#> 1 mu        4857.  4858.  32.0  32.4 4805. 4909.  1.00    945.      984.
#> 2 lambda[1] 5451.  5453. 447.  440.  4710. 6177.  1.00    887.      878.
#> 3 lambda[2] 1627.  1688. 776.  798.   241. 2844.  1.04     65.7     270.
if (requireNamespace("bayesplot", quietly = TRUE)) {
  bayesplot::mcmc_trace(x = fit$draws, pars = c("lambda[1]", "lambda[2]"))
}
```
