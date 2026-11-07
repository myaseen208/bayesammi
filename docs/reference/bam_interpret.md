# Interpretation of a Bayesian AMMI analysis

`bam_interpret()` turns a result of
[`bayes_ammi()`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md),
[`biplots()`](https://myaseen208.com/bayesammi/reference/biplots.md) or
[`ge_ammi()`](https://myaseen208.com/bayesammi/reference/ge_ammi.md)
into a plain-language interpretation: whether the chains have converged,
which genotypes and environments differ from the average (main effects
whose HPD interval excludes zero), how strong the bilinear terms are,
which genotypes and environments interact significantly and in which
direction, which are stable, and what to check before acting on the
results. Every statement is computed from the numbers of the result.

## Usage

``` r
bam_interpret(x, higher_better = TRUE, top = 5L, ...)

# S3 method for class 'bam_interpretation'
print(x, ...)

# S3 method for class 'bam_interpretation'
format(x, style = c("text", "markdown"), width = getOption("width", 80L), ...)
```

## Arguments

- x:

  a result of
  [`bayes_ammi()`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md),
  [`biplots()`](https://myaseen208.com/bayesammi/reference/biplots.md)
  or
  [`ge_ammi()`](https://myaseen208.com/bayesammi/reference/ge_ammi.md).

- higher_better:

  `TRUE` (default) when larger values of the trait are better.

- top:

  number of genotypes or environments named in lists (default 5).

- ...:

  unused.

- style:

  `"text"` (default) or `"markdown"`.

- width:

  line width of the text (default `getOption("width")`).

## Value

An object of class `"bam_interpretation"`: a list with `analysis`,
`title`, `headline`, `findings`, `recommendations` and `cautions`
(sentences), `facts` (the numbers and lists behind the sentences),
`guide` (how the method is read) and `sources` (the full references of
the method;
[`bam_cite()`](https://myaseen208.com/bayesammi/reference/bam_cite.md)
gives them as BibTeX). [`print()`](https://rdrr.io/r/base/print.html)
shows the text; `format(x, style = "markdown")` returns it as Markdown.

## Details

The rules:

- convergence: potential scale reduction factor below 1.1 (two or more
  chains), effective sample size of at least 400, Raftery-Lewis
  dependence factor below 5 and the Heidelberger-Welch tests passed
  (Crossa et al. 2011);

- main effects: a genotype or environment differs from the average when
  the HPD interval of its effect excludes zero;

- bilinear terms: a term whose HPD interval of \\\lambda_k\\ reaches
  down to near zero (lower limit below a tenth of the posterior mean) is
  weakly supported;

- interaction: a genotype (environment) interacts significantly on the
  first term when the HPD interval of its first score excludes zero;
  genotypes and environments with significant scores of the same sign
  interact positively (the genotype does better there than its main
  effect predicts), of opposite sign negatively;

- stability
  ([`biplots()`](https://myaseen208.com/bayesammi/reference/biplots.md)):
  a genotype is stable when its HPD region contains the origin.

## References

Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K.,
Liu, G. and Cornelius, P. L. (2011). Bayesian estimation of the additive
main effects and multiplicative interaction model. *Crop Science*,
**51**(4), 1458–1469.
[doi:10.2135/cropsci2010.06.0343](https://doi.org/10.2135/cropsci2010.06.0343)

## See also

[`bayes_ammi`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md),
[`biplots`](https://myaseen208.com/bayesammi/reference/biplots.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(Maiz)
fm1 <-
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
    , .model  = c("AMMI", "GGE")
    , .het    = FALSE
    , .waic   = TRUE
    )
txt <-
  bam_interpret(
      x             = fm1
    , higher_better = TRUE
    , top           = 5
    )
txt
#> Bayesian AMMI analysis
#> 
#> 2 bilinear term(s): lambda_1 = 5,451, lambda_2 = 1,627; convergence is not
#> yet established.
#> 
#> Findings:
#> - Convergence is not established: effective sample size as small as 128.
#> - 12 of 132 Heidelberger-Welch stationarity and half-width tests failed
#>   (about 5% fail by chance at level 0.05).
#> - Overall mean 4,857 (95% HPD 4,791 to 4,914); error variance 761,800, a
#>   coefficient of variation of 18.0%.
#> - Genotypes significantly above average (HPD of the effect excludes zero): 4,
#>   5 and 6; below average: 1, 2, 7 and 8.
#> - Environments with significantly high means: 4, 6, 11, 12, 13 and 2 more;
#>   low: 1, 2, 7, 8, 10 and 2 more.
#> - Singular values with 95% HPD intervals: lambda_1 = 5,451 (4,531 to 6,278);
#>   lambda_2 = 1,627 (55.09 to 2,850).
#> - On the first bilinear term, genotypes with significantly positive scores:
#>   4, 5, 6 and 9; negative: 1, 2, 3 and 8.
#> - Environments with significantly positive first scores: 4, 11, 12, 13, 15
#>   and 1 more; negative: 1, 3, 8 and 10.
#> - Genotypes and environments with significant first scores of the same sign
#>   interact positively: the genotype yields more there than its main effect
#>   predicts; of opposite sign, less.
#> 
#> Recommendations:
#> - Recommend the genotypes with significant first scores for the environments
#>   whose scores have the same sign (biplots() shows the regions).
#> - Draw biplots() to see which genotypes and environments form homogeneous
#>   groups (overlapping HPD regions).
#> 
#> Cautions:
#> - Run longer chains (and two or more chains) before relying on the estimates.
#> - The HPD interval of lambda_2 reaches down to near zero: the term is weakly
#>   supported and its scores are uncertain.
#> 
#> Sources:
#> - Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K., Liu,
#>   G. and Cornelius, P. L. (2011). Bayesian estimation of the additive main
#>   effects and multiplicative interaction model. Crop Science, 51(4),
#>   1458-1469. doi:10.2135/cropsci2010.06.0343
#> - Spiegelhalter, D. J., Best, N. G., Carlin, B. P. and van der Linde, A.
#>   (2002). Bayesian measures of model complexity and fit. Journal of the Royal
#>   Statistical Society B, 64(4), 583-639. doi:10.1111/1467-9868.00353
#> - Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and
#>   widely applicable information criterion in singular learning theory.
#>   Journal of Machine Learning Research, 11, 3571-3594.
cat(format(x = txt, style = "markdown"))
#> ## Bayesian AMMI analysis  2 bilinear term(s): lambda_1 = 5,451, lambda_2 = 1,627; convergence is not yet established.  ### Findings - Convergence is not established: effective sample size as small as 128. - 12 of 132 Heidelberger-Welch stationarity and half-width tests failed   (about 5% fail by chance at level 0.05). - Overall mean 4,857 (95% HPD 4,791 to 4,914); error variance 761,800, a   coefficient of variation of 18.0%. - Genotypes significantly above average (HPD of the effect excludes zero): 4,   5 and 6; below average: 1, 2, 7 and 8. - Environments with significantly high means: 4, 6, 11, 12, 13 and 2 more;   low: 1, 2, 7, 8, 10 and 2 more. - Singular values with 95% HPD intervals: lambda_1 = 5,451 (4,531 to 6,278);   lambda_2 = 1,627 (55.09 to 2,850). - On the first bilinear term, genotypes with significantly positive scores:   4, 5, 6 and 9; negative: 1, 2, 3 and 8. - Environments with significantly positive first scores: 4, 11, 12, 13, 15   and 1 more; negative: 1, 3, 8 and 10. - Genotypes and environments with significant first scores of the same sign   interact positively: the genotype yields more there than its main effect   predicts; of opposite sign, less.  ### Recommendations - Recommend the genotypes with significant first scores for the environments   whose scores have the same sign (biplots() shows the regions). - Draw biplots() to see which genotypes and environments form homogeneous   groups (overlapping HPD regions).  ### Cautions - Run longer chains (and two or more chains) before relying on the estimates. - The HPD interval of lambda_2 reaches down to near zero: the term is weakly   supported and its scores are uncertain.  ### Sources - Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K., Liu,   G. and Cornelius, P. L. (2011). Bayesian estimation of the additive main   effects and multiplicative interaction model. Crop Science, 51(4),   1458-1469. doi:10.2135/cropsci2010.06.0343 - Spiegelhalter, D. J., Best, N. G., Carlin, B. P. and van der Linde, A.   (2002). Bayesian measures of model complexity and fit. Journal of the Royal   Statistical Society B, 64(4), 583-639. doi:10.1111/1467-9868.00353 - Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and   widely applicable information criterion in singular learning theory.   Journal of Machine Learning Research, 11, 3571-3594.
```
