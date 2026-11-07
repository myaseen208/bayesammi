# Priors of the Bayesian stability models

`bam_stab_priors()` sets the prior distributions of
[`bayes_stability()`](https://myaseen208.com/bayesammi/reference/bayes_stability.md):
normal priors for the overall mean and the contrasts of the genotype
effects, environment effects (the environmental index) and slope
deviations, and scaled inverse chi-squared priors for the error variance
and the deviation variances.

## Usage

``` r
bam_stab_priors(
  mu_var = NULL,
  g_var = NULL,
  h_var = NULL,
  b_var = 100,
  nu = 1,
  s2 = NULL,
  nu_d = 1,
  s2_d = NULL
)
```

## Arguments

- mu_var, g_var, h_var:

  prior variances of the overall mean and of the contrasts of the
  genotype and environment effects (`NULL`: \\10^4\\ times the variance
  of the response).

- b_var:

  prior variance of the contrasts of the slope deviations \\b_i\\
  (default 100; the slopes are \\1 + b_i\\, without unit).

- nu, s2:

  degrees of freedom and scale of the prior of the error variance
  (`s2 = NULL`: the within-cell mean square).

- nu_d, s2_d:

  degrees of freedom and scale of the prior of the deviation variances
  (`s2_d = NULL`: the pooled least-squares deviation variance).

## Value

An object of class `c("bam_stab_priors", "list")`.

## Details

\$\$\mu \sim N(0, \sigma^2\_\mu), \quad \mathbf{K}\_G\mathbf{g} \sim
N(\mathbf{0}, \sigma^2_g\mathbf{I}), \quad \mathbf{K}\_E\mathbf{h} \sim
N(\mathbf{0}, \sigma^2_h\mathbf{I}), \quad \mathbf{K}\_G\mathbf{b} \sim
N(\mathbf{0}, \sigma^2_b\mathbf{I}),\$\$ \$\$\sigma^2 \sim
\mathrm{Inv}\textrm{-}\chi^2(\nu, S^2_0), \qquad \sigma^2\_{d,i} \sim
\mathrm{Inv}\textrm{-}\chi^2(\nu_d, S^2_d),\$\$ where \\\mathbf{K}\_n\\
are the orthonormal contrasts of
[`matrix_k()`](https://myaseen208.com/bayesammi/reference/matrix_k.md),
so that the effects and the slope deviations sum to zero. Variances left
`NULL` are set from the data: \\10^4\\ times the variance of the
response for the effects, and the least-squares estimates for \\S^2_0\\
and \\S^2_d\\, so that the priors carry the information of \\\nu\\ and
\\\nu_d\\ observations whatever the unit of the trait.

## References

Lian, L. and de los Campos, G. (2016). FW: an R package for
Finlay-Wilkinson regression that incorporates genomic/pedigree
information and covariance structures between environments. *G3: Genes,
Genomes, Genetics*, **6**(3), 589–597.
[doi:10.1534/g3.115.026328](https://doi.org/10.1534/g3.115.026328)

## See also

[`bayes_stability`](https://myaseen208.com/bayesammi/reference/bayes_stability.md)

Other stability:
[`bam_ammi_indices()`](https://myaseen208.com/bayesammi/reference/bam_ammi_indices.md),
[`bam_indices()`](https://myaseen208.com/bayesammi/reference/bam_indices.md),
[`bayes_stability()`](https://myaseen208.com/bayesammi/reference/bayes_stability.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
bam_stab_priors(
    mu_var = NULL
  , g_var  = NULL
  , h_var  = NULL
  , b_var  = 100
  , nu     = 1
  , s2     = NULL
  , nu_d   = 1
  , s2_d   = NULL
  )
#> Priors of the Bayesian stability models
#>   mu       ~ N(0, 1e4 x var(y))
#>   K_G g    ~ N(0, 1e4 x var(y) I)
#>   K_E h    ~ N(0, 1e4 x var(y) I)
#>   K_G b    ~ N(0, 100 I)   (slopes 1 + b_i)
#>   sigma^2  ~ Inv-chi^2(1, within-cell mean square)
#>   sigma_d^2 ~ Inv-chi^2(1, least-squares deviation variance)
```
