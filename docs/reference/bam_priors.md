# Prior distributions of the Bayesian AMMI model

`bam_priors()` sets the hyperparameters of the prior distributions of
[`bayes_ammi()`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md).
The defaults are the vague but proper priors of Crossa et al. (2011):
normal priors with mean zero and variance \\10^8\\ for the overall mean,
the main effects and the singular values, and a scaled inverse
chi-square prior with one degree of freedom for the error variance.

## Usage

``` r
bam_priors(
  mu_mean = 0,
  mu_var = 1e+08,
  tau_mean = 0,
  tau_var = 1e+08,
  delta_mean = 0,
  delta_var = 1e+08,
  lambda_mean = 0,
  lambda_var = 1e+08,
  nu = 1,
  s2 = NULL
)

# S3 method for class 'bam_priors'
print(x, ...)
```

## Arguments

- mu_mean, mu_var:

  prior mean and variance of the overall mean \\\mu\\.

- tau_mean, tau_var:

  prior mean (a scalar or a vector of length \\G\\) and variance of the
  genotype effects \\\tau_i\\.

- delta_mean, delta_var:

  prior mean (scalar or length \\E\\) and variance of the environment
  effects \\\delta_j\\.

- lambda_mean, lambda_var:

  prior mean and variance of the singular values \\\lambda_k\\.

- nu, s2:

  degrees of freedom \\\nu\\ and scale \\S^2_0\\ of the prior of the
  error variance (`s2 = NULL`: the least-squares residual mean square).

- x:

  a `"bam_priors"` object.

- ...:

  unused.

## Value

A list of class `"bam_priors"` with these elements.

## Details

With \\G\\ genotypes, \\E\\ environments and \\K\\ bilinear terms,
\$\$\mu \sim N(\mu_0, \sigma^2\_\mu), \quad \boldsymbol{\tau} \sim
N(\mathbf{m}\_\tau, \sigma^2\_\tau \mathbf{I}\_G), \quad
\boldsymbol{\delta} \sim N(\mathbf{m}\_\delta, \sigma^2\_\delta
\mathbf{I}\_E),\$\$ each main effect restricted to sum zero,
\$\$\lambda_k \sim N^{+}(\lambda_0, \sigma^2\_\lambda), \quad \lambda_1
\ge \lambda_2 \ge \cdots \ge \lambda_K \> 0, \qquad \sigma^2 \sim
\mathrm{Inv}\textrm{-}\chi^2(\nu, S^2_0),\$\$ and spherical uniform
priors for the singular vectors on the subspaces orthogonal to
\\\mathbf{1}\\ and to the other vectors. \\S^2_0\\ = `NULL` (default)
takes the residual mean square of the least-squares AMMI fit, so that
the prior is equally vague whatever the unit of the trait; with \\\nu =
1\\ it carries the information of a single observation. Crossa et al.
(2011) report \\S_0 = 10^8\\; their published posterior of \\\sigma^2\\
is reproduced with the vague default (see the vignette on validation).

## References

Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K.,
Liu, G. and Cornelius, P. L. (2011). Bayesian estimation of the additive
main effects and multiplicative interaction model. *Crop Science*,
**51**(4), 1458–1469.
[doi:10.2135/cropsci2010.06.0343](https://doi.org/10.2135/cropsci2010.06.0343)

Hobert, J. P. and Casella, G. (1996). The effect of improper priors on
Gibbs sampling in hierarchical linear mixed models. *Journal of the
American Statistical Association*, **91**(436), 1461–1473.

## See also

[`bayes_ammi`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
bam_priors(
    mu_mean     = 0
  , mu_var      = 1e8
  , tau_mean    = 0
  , tau_var     = 1e8
  , delta_mean  = 0
  , delta_var   = 1e8
  , lambda_mean = 0
  , lambda_var  = 1e8
  , nu          = 1
  , s2          = NULL
  )
#> Priors of the Bayesian AMMI model
#>   mu            ~ N(0e+00, 1e+08)
#>   tau_i         ~ N(0e+00, 1e+08), sum zero
#>   delta_j       ~ N(0e+00, 1e+08), sum zero
#>   lambda_k      ~ N+(0e+00, 1e+08), ordered
#>   sigma^2       ~ Inv-chi^2(nu = 1e+00, s2 = least-squares residual mean square)
#>   alpha_k, gamma_k uniform on the constrained unit spheres
```
