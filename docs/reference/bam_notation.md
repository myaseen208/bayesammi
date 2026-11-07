# Notation of bayesammi

The symbols used in the help pages and vignettes of bayesammi, with
their meaning and the argument or result component that holds them. The
notation follows one rule throughout: scalars are upper-case letters in
normal face (\\P\\, \\R\\), indices lower-case italic (\\i, j\\),
vectors lower-case bold (\\\mathbf{y}\\, \\\boldsymbol{\beta}\\) and
matrices upper-case bold (\\\mathbf{X}\\, \\\boldsymbol{\Sigma}\\).

## Usage

``` r
bam_notation(format = c("table", "path", "tex"))
```

## Arguments

- format:

  `"table"` (default: a `data.table` with columns `symbol` (LaTeX),
  `text` (plain text), `kind` (scalar, index, vector, matrix, set,
  function, operator), `meaning` and `argument`), `"path"` (the path of
  `notation.csv`) or `"tex"` (the path of `Notations.tex`, the notation
  appendix used in the articles).

## Value

See `format`.

## See also

Other bayesammi utilities:
[`bam_reference()`](https://myaseen208.com/bayesammi/reference/bam_reference.md),
[`bam_save()`](https://myaseen208.com/bayesammi/reference/bam_save.md),
[`bam_theme()`](https://myaseen208.com/bayesammi/reference/bam_theme.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
bam_notation(
    format = c("table", "path", "tex")[1]
  )
#>                                                                     symbol
#>                                                                     <char>
#>  1:                                                                      G
#>  2:                                                                      E
#>  3:                                                                      B
#>  4:                                                                      K
#>  5:                                                                      n
#>  6:                                                             i, j, l, k
#>  7:                                          \\mathbf{Y} = (\\bar{y}_{ij})
#>  8:                                                                   \\mu
#>  9:                                        \\boldsymbol{\\tau} = (\\tau_i)
#> 10:                                    \\boldsymbol{\\delta} = (\\delta_j)
#> 11:                                                             \\lambda_k
#> 12:                       \\boldsymbol{\\alpha}_k, \\boldsymbol{\\gamma}_k
#> 13:                                                              \\sigma^2
#> 14:                                                            \\sigma^2_j
#> 15:                                                             \\eta_{ij}
#> 16:                                                               g_i, h_j
#> 17:                                                                1 + b_i
#> 18:                                                                 d_{ij}
#> 19:                                                        \\sigma^2_{d,i}
#> 20:                                                          W_i, P_i, I_i
#> 21:                                                              \\theta_k
#> 22:                                                                 n_{ij}
#> 23:                                                                      S
#> 24:                                          \\mathrm{DIC}, \\mathrm{WAIC}
#> 25:                                                          \\mathbf{K}_n
#> 26: \\sigma^2_\\mu, \\sigma^2_\\tau, \\sigma^2_\\delta, \\sigma^2_\\lambda
#> 27:                                                            \\nu, S^2_0
#> 28:                                                          \\mathrm{SSW}
#> 29:                            \\mathbf{A}^\\top, \\mathbf{1}, \\mathbf{I}
#>                                                                     symbol
#>                                                                     <char>
#>                            text     kind
#>                          <char>   <char>
#>  1:                           G   scalar
#>  2:                           E   scalar
#>  3:                           B   scalar
#>  4:                           K   scalar
#>  5:                           n   scalar
#>  6:                     i j l k    index
#>  7:                           Y   matrix
#>  8:                          mu   scalar
#>  9:                         tau   vector
#> 10:                       delta   vector
#> 11:                    lambda_k   scalar
#> 12:             alpha_k gamma_k   vector
#> 13:                     sigma^2   scalar
#> 14:                   sigma_j^2   scalar
#> 15:                      eta_ij   scalar
#> 16:                     g_i h_j   scalar
#> 17:                     1 + b_i   scalar
#> 18:                        d_ij   scalar
#> 19:                 sigma_d_i^2   scalar
#> 20:                 W_i P_i I_i   scalar
#> 21:                     theta_k   scalar
#> 22:                        n_ij   scalar
#> 23:                           S   scalar
#> 24:                    DIC WAIC   scalar
#> 25:                         K_n   matrix
#> 26: v_mu v_tau v_delta v_lambda   scalar
#> 27:                     nu s0^2   scalar
#> 28:                         SSW   scalar
#> 29:                     A^T 1 I operator
#>                            text     kind
#>                          <char>   <char>
#>                                                                                          meaning
#>                                                                                           <char>
#>  1:                                                                          number of genotypes
#>  2:                                                                       number of environments
#>  3: number of replicates in each environment (largest number of plots in a cell when unbalanced)
#>  4:                                                                     number of bilinear terms
#>  5:                                                                   number of observations BGE
#>  6:                                              genotype, environment, replicate, bilinear term
#>  7:                                                       table of genotype-by-environment means
#>  8:                                                                                 overall mean
#>  9:                                                             genotype main effects (sum zero)
#> 10:                                                          environment main effects (sum zero)
#> 11:                                                            singular value of bilinear term k
#> 12:                               genotype and environment score vectors of term k (unit length)
#> 13:                                                                               error variance
#> 14:                                                error variance of environment j (.het = TRUE)
#> 15:                                                expected value of genotype i in environment j
#> 16:                              genotype effect and environmental index of the stability models
#> 17:                                               slope of genotype i on the environmental index
#> 18:                                 deviation of genotype i from its regression in environment j
#> 19:                    deviation (Eberhart-Russell) or stability (Shukla) variance of genotype i
#> 20:                 Wricke's ecovalence, Lin and Binns' superiority, Annicchiarico's reliability
#> 21:                                         share of term k in the fitted interaction (Za, WAAS)
#> 22:                                      number of observed plots of genotype i in environment j
#> 23:                                                               number of kept posterior draws
#> 24:                                          deviance and widely applicable information criteria
#> 25:                                     orthonormal contrasts, (n - 1) x n, rows orthogonal to 1
#> 26:                                                                              prior variances
#> 27:                              degrees of freedom and scale of the prior of the error variance
#> 28:                                                                   within-cell sum of squares
#> 29:                                                   transpose, vector of ones, identity matrix
#>                                                                                          meaning
#>                                                                                           <char>
#>     argument
#>       <char>
#>  1:     .gen
#>  2:     .env
#>  3:     .rep
#>  4:       .k
#>  5:         
#>  6:         
#>  7:       .y
#>  8:         
#>  9:         
#> 10:         
#> 11:         
#> 12:         
#> 13:         
#> 14:     .het
#> 15:         
#> 16:         
#> 17:         
#> 18:         
#> 19:   .model
#> 20:         
#> 21:         
#> 22:         
#> 23:   .nIter
#> 24:    .waic
#> 25:         
#> 26:  .priors
#> 27:  .priors
#> 28:         
#> 29:         
#>     argument
#>       <char>
```
