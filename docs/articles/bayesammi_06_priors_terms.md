# Priors and the number of bilinear terms

\\ \newcommand{\vy}{\mathbf{y}} \newcommand{\vx}{\mathbf{x}}
\newcommand{\vb}{\mathbf{b}} \newcommand{\ve}{\mathbf{e}}
\newcommand{\vone}{\mathbf{1}} \newcommand{\vbeta}{\boldsymbol{\beta}}
\newcommand{\vmu}{\boldsymbol{\mu}} \newcommand{\mX}{\mathbf{X}}
\newcommand{\mZ}{\mathbf{Z}} \newcommand{\mV}{\mathbf{V}}
\newcommand{\mI}{\mathbf{I}} \newcommand{\mJ}{\mathbf{J}}
\newcommand{\mSigma}{\boldsymbol{\Sigma}} \newcommand{\mY}{\mathbf{Y}}
\newcommand{\mK}{\mathbf{K}} \newcommand{\mH}{\mathbf{H}}
\newcommand{\va}{\boldsymbol{\alpha}}
\newcommand{\vg}{\boldsymbol{\gamma}}
\newcommand{\vt}{\boldsymbol{\tau}}
\newcommand{\vd}{\boldsymbol{\delta}}
\newcommand{\vtheta}{\boldsymbol{\theta}} \newcommand{\vv}{\mathbf{v}}
\newcommand{\E}{\mathrm{E}} \newcommand{\Var}{\mathrm{Var}}
\newcommand{\Cov}{\mathrm{Cov}} \newcommand{\tr}{\operatorname{tr}} \\

## 1 Abstract

The defaults of
[`bam_priors()`](https://myaseen208.com/bayesammi/reference/bam_priors.md)
are vague: the data dominate the posterior. This vignette shows how to
set informative priors from earlier trials and how the posterior of the
last singular value guides the choice of the number of bilinear terms
\\K\\.

## 2 Vague priors

``` r

f0 <-
  bayes_ammi(
      .data   = Maiz
    , .y      = y
    , .gen    = entry
    , .env    = site
    , .rep    = rep
    , .nIter  = 3000
    , .k      = 2
    , .burnin = 500
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
    , .seed   = 6
    , .align  = TRUE
    , .model  = "AMMI"
    , .het    = FALSE
    , .waic   = TRUE
    )
f0$summary[block %in% c("lambda", "sigma2"), list(parameter, mean, sd, lower, upper)]
```

       parameter       mean         sd        lower      upper
          <char>      <num>      <num>        <num>      <num>
    1:    sigma2 763107.479 43844.8853 686109.01145 857075.929
    2: lambda[1]   5462.493   463.4891   4552.72387   6357.336
    3: lambda[2]   1682.008   783.5644     93.71569   2947.843

## 3 Informative priors

Suppose earlier trials in the region give an error variance of about \\8
\times 10^5\\, with the information of 50 plots, and a first singular
value near 5000:

``` r

f1 <-
  bayes_ammi(
      .data   = Maiz
    , .y      = y
    , .gen    = entry
    , .env    = site
    , .rep    = rep
    , .nIter  = 3000
    , .k      = 2
    , .burnin = 500
    , .thin   = 1
    , .chains = 1
    , .priors = bam_priors(
          mu_mean     = 0
        , mu_var      = 1e+8
        , tau_mean    = 0
        , tau_var     = 1e+8
        , delta_mean  = 0
        , delta_var   = 1e+8
        , lambda_mean = 5000
        , lambda_var  = 250^2
        , nu          = 50
        , s2          = 8e5
        )
    , .level  = 0.95
    , .seed   = 6
    , .align  = TRUE
    , .model  = "AMMI"
    , .het    = FALSE
    , .waic   = TRUE
    )
pm <- rbind(vague       = f0$summary[block %in% c("lambda", "sigma2"), mean],
            informative = f1$summary[block %in% c("lambda", "sigma2"), mean])
colnames(pm) <- f0$summary[block %in% c("lambda", "sigma2"), parameter]
round(pm)
```

                sigma2 lambda[1] lambda[2]
    vague       763107      5462      1682
    informative 767603      5077      4449

The prior on \\\lambda_k\\ applies to every term. An informative prior
suited to the first term therefore also pulls the second towards 5000
(here from 1682 to 4449), a reason to keep `lambda_var` large unless the
prior information is about all terms.

## 4 How many terms?

The paper fitted \\K = 2\\ terms. A term that the data do not support
has a singular value whose posterior reaches down to zero:

``` r

fits <-
  lapply(
      X   = 1:3
    , FUN = function(K) {
        bayes_ammi(
            .data   = Maiz
          , .y      = y
          , .gen    = entry
          , .env    = site
          , .rep    = rep
          , .nIter  = 3000
          , .k      = K
          , .burnin = 500
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
          , .seed   = 6
          , .align  = TRUE
          , .model  = "AMMI"
          , .het    = FALSE
          , .waic   = TRUE
          )
      }
    )
data.table::rbindlist(l = lapply(X = fits, FUN = function(f) {
  f$summary[block == "lambda", list(K = f$settings$k, parameter, mean, sd, lower, upper)]
}))
```

           K parameter      mean       sd        lower    upper
       <int>    <char>     <num>    <num>        <num>    <num>
    1:     1 lambda[1] 5451.2684 454.4051 4560.3608937 6317.846
    2:     2 lambda[1] 5462.4934 463.4891 4552.7238741 6357.336
    3:     2 lambda[2] 1682.0078 783.5644   93.7156918 2947.843
    4:     3 lambda[1] 5477.0740 454.1756 4520.6161688 6300.692
    5:     3 lambda[2] 1942.2296 606.0342  764.5888482 3112.183
    6:     3 lambda[3]  882.9223 565.8721    0.2552065 1868.912

With \\K = 2\\ the HPD interval of \\\lambda_2\\ reaches down to 94, and
with \\K = 3\\ that of \\\lambda_3\\ reaches down to 0.26. In both cases
the last term is barely supported, whereas the interval of \\\lambda_1\\
stays far from zero. With a third term in the model, \\\lambda_2\\ moves
away from zero because the order \\\lambda_2 \ge \lambda_3\\ gives the
third term the weakest pattern. The least-squares singular values (5923,
3070, 2552) do not show any of this, because they are biased upwards.
Cornelius et al. (1996) found two terms significant for these data, and
Liu (2001, unpublished) suggested two or three.

## 5 The posterior of the error variance with more terms

``` r

vapply(
    X         = fits
  , FUN       = function(f) f$summary[parameter == "sigma2", mean]
  , FUN.VALUE = numeric(length = 1)
  )
```

    [1] 781037.5 763107.5 753356.4

Each further term takes \\G + E - 1 - 2k\\ degrees of freedom from the
residual, so the error variance decreases only slowly once the terms are
weak.

## 6 Information criteria

WAIC compares the three fits by their estimated out-of-sample fit
(vignette 9):

``` r

names(fits) <- paste0("K", 1:3)
cp <- do.call(what = bam_compare, args = fits)
cp
```

    Comparison of Bayesian models (smaller is better)

        model terms variances      DIC      pD     WAIC  p_waic      se
       <char> <int>    <char>    <num>   <num>    <num>   <num>   <num>
    1:     K3     3    common 11860.21 74.7002 11864.35 72.1106 44.5008
    2:     K2     2    common 11862.78 68.2585 11866.86 66.4888 44.2687
    3:     K1     1    common 11864.96 53.7269 11868.55 53.3171 43.6160
       lambda_last_lower lambda_last_mean  dWAIC dWAIC_se
                   <num>            <num>  <num>    <num>
    1:            0.2552         882.9223 0.0000   0.0000
    2:           93.7157        1682.0078 2.5098   1.6986
    3:         4560.3609        5451.2684 4.2048   5.1263

    Sources: Spiegelhalter et al. (2002); Watanabe (2010); Vehtari et al. (2017).
    Full references: bam_cite(x).

WAIC prefers K3. `dWAIC` is the difference from the best model and
`dWAIC_se` its standard error, computed from the plot-by-plot
differences between the two models; a difference of less than about two
of its standard errors does not separate the models.

## 7 References

Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K.,
Liu, G. and Cornelius, P. L. (2011). Bayesian estimation of the additive
main effects and multiplicative interaction model. *Crop Science*,
51(4), 1458–1469. doi:10.2135/cropsci2010.06.0343

Cornelius, P. L., Crossa, J. and Seyedsadr, M. S. (1996). Statistical
tests and estimators of multiplicative models for
genotype-by-environment interaction. In M. S. Kang and H. G. Gauch
(eds), *Genotype-by-Environment Interaction*. CRC Press, Boca Raton,
199–234.

Cornelius, P. L. and Crossa, J. (1999). Prediction assessment of
shrinkage estimators of multiplicative models for multi-environment
cultivar trials. *Crop Science*, 39(4), 998–1009.

Hobert, J. P. and Casella, G. (1996). The effect of improper priors on
Gibbs sampling in hierarchical linear mixed models. *Journal of the
American Statistical Association*, 91(436), 1461–1473.

Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation
and widely applicable information criterion in singular learning theory.
*Journal of Machine Learning Research*, 11, 3571–3594.

Gelman, A., Hwang, J. and Vehtari, A. (2014). Understanding predictive
information criteria for Bayesian models. *Statistics and Computing*,
24(6), 997–1016.
