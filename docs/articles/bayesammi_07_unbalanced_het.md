# Unbalanced trials and environment-specific variances

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

Real trials lose plots, and some genotypes are not grown in every
environment. The error variance also differs from one site to the next.
This vignette shows how
[`bayes_ammi()`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md)
handles both: missing plots and empty cells by data augmentation, and a
separate error variance for every environment with `.het = TRUE`. It
ends with a comparison of the two variance models on `Maiz`.

## 2 Missing plots and empty cells

### 2.1 The idea

With \\B\\ the largest number of plots in a cell, a cell \\(i, j)\\ with
\\n\_{ij} \< B\\ plots has \\B - n\_{ij}\\ missing plots. A cell with
\\n\_{ij} = 0\\ is an empty cell. At every iteration the sampler draws
the missing plots from the model, \\y^{\mathrm{mis}}\_{ijl} \mid \vtheta
\sim N(\eta\_{ij}, \sigma^2_j), \qquad \eta\_{ij} = \mu + \tau_i +
\delta_j + \sum\_{k=1}^{K}\lambda_k\alpha\_{ik}\gamma\_{jk},\\ and then
draws the parameters from the conditional distributions of the balanced
design (vignette 2) given the completed data. This is the data
augmentation of Tanner and Wong (1987). The draws of the parameters come
from their posterior given the observed plots only, so nothing is
invented: a missing plot adds no information, and an empty cell is
predicted from the row and column structure of the model. The missing
plots are assumed to be missing at random (Little and Rubin 2019).

The sampler needs only the sums, sums of squares and counts of the
observed plots in every cell, so the completed data are never stored.

### 2.2 An example

From `Maiz` we remove 60 plots at random and every plot of three
genotype-environment combinations:

``` r

set.seed(seed = 21)
mz   <- Maiz[-sample(x = nrow(Maiz), size = 60), ]
gone <- (mz$entry == 2 & mz$site == 5) | (mz$entry == 7 & mz$site == 11) | (mz$entry == 9 & mz$site == 18)
mz   <- mz[!gone, ]
table(plots = table(factor(mz$entry):factor(mz$site)))
```

    plots
      0   1   2   3   4
      3   1   5  47 124 

``` r

fu <-
  bayes_ammi(
      .data   = mz
    , .y      = y
    , .gen    = entry
    , .env    = site
    , .rep    = rep
    , .nIter  = 3000
    , .k      = 2
    , .burnin = 500
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
    , .seed   = 3
    , .align  = TRUE
    , .model  = "AMMI"
    , .het    = FALSE
    , .waic   = TRUE
    )
fu
```

    Bayesian AMMI with 2 multiplicative term(s): 9 genotypes, 20 environments, 648 plots (3 empty cells)
    2 chain(s) x 3000 kept draws (iterations 3000 after a burn-in of 500, thinning 1); 1.9 s

    Posterior means, standard deviations and 95% HPD intervals
        parameter       mean       sd      lower     upper   ess   rhat
           <char>      <num>    <num>      <num>     <num> <num>  <num>
     1:        mu   4865.378    34.85  4.797e+03   4931.90  4639 1.0023
     2:    sigma2 771133.151 45862.81  6.884e+05 866667.90  3493 1.0004
     3:    tau[1]   -235.197    95.45 -4.318e+02    -56.12  5393 1.0003
     4:    tau[2]   -287.087    99.35 -4.902e+02   -100.13  4427 1.0000
     5:    tau[3]    -16.772    97.76 -2.130e+02    165.23  5358 1.0002
     6:    tau[4]    368.080    96.72  1.815e+02    557.44  5813 1.0006
     7:    tau[5]    395.278    95.27  2.256e+02    600.17  5604 1.0004
     8:    tau[6]    439.592    99.69  2.415e+02    628.05  4629 1.0000
     9:    tau[7]   -197.581   105.05 -4.019e+02     10.22  3681 1.0001
    10:    tau[8]   -528.930   100.93 -7.224e+02   -328.54  4393 1.0000
    11:    tau[9]     62.618    99.87 -1.305e+02    260.72  4541 1.0007
    12:  delta[1]  -1246.974   147.62 -1.550e+03   -973.11  5313 0.9998
    13:  delta[2]   -603.360   155.02 -9.078e+02   -307.74  4487 1.0008
    14:  delta[3]    287.881   144.60 -8.507e+00    548.59  5619 0.9999
    15:  delta[4]    436.272   151.48  1.404e+02    727.13  4672 1.0000
    16:  delta[5]    -21.712   154.67 -3.280e+02    270.00  4102 0.9999
    17:  delta[6]   1476.541   147.84  1.181e+03   1755.92  5320 1.0010
    18:  delta[7]  -1603.056   145.40 -1.891e+03  -1313.72  5601 1.0000
    19:  delta[8]   -820.386   145.03 -1.102e+03   -538.55  5335 1.0000
    20:  delta[9]     56.727   146.65 -2.301e+02    345.23  5396 1.0000
    21: delta[10]  -1931.169   143.26 -2.212e+03  -1658.16  5660 1.0009
    22: delta[11]    429.876   168.20  9.977e+01    752.36  3441 1.0018
    23: delta[12]   2682.222   155.78  2.381e+03   2990.30  4191 1.0004
    24: delta[13]   1547.696   157.58  1.240e+03   1852.64  3458 1.0004
    25: delta[14]   1188.212   146.55  9.058e+02   1478.69  4743 0.9999
    26: delta[15]    213.078   146.76 -7.575e+01    495.96  4750 0.9999
    27: delta[16]    542.632   145.82  2.585e+02    829.31  5435 1.0004
    28: delta[17]    -60.923   156.00 -3.637e+02    241.27  4233 1.0003
    29: delta[18]   -408.313   163.28 -7.449e+02   -105.42  3649 1.0002
    30: delta[19]  -2172.886   151.46 -2.468e+03  -1880.16  4538 1.0019
    31: delta[20]      7.641   151.29 -2.878e+02    301.10  4954 1.0000
    32: lambda[1]   5106.531   477.92  4.121e+03   6013.13  3157 1.0016
    33: lambda[2]   1339.303   806.12  5.931e-02   2705.39   639 1.0065
        parameter       mean       sd      lower     upper   ess   rhat
           <char>      <num>    <num>      <num>     <num> <num>  <num>

    Scores of the bilinear terms: x$summary (blocks alpha, gamma); biplot coordinates: x$scores

    Sources: Crossa et al. (2011); Spiegelhalter et al. (2002); Watanabe (2010).
    Full references: bam_cite(x).

The print shows the number of empty cells.
[`bam_predict()`](https://myaseen208.com/bayesammi/reference/bam_predict.md)
gives the posterior of every cell mean \\\eta\_{ij}\\, including the
empty cells:

``` r

pu <-
  bam_predict(
      fit    = fu
    , .level = 0.95
    )
pu$cells[n == 0]
```

          Gen    Env     n observed     mean       sd    lower    upper
       <fctr> <fctr> <int>    <num>    <num>    <num>    <num>    <num>
    1:      2      5     0       NA 4369.859 279.9944 3794.720 4891.308
    2:      7     11     0       NA 4903.347 322.4456 4269.447 5523.364
    3:      9     18     0       NA 4761.894 311.5227 4165.312 5390.920

``` r

Maiz_cells <- data.table::as.data.table(x = Maiz)[, list(full = mean(y)), by = list(
    Gen = factor(entry)
  , Env = factor(site)
  )]
merge(x = pu$cells[n == 0], y = Maiz_cells, by = c("Gen", "Env"))[, list(
    Gen
  , Env
  , mean
  , lower
  , upper
  , full
  )]
```

    Key: <Gen, Env>
          Gen    Env     mean    lower    upper     full
       <fctr> <fctr>    <num>    <num>    <num>    <num>
    1:      2      5 4369.859 3794.720 4891.308 5201.475
    2:      7     11 4903.347 4269.447 5523.364 5051.575
    3:      9     18 4761.894 4165.312 5390.920 5244.400

The column `full` holds the cell means of the complete data, which the
fit never saw. 2 of the 3 HPD intervals contain them. The miss, genotype
2 in environment 5, is a cell whose interaction is not carried by the
two bilinear terms: the model predicts an empty cell only from the
pattern of the rest of the trial. The posterior standard deviations of
the empty cells are on average 1.1 times those of the cells with all
four plots.

``` r

plot(x = pu, which = "cells")
```

![](data:image/webp;base64,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)

## 3 Environment-specific error variances

### 3.1 The model

With `.het = TRUE` the plots of environment \\j\\ have their own
variance, \\e\_{ijl} \sim N(0, \sigma^2_j)\\, and every \\\sigma^2_j\\
has the scaled inverse chi-squared prior of the common variance. Write
\\w_j = 1/\sigma^2_j\\. The conditional distributions change as follows.

- \\\sigma^2_j\\: scaled inverse chi-squared with \\\nu + BG\\ degrees
  of freedom and the residual sum of squares of environment \\j\\.
- \\\tau_i\\: the data precision becomes \\B\sum_j w_j\\, the same for
  every genotype.
- \\\mu\\ and \\\delta_j\\: no longer independent, so \\(\mu,
  \mathbf{d})\\ is drawn jointly from its multivariate normal
  conditional.
- \\\lambda_k\\: mean and variance weighted by \\w_j\\.
- \\\va_k\\: still von Mises-Fisher, because the quadratic term \\\sum_j
  w_j\gamma\_{jk}^2\\ does not depend on \\\va_k\\.
- \\\vg_k\\: the quadratic term \\\sum_j w_j\gamma\_{jk}^2\\ now varies
  on the sphere. The conditional is a Bingham-von Mises-Fisher
  distribution (Hoff 2009), drawn one coordinate at a time in the
  eigenbasis of the quadratic form by slice sampling.

Perez-Elizalde et al. (2012) give the general linear-bilinear model with
these conditionals.

### 3.2 The fit

``` r

fh <-
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
    , .seed   = 4
    , .align  = TRUE
    , .model  = "AMMI"
    , .het    = TRUE
    , .waic   = TRUE
    )
fh$summary[block == "sigma2", list(parameter, mean, lower, upper, rhat)]
```

         parameter       mean      lower     upper      rhat
            <char>      <num>      <num>     <num>     <num>
     1:  sigma2[1]  899051.19  511891.64 1349890.6 1.0006236
     2:  sigma2[2]  359098.16  192026.52  555745.0 0.9998866
     3:  sigma2[3] 1358450.15  751176.80 2064839.8 1.0065478
     4:  sigma2[4]  265204.01  143568.41  406218.9 0.9998684
     5:  sigma2[5]  486396.32  270327.39  743027.4 1.0003240
     6:  sigma2[6]  542883.85  302578.54  865847.9 1.0000071
     7:  sigma2[7]  485314.39  269353.19  740659.9 1.0036182
     8:  sigma2[8] 1083672.16  514349.66 1779086.2 1.0011059
     9:  sigma2[9] 1735925.45 1033048.68 2624506.7 0.9998791
    10: sigma2[10]   90216.12   47092.66  139392.2 0.9998561
    11: sigma2[11]  987914.43  558046.74 1522479.9 1.0007315
    12: sigma2[12]  639387.39  368112.12  997579.5 1.0002132
    13: sigma2[13]  912822.78  513857.36 1405516.9 0.9999362
    14: sigma2[14] 1048576.95  599100.73 1565103.0 1.0005215
    15: sigma2[15]  582948.13  304895.26  882299.0 1.0003107
    16: sigma2[16]  301506.87  170387.93  455563.0 1.0000115
    17: sigma2[17]  589627.10  309397.79  904723.1 1.0035253
    18: sigma2[18] 1369988.68  742355.05 2056858.1 1.0005439
    19: sigma2[19] 2031831.15 1162647.09 3044867.7 0.9998868
    20: sigma2[20]  392824.32  208244.34  594441.1 1.0014998
         parameter       mean      lower     upper      rhat
            <char>      <num>      <num>     <num>     <num>

The posterior means range from 90,216 to 2,031,831, a ratio of 23. One
common variance cannot describe both extremes.

``` r

ggplot2::ggplot(
    data    = fh$summary[block == "sigma2"]
  , mapping = ggplot2::aes(x = reorder(label, mean), y = mean, ymin = lower, ymax = upper)
  ) +
  ggplot2::geom_pointrange(colour = "#2c7fb8") +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 6), labels = scales::label_comma()) +
  ggplot2::labs(x = "Environment", y = "Error variance (posterior mean, 95% HPD)") +
  ggplot2::theme_bw()
```

![](data:image/webp;base64,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)

## 4 Which variance model?

``` r

fc <-
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
    , .seed   = 4
    , .align  = TRUE
    , .model  = "AMMI"
    , .het    = FALSE
    , .waic   = TRUE
    )
cp <-
  bam_compare(
      common          = fc
    , environment_var = fh
    )
cp
```

    Comparison of Bayesian models (smaller is better)

                 model terms       variances      DIC      pD     WAIC  p_waic
                <char> <int>          <char>    <num>   <num>    <num>   <num>
    1: environment_var     2 per environment 11735.63 85.6869 11741.45 82.0777
    2:          common     2          common 11862.02 67.9290 11865.68 65.9309
            se lambda_last_lower lambda_last_mean    dWAIC dWAIC_se
         <num>             <num>            <num>    <num>    <num>
    1: 41.9641          910.1783         2093.077   0.0000   0.0000
    2: 44.2329           63.7490         1686.287 124.2365  22.4796

    Sources: Spiegelhalter et al. (2002); Watanabe (2010); Vehtari et al. (2017).
    Full references: bam_cite(x).

The environment-specific variances lower WAIC by 124, 5.5 times the
standard error of the difference. They also sharpen the interaction: the
lower HPD limit of \\\lambda_2\\ rises from 64 to 910, because the plots
of the noisy environments no longer blur the second term. The posterior
predictive check of the common-variance model tells the same story
through the ratio of the largest to the smallest environment variance:

``` r

pp <-
  bam_ppc(
      fit    = fc
    , .draws = 500
    , .seed  = 1
    )
pp
```

    Posterior predictive checks (500 replicated trials)

                                     statistic  observed  mean_rep p_value
                                        <char>     <num>     <num>   <num>
    1:                  variance of cell means 1.842e+06 1.901e+06   0.682
    2:                       largest cell mean 8.385e+03 8.608e+03   0.728
    3:                      smallest cell mean 1.632e+03 1.743e+03   0.672
    4:              variance of genotype means 1.247e+05 1.429e+05   0.662
    5:                    first singular value 5.923e+03 5.902e+03   0.482
    6:                    within-cell variance 7.559e+05 7.660e+05   0.526
    7: largest / smallest environment variance 2.582e+01 2.924e+00   0.000

    Sources: Gelman et al. (1996). Full references: bam_cite(x).

The observed ratio is 26, against 2.9 on average in the replicated
trials of the common-variance model (posterior predictive p-value
0.000). The other statistics are reproduced well. Vignette 9 returns to
these tools.

## 5 References

Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K.,
Liu, G. and Cornelius, P. L. (2011). Bayesian estimation of the additive
main effects and multiplicative interaction model. *Crop Science*,
51(4), 1458–1469. doi:10.2135/cropsci2010.06.0343

Perez-Elizalde, S., Jarquin, D. and Crossa, J. (2012). A general
Bayesian estimation method of linear-bilinear models applied to plant
breeding trials with genotype x environment interaction. *Journal of
Agricultural, Biological, and Environmental Statistics*, 17, 15–37.
doi:10.1007/s13253-011-0063-9

Hoff, P. D. (2009). Simulation of the matrix Bingham-von Mises-Fisher
distribution, with applications to multivariate and relational data.
*Journal of Computational and Graphical Statistics*, 18(2), 438–456.

Tanner, M. A. and Wong, W. H. (1987). The calculation of posterior
distributions by data augmentation. *Journal of the American Statistical
Association*, 82(398), 528–540.

Little, R. J. A. and Rubin, D. B. (2019). *Statistical Analysis with
Missing Data*, 3rd edn. Wiley, Hoboken.
