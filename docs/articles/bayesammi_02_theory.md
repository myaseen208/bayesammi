# Model, priors and the Gibbs sampler

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

This vignette states the AMMI model, the prior distributions and the
conditional posterior distributions that
[`bayes_ammi()`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md)
samples. It follows Crossa et al. (2011), with each conditional written
out for the balanced design. Each formula is followed by a check against
the package.

## 2 Model

For genotype \\i = 1, \ldots, G\\ in replicate \\l = 1, \ldots, B\\ of
environment \\j = 1, \ldots, E\\, \\y\_{ijl} = \mu + \tau_i + \delta_j +
\sum\_{k=1}^{K}\lambda_k\alpha\_{ik}\gamma\_{jk} + e\_{ijl}, \qquad
e\_{ijl} \sim N(0, \sigma^2),\\ with \\\sum_i\tau_i = \sum_j\delta_j =
0\\ and \\\lambda_1 \ge \cdots \ge \lambda_K \> 0\\. The score vectors
\\\va_k = (\alpha\_{ik})\\ and \\\vg_k = (\gamma\_{jk})\\ have unit
length and are orthogonal to \\\vone\\ and to the other vectors of their
kind. With \\n = BGE\\ observations, the cell means \\\mY =
(\bar{y}\_{ij})\\ and the within-cell sum of squares \\\mathrm{SSW} =
\sum\_{ijl}(y\_{ijl} - \bar{y}\_{ij})^2\\ are sufficient statistics when
the design is balanced: \\\sum\_{ijl}(y\_{ijl} - \eta\_{ij})^2 =
\mathrm{SSW} + B\sum\_{ij}(\bar{y}\_{ij} - \eta\_{ij})^2.\\

## 3 Priors

[`bam_priors()`](https://myaseen208.com/bayesammi/reference/bam_priors.md)
sets the vague but proper priors of the paper: \\\mu \sim N(\mu_0,
\sigma^2\_\mu), \quad \vt \sim N(\mathbf{m}\_\tau, \sigma^2\_\tau\mI),
\quad \vd \sim N(\mathbf{m}\_\delta, \sigma^2\_\delta\mI), \quad
\lambda_k \sim N^+(\lambda_0, \sigma^2\_\lambda),\\ \\\sigma^2 \sim
\mathrm{Inv}\textrm{-}\chi^2(\nu, S^2_0),\\ with uniform priors for
\\\va_k\\ and \\\vg_k\\ on their constrained spheres. The variances are
\\10^8\\ and \\\nu = 1\\ by default. \\S^2_0\\ defaults to the
least-squares residual mean square, so the prior carries the information
of one observation whatever the unit of the trait. Hobert and Casella
(1996) show why \\\nu\\ and \\S^2_0\\ must be positive.

``` r

bam_priors(
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
```

    Priors of the Bayesian AMMI model
      mu            ~ N(0e+00, 1e+08)
      tau_i         ~ N(0e+00, 1e+08), sum zero
      delta_j       ~ N(0e+00, 1e+08), sum zero
      lambda_k      ~ N+(0e+00, 1e+08), ordered
      sigma^2       ~ Inv-chi^2(nu = 1e+00, s2 = least-squares residual mean square)
      alpha_k, gamma_k uniform on the constrained unit spheres

## 4 Conditional distributions

Write \\c_x = \sigma^2/\sigma^2_x\\ for the ratio of the error variance
to a prior variance.

### 4.1 Overall mean and main effects

The main effects are written in orthonormal contrasts, \\\vt =
\mK_G^\top\mathbf{t}\\ and \\\vd = \mK_E^\top\mathbf{d}\\, where
\\\mK_n\\ is the \\(n-1) \times n\\ matrix of
[`matrix_k()`](https://myaseen208.com/bayesammi/reference/matrix_k.md),
\\\mK_n\mK_n^\top = \mI\\ and \\\mK_n\vone = \mathbf{0}\\. The
interaction terms have zero row and column sums. In a balanced design
\\\mu\\, \\\mathbf{t}\\ and \\\mathbf{d}\\ are therefore conditionally
independent: \\\mu \sim N\left(\frac{n\bar{y} + c\_\mu\mu_0}{n +
c\_\mu}, \frac{\sigma^2}{n + c\_\mu}\right), \qquad \mathbf{t} \sim
N\left(\frac{BE\\\mK_G\bar{\vy}\_{G} +
c\_\tau\mK_G\mathbf{m}\_\tau}{BE + c\_\tau}, \frac{\sigma^2}{BE +
c\_\tau}\mI\right),\\ and \\\mathbf{d}\\ in the same way, with \\BG\\
and the environment means \\\bar{\vy}\_E\\. This is equation 5 of the
paper with \\(\mathbf{W}^\top\mathbf{W} + \mathbf{D})\\ evaluated.

``` r

K <-
  matrix_k(
      n = 4
    )
round(x = K %*% t(K), digits = 12)
```

         [,1] [,2] [,3]
    [1,]    1    0    0
    [2,]    0    1    0
    [3,]    0    0    1

``` r

round(x = K %*% rep(1, 4), digits = 12)
```

         [,1]
    [1,]    0
    [2,]    0
    [3,]    0

### 4.2 Singular values

Since \\\va_k \perp \vone, \va\_{k^\*}\\ and \\\vg_k \perp \vone,
\vg\_{k^\*}\\, the data enter the conditional of \\\lambda_k\\ only
through \\\va_k^\top\mY\vg_k\\, and \\\mathbf{X}\_3^\top\mathbf{X}\_3 =
B\\ in the paper’s notation: \\\lambda_k \sim
N\left(\frac{B\\\va_k^\top\mY\vg_k + c\_\lambda\lambda_0}{B +
c\_\lambda}, \frac{\sigma^2}{B + c\_\lambda}\right) \quad\text{truncated
to } (\lambda\_{k+1}, \lambda\_{k-1}),\\ with \\\lambda_0 = \infty\\
above the first and \\\lambda\_{K+1} = 0\\ below the last (equation 7).
The truncation keeps the order, and a draw that would leave the interval
is drawn from its tail by the inverse distribution function.

### 4.3 Score vectors

As a function of \\\va_k\\ the likelihood is
\\\exp\\(B\lambda_k/\sigma^2)\\\va_k^\top\mY\vg_k\\\\. This is a von
Mises-Fisher density on the unit sphere of the subspace orthogonal to
\\\vone\\ and to the other genotype vectors (equations 8 and 9). With
\\\mH_k\\ an orthonormal basis of that subspace (\\G - K\\ columns),
\\\va_k = \mH_k\mathbf{a}\\ and \\\mathbf{a} \sim
\mathrm{vMF}\left(\frac{B\lambda_k}{\sigma^2}\mH_k^\top\mY\vg_k\right),\\
drawn by the algorithm of Wood (1994). \\\vg_k\\ is drawn in the same
way from \\\mY^\top\va_k\\.

``` r

set.seed(seed = 1)
X <- t(replicate(n = 2000, expr = bayesammi:::.bam_rvmf(v = c(10, 0, 0))))
c(mean_first = mean(X[, 1]),
  theory     = besselI(x = 10, nu = 1.5) / besselI(x = 10, nu = 0.5))
```

    mean_first     theory
     0.8983546  0.9000000 

### 4.4 Error variance

\\\sigma^2 \sim \frac{\nu S^2_0 + \mathrm{SSW} +
B\sum\_{ij}(\bar{y}\_{ij} - \hat{\eta}\_{ij})^2}{\chi^2\_{\nu + n}}
\qquad \text{(equation 10)}.\\

## 5 Signs and identifiability

Replacing \\(\va_k, \vg_k)\\ by \\(-\va_k, -\vg_k)\\ leaves the model
unchanged, so the posterior has \\2^K\\ mirror images. Within one
iteration the sign of \\\va_k\\ follows that of \\\vg_k\\. Over a long
run, though, a weak term can move from one mirror image to the other,
and its posterior means then shrink towards zero. With `.align = TRUE`
every draw takes the orientation that agrees with the least-squares
vectors. This chooses one of the mirror images and leaves the model
unchanged.

## 6 The sampler

Every iteration draws \\\mu\\, \\\vt\\, \\\vd\\; then, for \\k = 1,
\ldots, K\\, \\\lambda_k\\, \\\va_k\\ and \\\vg_k\\; and finally
\\\sigma^2\\. The first chain starts at the least-squares fit and the
other chains at dispersed values, so that the Gelman-Rubin diagnostic
can detect a chain that has not mixed.

The sampler is written in C++ (**Rcpp** and **RcppArmadillo**) and takes
its random numbers from R, so `.seed` or
[`set.seed()`](https://rdrr.io/r/base/Random.html) reproduces a run. A
plain R version of the same sampler is kept in the package as a
reference; on a balanced trial with one error variance both return the
same draws (vignette 13).

## 7 Extensions

The same scheme carries over to three extensions, each treated in its
own vignette.

- **Unbalanced trials** (vignette 7). Missing plots and empty cells are
  drawn from \\N(\eta\_{ij}, \sigma^2_j)\\ at every iteration (data
  augmentation; Tanner and Wong 1987), after which the conditionals
  above apply to the completed data.
- **Environment-specific variances** (vignette 7). With \\e\_{ijl} \sim
  N(0, \sigma^2_j)\\, every \\\sigma^2_j\\ has a scaled inverse
  chi-squared conditional, \\\mu\\ and \\\vd\\ are drawn jointly, and
  the conditional of \\\vg_k\\ becomes a Bingham-von Mises-Fisher
  distribution (Hoff 2009), drawn by slice sampling one coordinate at a
  time.
- **The GGE model** (vignette 8). The genotype effects are dropped and
  \\\vg_k\\ is no longer orthogonal to \\\vone\\ (Perez-Elizalde et
  al. 2012).

The pointwise log-likelihood of every plot is accumulated during the
run, which gives DIC and WAIC without storing the draws (vignette 9).

## 8 References

Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K.,
Liu, G. and Cornelius, P. L. (2011). Bayesian estimation of the additive
main effects and multiplicative interaction model. *Crop Science*,
51(4), 1458–1469. doi:10.2135/cropsci2010.06.0343

Viele, K. and Srinivasan, C. (2000). Parsimonious estimation of
multiplicative interaction in analysis of variance using
Kullback-Leibler information. *Journal of Statistical Planning and
Inference*, 84(1–2), 201–219.

Wood, A. T. A. (1994). Simulation of the von Mises Fisher distribution.
*Communications in Statistics - Simulation and Computation*, 23(1),
157–164.

Hobert, J. P. and Casella, G. (1996). The effect of improper priors on
Gibbs sampling in hierarchical linear mixed models. *Journal of the
American Statistical Association*, 91(436), 1461–1473.

Mardia, K. V. and Jupp, P. E. (2000). *Directional Statistics*. Wiley,
Chichester.

Tanner, M. A. and Wong, W. H. (1987). The calculation of posterior
distributions by data augmentation. *Journal of the American Statistical
Association*, 82(398), 528–540.

Hoff, P. D. (2009). Simulation of the matrix Bingham-von Mises-Fisher
distribution, with applications to multivariate and relational data.
*Journal of Computational and Graphical Statistics*, 18(2), 438–456.

Perez-Elizalde, S., Jarquin, D. and Crossa, J. (2012). A general
Bayesian estimation method of linear-bilinear models applied to plant
breeding trials with genotype x environment interaction. *Journal of
Agricultural, Biological, and Environmental Statistics*, 17, 15–37.
doi:10.1007/s13253-011-0063-9
