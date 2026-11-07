# Changelog

## bayesammi 1.0.0

### Changes visible to users

- Every function keeps its name, its arguments and the names of its
  results, so existing scripts run unchanged. Results print, summarise
  and plot ([`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html)), and any
  component converts to a table with
  `data.table::as.data.table(x, component = "...")`.
- [`bayes_ammi()`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md):
  new arguments `.k` (number of bilinear terms, default 2 as in Crossa
  et al. 2011; earlier versions fitted all min(G, E) - 1), `.burnin`,
  `.thin`, `.chains`, `.priors`
  ([`bam_priors()`](https://myaseen208.com/bayesammi/reference/bam_priors.md)),
  `.level`, `.seed` (the random number state of the session is restored)
  and `.align` (sign of every term); `.nIter` now has a default (1000).
  New components `summary` (posterior mean, SD, Monte Carlo error, HPD
  interval, effective sample size, potential scale reduction factor of
  every parameter), `scores`, `draws`
  ([`coda::mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html)),
  `diagnostics` (Gelman-Rubin, Raftery-Lewis, Heidelberger-Welch), `ls`,
  `settings`.
- `tau1` now holds the draws of the error variance sigma^2 (it held a
  constant precision).
- [`biplots()`](https://myaseen208.com/bayesammi/reference/biplots.md):
  HPD regions of the biplot coordinates lambda_k^(1/2) alpha_ik of the
  paper, at one or more levels (`pb = c(0.05, 0.01)`); new arguments
  `.which`, `.grid`; `biplot_data` gains `axis1_lower`, `axis1_upper`,
  `axis1_significant` and `group`.
- New
  [`bam_priors()`](https://myaseen208.com/bayesammi/reference/bam_priors.md),
  [`bam_interpret()`](https://myaseen208.com/bayesammi/reference/bam_interpret.md);
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of a
  [`bayes_ammi()`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md)
  result draws traces, histograms, the biplot of posterior means and the
  main effects with HPD intervals.
- Genotypes and environments keep their labels in all tables.

### Beyond the article

- Unbalanced trials: missing plots and whole genotype-environment cells
  are handled by data augmentation inside the sampler (earlier versions
  required balanced data).
- [`bayes_ammi()`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md):
  new arguments `.model` (`"AMMI"` or `"GGE"`, the sites regression
  model of Perez-Elizalde et al. 2012), `.het` (an error variance for
  every environment; the environment scores are then drawn from their
  Bingham-von Mises-Fisher conditional) and `.waic`; new component `fit`
  (DIC, WAIC and the contribution of every plot to WAIC).
- Decisions:
  [`bam_predict()`](https://myaseen208.com/bayesammi/reference/bam_predict.md)
  (posterior of every cell and genotype mean, including empty cells),
  [`bam_winners()`](https://myaseen208.com/bayesammi/reference/bam_winners.md)
  (probability that each genotype wins in each environment, groups of
  environments),
  [`bam_superiority()`](https://myaseen208.com/bayesammi/reference/bam_superiority.md)
  (probabilities of superior performance and stability, Dias et
  al. 2022),
  [`bam_ranks()`](https://myaseen208.com/bayesammi/reference/bam_ranks.md)
  and
  [`bam_stability()`](https://myaseen208.com/bayesammi/reference/bam_stability.md)
  (Bayesian AMMI stability value, interaction variance, yield stability
  index).
- Model choice:
  [`bam_compare()`](https://myaseen208.com/bayesammi/reference/bam_compare.md)
  (DIC, WAIC and the standard error of WAIC differences) and
  [`bam_ppc()`](https://myaseen208.com/bayesammi/reference/bam_ppc.md)
  (posterior predictive checks).
- `as_draws_df()`, `as_draws()` and `as.mcmc.list()` methods, for the
  **posterior**, **bayesplot** and **coda** packages.
- [`bam_app()`](https://myaseen208.com/bayesammi/reference/bam_app.md)
  and
  [`bam_app_object()`](https://myaseen208.com/bayesammi/reference/bam_app.md):
  a built-in **shiny** application with tabs for the interpretation, the
  plot, the table and the printed output, each with copy buttons, and
  the R code of the fit.
- All plots are **ggplot2** with
  [`scales::pretty_breaks()`](https://scales.r-lib.org/reference/pretty_breaks.html).
  [`bam_save()`](https://myaseen208.com/bayesammi/reference/bam_save.md)
  writes lossless WebP; the HTML vignettes and the application use
  lossless WebP images (about half the size of PNG).
- Bayesian stability analysis:
  [`bayes_stability()`](https://myaseen208.com/bayesammi/reference/bayes_stability.md)
  and
  [`bam_stab_priors()`](https://myaseen208.com/bayesammi/reference/bam_stab_priors.md)
  fit the Finlay-Wilkinson (1963), Eberhart-Russell (1966) and
  Shukla (1972) models by a partially collapsed Gibbs sampler (compiled,
  with an R reference giving the same draws), with the environmental
  index estimated jointly, unbalanced data, slopes, deviation and
  stability variances, the safety-first index of Eskridge (1990), DIC
  and WAIC. The decision tools,
  [`bam_compare()`](https://myaseen208.com/bayesammi/reference/bam_compare.md)
  and
  [`bam_ppc()`](https://myaseen208.com/bayesammi/reference/bam_ppc.md)
  accept these fits.
- [`bam_indices()`](https://myaseen208.com/bayesammi/reference/bam_indices.md):
  posterior distributions of Wricke’s ecovalence, Shukla’s stability
  variance, Lin and Binns’ superiority, Annicchiarico’s reliability and
  Eskridge’s safety-first index, from any fit.
  [`bam_ammi_indices()`](https://myaseen208.com/bayesammi/reference/bam_ammi_indices.md):
  ASV, SIPC, EV, Za, WAAS and WAASY from AMMI or GGE fits. Each index
  has HPD intervals, probabilities of being best and expected ranks.
- Citations: every printed result and interpretation ends with the
  sources of its method;
  [`bam_cite()`](https://myaseen208.com/bayesammi/reference/bam_cite.md)
  gives them as text or BibTeX; `citation("bayesammi")` lists the
  package and every method; `inst/REFERENCES.bib` has all references.
  The application shows the sources under each interpretation.
- The application gains a Stability menu (regression models, stability
  indices, AMMI indices).
- Seven new vignettes: unbalanced trials and environment-specific
  variances, the GGE model, decision tools, the Bayesian stability
  models, posterior stability indices, the application, and the
  references of the methods.

### Corrections

- [`bayes_ammi()`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md):
  the singular values were computed once before the loop and never
  sampled (`lambdas1` was constant); they are now drawn every iteration
  from their truncated normal conditional.
- [`bayes_ammi()`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md):
  the residuals for sigma^2 compared every cell with the environment of
  a data row instead of its own environment (with `Maiz`, always site
  1); the within-cell variances were summed in place of the sums of
  squares; the covariance matrices of the main effects were square roots
  of variances; the environment effects used B E in place of B G.
- [`bayes_ammi()`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md):
  the scores were drawn on the unconstrained Stiefel manifold and then
  orthogonalised; they are now von Mises-Fisher draws on the constrained
  spheres (Crossa et al. 2011, eq. 8-9).
- [`bayes_ammi()`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md):
  the prior of sigma^2 was improper (0 degrees of freedom); it is now
  proper.
- [`biplots()`](https://myaseen208.com/bayesammi/reference/biplots.md):
  the contours were of the singular vectors with a fixed bandwidth of
  0.2, and “stable” was decided by the bounding box of a contour.
- `Maiz` help page: it described another data set; the data have 720
  rows (9 genotypes, 20 sites, 4 replicates).

### Internal

- The sampler is written in C++ (**Rcpp**, **RcppArmadillo**) and draws
  its random numbers from R; a plain R version is kept as a reference,
  and both give the same draws from the same seed.
- The tidyverse, `rlang`, `rstiefel`, `tmvtnorm`, `mvtnorm` and `ks` are
  no longer needed; `lme4` is suggested (only
  [`ge_model()`](https://myaseen208.com/bayesammi/reference/ge_model.md)
  uses it;
  [`ge_ammi()`](https://myaseen208.com/bayesammi/reference/ge_ammi.md)
  computes the error variance in closed form). The diagnostics use
  `coda`.
- Tests against Tables 1 and B1 of Crossa et al. (2011), simulated data,
  the building blocks of the sampler and the results of version 0.3.0
  for the tables.

## bayesammi 0.3.0

CRAN release: 2024-11-22

- Added the
  [`biplots()`](https://myaseen208.com/bayesammi/reference/biplots.md)
  function.

## bayesammi 0.2.0

CRAN release: 2024-10-15

- Updated the DESCRIPTION file to meet the requirements of CRAN.

## bayesammi 0.1.0

CRAN release: 2018-04-12

- First version of bayesammi.
