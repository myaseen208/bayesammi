# Package index

## Bayesian AMMI

Priors, the Gibbs sampler of Crossa et al. (2011) and the biplot with
HPD regions.

- [`bam_priors()`](https://myaseen208.com/bayesammi/reference/bam_priors.md)
  [`print(`*`<bam_priors>`*`)`](https://myaseen208.com/bayesammi/reference/bam_priors.md)
  : Prior distributions of the Bayesian AMMI model
- [`bayes_ammi()`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md)
  : Bayesian estimation of the AMMI and GGE models
- [`biplots()`](https://myaseen208.com/bayesammi/reference/biplots.md) :
  Biplots with highest posterior density regions

## Decisions for breeders

Predictions, winners, probabilities of superior performance and
stability, ranks and stability measures.

- [`bam_predict()`](https://myaseen208.com/bayesammi/reference/bam_predict.md)
  : Predictions of every genotype-environment cell
- [`bam_winners()`](https://myaseen208.com/bayesammi/reference/bam_winners.md)
  : Which genotype wins where: posterior probabilities
- [`bam_superiority()`](https://myaseen208.com/bayesammi/reference/bam_superiority.md)
  : Probabilities of superior performance and stability
- [`bam_ranks()`](https://myaseen208.com/bayesammi/reference/bam_ranks.md)
  : Posterior distribution of the ranks of the genotypes
- [`bam_stability()`](https://myaseen208.com/bayesammi/reference/bam_stability.md)
  : Bayesian stability of the genotypes

## Bayesian stability analysis

Finlay-Wilkinson, Eberhart-Russell and Shukla models; posterior
distributions of stability indices.

- [`bayes_stability()`](https://myaseen208.com/bayesammi/reference/bayes_stability.md)
  : Bayesian stability analysis: Finlay-Wilkinson, Eberhart-Russell and
  Shukla models
- [`bam_stab_priors()`](https://myaseen208.com/bayesammi/reference/bam_stab_priors.md)
  : Priors of the Bayesian stability models
- [`bam_indices()`](https://myaseen208.com/bayesammi/reference/bam_indices.md)
  : Posterior distributions of classical stability indices
- [`bam_ammi_indices()`](https://myaseen208.com/bayesammi/reference/bam_ammi_indices.md)
  : Posterior distributions of AMMI-based stability indices

## Model choice

DIC and WAIC, posterior predictive checks and draws for posterior,
bayesplot and coda.

- [`bam_compare()`](https://myaseen208.com/bayesammi/reference/bam_compare.md)
  : Compare Bayesian AMMI and GGE models
- [`bam_ppc()`](https://myaseen208.com/bayesammi/reference/bam_ppc.md) :
  Posterior predictive checks
- [`as_draws_df(`*`<bam_bayes_ammi>`*`)`](https://myaseen208.com/bayesammi/reference/bam-draws.md)
  [`as_draws(`*`<bam_bayes_ammi>`*`)`](https://myaseen208.com/bayesammi/reference/bam-draws.md)
  [`as.mcmc.list(`*`<bam_bayes_ammi>`*`)`](https://myaseen208.com/bayesammi/reference/bam-draws.md)
  : The draws of a fit for posterior, bayesplot and coda

## Application

- [`bam_app()`](https://myaseen208.com/bayesammi/reference/bam_app.md)
  [`bam_app_object()`](https://myaseen208.com/bayesammi/reference/bam_app.md)
  : Bayesian AMMI and GGE analysis in a web browser

## Interpretation and references

- [`bam_interpret()`](https://myaseen208.com/bayesammi/reference/bam_interpret.md)
  [`print(`*`<bam_interpretation>`*`)`](https://myaseen208.com/bayesammi/reference/bam_interpret.md)
  [`format(`*`<bam_interpretation>`*`)`](https://myaseen208.com/bayesammi/reference/bam_interpret.md)
  : Interpretation of a Bayesian AMMI analysis
- [`bam_cite()`](https://myaseen208.com/bayesammi/reference/bam_cite.md)
  : References of a result or of the package

## Tables and least-squares AMMI

Cell means and variances, main and interaction effects, least-squares
AMMI and the mixed model.

- [`ge_mean()`](https://myaseen208.com/bayesammi/reference/ge_mean.md) :
  Genotype-by-environment means
- [`ge_var()`](https://myaseen208.com/bayesammi/reference/ge_var.md) :
  Genotype-by-environment variances
- [`g_eff()`](https://myaseen208.com/bayesammi/reference/g_eff.md) :
  Genotype effects
- [`e_eff()`](https://myaseen208.com/bayesammi/reference/e_eff.md) :
  Environment effects
- [`ge_eff()`](https://myaseen208.com/bayesammi/reference/ge_eff.md) :
  Genotype-by-environment interaction effects
- [`ge_ammi()`](https://myaseen208.com/bayesammi/reference/ge_ammi.md) :
  Least-squares AMMI model
- [`ge_model()`](https://myaseen208.com/bayesammi/reference/ge_model.md)
  : Linear mixed model of a multi-environment trial

## Matrices

- [`matrix_k()`](https://myaseen208.com/bayesammi/reference/matrix_k.md)
  : Orthonormal contrast matrix
- [`orthnorm()`](https://myaseen208.com/bayesammi/reference/orthnorm.md)
  : Gram-Schmidt orthonormalisation

## Data

- [`Maiz`](https://myaseen208.com/bayesammi/reference/Maiz.md) : Maize
  multi-environment trial

## Results and utilities

- [`print(`*`<bam_result>`*`)`](https://myaseen208.com/bayesammi/reference/bam-methods.md)
  [`summary(`*`<bam_result>`*`)`](https://myaseen208.com/bayesammi/reference/bam-methods.md)
  [`print(`*`<summary.bam_result>`*`)`](https://myaseen208.com/bayesammi/reference/bam-methods.md)
  [`format(`*`<bam_result>`*`)`](https://myaseen208.com/bayesammi/reference/bam-methods.md)
  [`as.data.table(`*`<bam_result>`*`)`](https://myaseen208.com/bayesammi/reference/bam-methods.md)
  [`as.data.frame(`*`<bam_result>`*`)`](https://myaseen208.com/bayesammi/reference/bam-methods.md)
  [`plot(`*`<bam_result>`*`)`](https://myaseen208.com/bayesammi/reference/bam-methods.md)
  : Methods for bayesammi results
- [`bam_theme()`](https://myaseen208.com/bayesammi/reference/bam_theme.md)
  : ggplot2 theme of bayesammi
- [`bam_save()`](https://myaseen208.com/bayesammi/reference/bam_save.md)
  : Save a figure (vector PDF by default)
- [`bam_notation()`](https://myaseen208.com/bayesammi/reference/bam_notation.md)
  : Notation of bayesammi
- [`bam_reference()`](https://myaseen208.com/bayesammi/reference/bam_reference.md)
  : Published values used by the tests and vignettes
- [`bayesammi`](https://myaseen208.com/bayesammi/reference/bayesammi-package.md)
  [`bayesammi-package`](https://myaseen208.com/bayesammi/reference/bayesammi-package.md)
  : bayesammi: Bayesian Estimation of the Additive Main Effects and
  Multiplicative Interaction Model
