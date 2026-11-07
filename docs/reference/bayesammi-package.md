# bayesammi: Bayesian Estimation of the Additive Main Effects and Multiplicative Interaction Model

bayesammi estimates the additive main effects and multiplicative
interaction (AMMI) model of multi-environment trials in plant breeding
by Bayesian methods, with the Gibbs sampler of Crossa et al. (2011).
Every parameter has a posterior distribution, so the AMMI biplot gets
statistical inference: highest posterior density (HPD) regions show
which genotypes and environments interact significantly and which form
homogeneous groups. Beyond the article, the package fits unbalanced
trials, environment-specific error variances and the GGE (sites
regression) model, turns the posterior into decisions (winners,
probabilities of superior performance and stability, ranks, stability
measures, predictions), compares and checks models, and runs every
analysis in a built-in application. It also fits the Finlay-Wilkinson,
Eberhart-Russell and Shukla stability models by Gibbs sampling and gives
posterior distributions of classical and AMMI-based stability indices.
The package is for breeders and agronomists who analyse
genotype-by-environment interaction and want uncertainty statements
about their conclusions.

## Details

**Functions**, in the order of an analysis:

- [`ge_mean`](https://myaseen208.com/bayesammi/reference/ge_mean.md),
  [`ge_var`](https://myaseen208.com/bayesammi/reference/ge_var.md),
  [`g_eff`](https://myaseen208.com/bayesammi/reference/g_eff.md),
  [`e_eff`](https://myaseen208.com/bayesammi/reference/e_eff.md),
  [`ge_eff`](https://myaseen208.com/bayesammi/reference/ge_eff.md): cell
  means and variances, main effects, interaction effects and their
  singular value decomposition.

- [`ge_ammi`](https://myaseen208.com/bayesammi/reference/ge_ammi.md),
  [`ge_model`](https://myaseen208.com/bayesammi/reference/ge_model.md):
  least-squares AMMI and the mixed model of the error variance.

- [`bam_priors`](https://myaseen208.com/bayesammi/reference/bam_priors.md),
  [`bayes_ammi`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md):
  priors and the Gibbs sampler of the AMMI or GGE model, with a common
  or environment-specific error variances and missing plots or cells
  (chains, burn-in, thinning, posterior summaries with HPD intervals,
  convergence diagnostics, DIC and WAIC).

- [`biplots`](https://myaseen208.com/bayesammi/reference/biplots.md):
  the biplot with HPD regions.

- [`bam_predict`](https://myaseen208.com/bayesammi/reference/bam_predict.md),
  [`bam_winners`](https://myaseen208.com/bayesammi/reference/bam_winners.md),
  [`bam_superiority`](https://myaseen208.com/bayesammi/reference/bam_superiority.md),
  [`bam_ranks`](https://myaseen208.com/bayesammi/reference/bam_ranks.md),
  [`bam_stability`](https://myaseen208.com/bayesammi/reference/bam_stability.md):
  decisions for breeders.

- [`bam_compare`](https://myaseen208.com/bayesammi/reference/bam_compare.md),
  [`bam_ppc`](https://myaseen208.com/bayesammi/reference/bam_ppc.md),
  [`bam-draws`](https://myaseen208.com/bayesammi/reference/bam-draws.md):
  model choice, posterior predictive checks, and the draws for
  posterior, bayesplot and coda.

- [`bam_stab_priors`](https://myaseen208.com/bayesammi/reference/bam_stab_priors.md),
  [`bayes_stability`](https://myaseen208.com/bayesammi/reference/bayes_stability.md):
  Bayesian Finlay-Wilkinson, Eberhart-Russell and Shukla stability
  models;
  [`bam_indices`](https://myaseen208.com/bayesammi/reference/bam_indices.md),
  [`bam_ammi_indices`](https://myaseen208.com/bayesammi/reference/bam_ammi_indices.md):
  posterior distributions of classical and AMMI-based stability indices.

- [`bam_cite`](https://myaseen208.com/bayesammi/reference/bam_cite.md):
  the references of the methods behind a result (also printed under
  every result).

- [`bam_app`](https://myaseen208.com/bayesammi/reference/bam_app.md):
  the built-in application.

- [`bam_interpret`](https://myaseen208.com/bayesammi/reference/bam_interpret.md):
  a plain-language interpretation of a result.

- [`matrix_k`](https://myaseen208.com/bayesammi/reference/matrix_k.md),
  [`orthnorm`](https://myaseen208.com/bayesammi/reference/orthnorm.md):
  orthonormal contrasts and Gram-Schmidt orthonormalisation.

Every analysis returns the named list documented on its help page, with
a class and methods (`print`, `summary`, `plot`, `as.data.table`); see
[`bam-methods`](https://myaseen208.com/bayesammi/reference/bam-methods.md).
Utilities:
[`bam_theme`](https://myaseen208.com/bayesammi/reference/bam_theme.md),
[`bam_save`](https://myaseen208.com/bayesammi/reference/bam_save.md),
[`bam_notation`](https://myaseen208.com/bayesammi/reference/bam_notation.md),
[`bam_reference`](https://myaseen208.com/bayesammi/reference/bam_reference.md).

**Data.** [`Maiz`](https://myaseen208.com/bayesammi/reference/Maiz.md):
grain yield of 9 maize genotypes in 20 sites with 4 replicates (the
example of Crossa et al. 2011).

**Vignettes.** Start with
[`vignette("bayesammi_01_getting_started", package = "bayesammi")`](https://myaseen208.com/bayesammi/articles/bayesammi_01_getting_started.md);
`browseVignettes("bayesammi")` lists all of them, each as HTML and as
PDF.

## Notation

Scalars are upper-case letters in normal face, indices lower-case
italic, vectors lower-case bold and matrices upper-case bold; Greek
letters follow the same rule. The full table is
[`bam_notation()`](https://myaseen208.com/bayesammi/reference/bam_notation.md).

|  |  |
|----|----|
| \\G\\, \\E\\, \\B\\, \\K\\ | numbers of genotypes, environments, replicates, bilinear terms |
| \\\mathbf{Y}\\ | \\G \times E\\ table of cell means |
| \\\mu\\, \\\boldsymbol{\tau}\\, \\\boldsymbol{\delta}\\ | overall mean, genotype and environment effects |
| \\\lambda_k\\, \\\boldsymbol{\alpha}\_k\\, \\\boldsymbol{\gamma}\_k\\ | singular value and genotype and environment score vectors of term \\k\\ |
| \\\sigma^2\\, \\\sigma^2_j\\ | common and environment-specific error variance |
| \\\eta\_{ij}\\ | expected value of genotype \\i\\ in environment \\j\\ |
| \\g_i\\, \\h_j\\, \\1 + b_i\\ | genotype effect, environmental index and slope of the stability models |
| \\d\_{ij}\\, \\\sigma^2\_{d,i}\\ | deviation from the regression and its (genotype-specific) variance |

## References

Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K.,
Liu, G. and Cornelius, P. L. (2011). Bayesian estimation of the additive
main effects and multiplicative interaction model. *Crop Science*,
**51**(4), 1458–1469.
[doi:10.2135/cropsci2010.06.0343](https://doi.org/10.2135/cropsci2010.06.0343)

Perez-Elizalde, S., Jarquin, D. and Crossa, J. (2012). A general
Bayesian estimation method of linear-bilinear models applied to plant
breeding trials with genotype x environment interaction. *Journal of
Agricultural, Biological, and Environmental Statistics*, **17**(1),
15–37.
[doi:10.1007/s13253-011-0063-9](https://doi.org/10.1007/s13253-011-0063-9)

Viele, K. and Srinivasan, C. (2000). Parsimonious estimation of
multiplicative interaction in analysis of variance using
Kullback-Leibler information. *Journal of Statistical Planning and
Inference*, **84**(1–2), 201–219.

## See also

[`bam-methods`](https://myaseen208.com/bayesammi/reference/bam-methods.md),
[`bam_notation`](https://myaseen208.com/bayesammi/reference/bam_notation.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>), ORCID 0000-0002-5923-1714.
