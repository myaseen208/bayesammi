# Progress of bayesammi 1.0.0

Phase: P5 done (cores, interface, help pages, tests); next P6 vignettes.

| function | state |
|---|---|
| bayes_ammi | rewritten Gibbs sampler (Crossa et al. 2011 eq. 5-10); chains, burn-in, thinning, priors, HPD, diagnostics |
| biplots | HPD regions of lambda^(1/2) alpha by kernel density; stable = region contains origin |
| bam_priors, bam_interpret | new |
| ge_mean, ge_var, g_eff, e_eff, ge_eff, ge_ammi, ge_model, matrix_k, orthnorm | same results (legacy snapshot), data.table/closed forms |

Decisions: .k = 2 by default (paper); signs aligned to least squares (.align); s2 prior scale = least-squares residual
mean square (published sigma^2 reproduced; s0 = 10^8 of the paper gives 9.2e5); lme4 in Suggests (ge_model only).

Last tests: all green; Tables 1 and B1 reproduced within 0.25 posterior SD (test-crossa2011.R).
