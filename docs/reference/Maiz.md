# Maize multi-environment trial

Grain yield of 9 maize (*Zea mays* L.) genotypes in 20 sites, a
randomized complete block design with 4 replicates in every site: the
example of Crossa et al. (2011), also used by Burgueno et al. (2008).

## Usage

``` r
Maiz
```

## Format

A `data.frame` with 720 rows and 7 columns:

- `site`:

  site (environment), 1 to 20.

- `entry`:

  genotype, 1 to 9.

- `fold10`, `fold5`, `fold3`:

  fold numbers of 10-, 5- and 3-fold cross-validation (by genotype).

- `rep`:

  replicate within the site, 1 to 4.

- `y`:

  grain yield (kg/ha).

## Source

International Maize and Wheat Improvement Center (CIMMYT).

## Details

The grand mean is 4858 kg/ha; the least-squares singular values of the
interaction are 5923 and 3070 for the first two terms (Crossa et al.
2011).

**Change from the earlier versions.** The help page described another
data set (1320 rows; columns Gen, Institute, Rep, Block, Env, Yield).

## References

Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K.,
Liu, G. and Cornelius, P. L. (2011). Bayesian estimation of the additive
main effects and multiplicative interaction model. *Crop Science*,
**51**(4), 1458–1469.
[doi:10.2135/cropsci2010.06.0343](https://doi.org/10.2135/cropsci2010.06.0343)

Burgueno, J., Crossa, J., Cornelius, P. L. and Yang, R.-C. (2008). Using
factor analytic models for joining environments and genotypes without
crossover genotype x environment interaction. *Crop Science*, **48**(4),
1291–1305.

## Author

Muhammad Yaseen (<myaseen208@gmail.com>), Jose Crossa
(<j.crossa@cgiar.org>)

## Examples

``` r
data(Maiz)
str(Maiz)
#> 'data.frame':    720 obs. of  7 variables:
#>  $ site  : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ entry : int  1 1 1 1 2 2 2 2 3 3 ...
#>  $ fold10: int  1 1 1 1 2 2 2 2 3 3 ...
#>  $ rep   : int  1 2 3 4 1 2 3 4 1 2 ...
#>  $ y     : num  3700 2948 3745 4093 3801 ...
#>  $ fold5 : int  1 1 1 1 2 2 2 2 3 3 ...
#>  $ fold3 : int  1 1 1 1 2 2 2 2 3 3 ...
ge_mean(
    .data = Maiz
  , .y    = y
  , .gen  = entry
  , .env  = site
  )$grand_mean
#> [1] 4858.138
```
