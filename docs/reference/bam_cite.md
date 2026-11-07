# References of a result or of the package

`bam_cite()` gives the references of the methods behind a result of the
package, or the references of every method in the package, as text or
BibTeX. `citation("bayesammi")` cites the package and its main methods.

## Usage

``` r
bam_cite(x = NULL, style = c("text", "bibtex"))
```

## Arguments

- x:

  a result of a function of the package, or `NULL` (default) for every
  reference of the package.

- style:

  `"text"` (default) or `"bibtex"`.

## Value

An object of class `bibentry` (printed in the chosen style).

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
bam_cite(
    x     = NULL
  , style = "text"
  )
#> Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K., Liu, G. and Cornelius, P. L. (2011). Bayesian estimation of the additive main effects and multiplicative interaction model. Crop Science, 51(4), 1458-1469. doi:10.2135/cropsci2010.06.0343
#> 
#> Perez-Elizalde, S., Jarquin, D. and Crossa, J. (2012). A general Bayesian estimation method of linear-bilinear models applied to plant breeding trials with genotype x environment interaction. Journal of Agricultural, Biological, and Environmental Statistics, 17(1), 15-37. doi:10.1007/s13253-011-0063-9
#> 
#> Viele, K. and Srinivasan, C. (2000). Parsimonious estimation of multiplicative interaction in analysis of variance using Kullback-Leibler information. Journal of Statistical Planning and Inference, 84(1-2), 201-219.
#> 
#> Hoff, P. D. (2009). Simulation of the matrix Bingham-von Mises-Fisher distribution, with applications to multivariate and relational data. Journal of Computational and Graphical Statistics, 18(2), 438-456. doi:10.1198/jcgs.2009.07177
#> 
#> Yan, W., Hunt, L. A., Sheng, Q. and Szlavnics, Z. (2000). Cultivar evaluation and mega-environment investigation based on the GGE biplot. Crop Science, 40(3), 597-605. doi:10.2135/cropsci2000.403597x
#> 
#> Dias, K. O. G., dos Santos, J. P. R., Krause, M. D., Piepho, H.-P., Guimaraes, L. J. M., Pastina, M. M. and Garcia, A. A. F. (2022). Leveraging probability concepts for cultivar recommendation in multi-environment trials. Theoretical and Applied Genetics. doi:10.1007/s00122-022-04041-y
#> 
#> Finlay, K. W. and Wilkinson, G. N. (1963). The analysis of adaptation in a plant-breeding programme. Australian Journal of Agricultural Research, 14(6), 742-754. doi:10.1071/AR9630742
#> 
#> Eberhart, S. A. and Russell, W. A. (1966). Stability parameters for comparing varieties. Crop Science, 6(1), 36-40. doi:10.2135/cropsci1966.0011183X000600010011x
#> 
#> Shukla, G. K. (1972). Some statistical aspects of partitioning genotype-environmental components of variability. Heredity, 29(2), 237-245. doi:10.1038/hdy.1972.87
#> 
#> Wricke, G. (1962). Uber eine Methode zur Erfassung der okologischen Streubreite in Feldversuchen. Zeitschrift fur Pflanzenzuchtung, 47, 92-96.
#> 
#> Lin, C. S. and Binns, M. R. (1988). A superiority measure of cultivar performance for cultivar x location data. Canadian Journal of Plant Science, 68(1), 193-198. doi:10.4141/cjps88-018
#> 
#> Annicchiarico, P. (1992). Cultivar adaptation and recommendation from alfalfa trials in Northern Italy. Journal of Genetics and Breeding, 46, 269-278.
#> 
#> Eskridge, K. M. (1990). Selection of stable cultivars using a safety-first rule. Crop Science, 30, 369-374.
#> 
#> Edwards, J. W. and Jannink, J.-L. (2006). Bayesian modeling of heterogeneous error and genotype x environment interaction variances. Crop Science, 46(2), 820. doi:10.2135/cropsci2005.0164
#> 
#> Lian, L. and de los Campos, G. (2016). FW: an R package for Finlay-Wilkinson regression that incorporates genomic/pedigree information and covariance structures between environments. G3: Genes, Genomes, Genetics, 6(3), 589-597. doi:10.1534/g3.115.026328
#> 
#> Oliveira, T. R. A., Carvalho, H. W. L., Nascimento, M., Costa, E. F. N., Amaral Junior, A. T., Gravina, G. A. and Carvalho Filho, J. L. S. (2018). The Eberhart and Russell's Bayesian method used as an instrument to select maize hybrids. Euphytica, 214, 64. doi:10.1007/s10681-018-2146-y
#> 
#> Purchase, J. L., Hatting, H. and van Deventer, C. S. (2000). Genotype x environment interaction of winter wheat (Triticum aestivum L.) in South Africa: II. Stability analysis of yield performance. South African Journal of Plant and Soil, 17(3), 101-107.
#> 
#> Sneller, C. H., Kilgore-Norquest, L. and Dombek, D. (1997). Repeatability of yield stability statistics in soybean. Crop Science, 37(2), 383-390.
#> 
#> Zobel, R. W. (1994). Stress resistance and root systems. In Proceedings of the Workshop on Adaptation of Plants to Soil Stress, INTSORMIL Publication 94-2, 80-99. Institute of Agriculture and Natural Resources, University of Nebraska-Lincoln.
#> 
#> Zali, H., Farshadfar, E., Sabaghpour, S. H. and Karimizadeh, R. (2012). Evaluation of genotype x environment interaction in chickpea using measures of stability from AMMI model. Annals of Biological Research, 3(7), 3126-3136.
#> 
#> Olivoto, T., Lucio, A. D. C., da Silva, J. A. G., Marchioro, V. S., de Souza, V. Q. and Jost, E. (2019). Mean performance and stability in multi-environment trials I: combining features of AMMI and BLUP techniques. Agronomy Journal, 111(6), 2949-2960. doi:10.2134/agronj2019.03.0220
#> 
#> Farshadfar, E. (2008). Incorporation of AMMI stability value and grain yield in a single non-parametric index (GSI) in bread wheat. Pakistan Journal of Biological Sciences, 11(14), 1791-1796.
#> 
#> Spiegelhalter, D. J., Best, N. G., Carlin, B. P. and van der Linde, A. (2002). Bayesian measures of model complexity and fit. Journal of the Royal Statistical Society B, 64(4), 583-639. doi:10.1111/1467-9868.00353
#> 
#> Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. Journal of Machine Learning Research, 11, 3571-3594.
#> 
#> Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. Statistics and Computing, 27(5), 1413-1432. doi:10.1007/s11222-016-9696-4
#> 
#> Gelman, A., Meng, X.-L. and Stern, H. (1996). Posterior predictive assessment of model fitness via realized discrepancies. Statistica Sinica, 6(4), 733-807.
#> 
#> Gabry, J., Simpson, D., Vehtari, A., Betancourt, M. and Gelman, A. (2019). Visualization in Bayesian workflow. Journal of the Royal Statistical Society A, 182(2), 389-402. doi:10.1111/rssa.12378 
data(Maiz)
g <-
  ge_ammi(
      .data = Maiz
    , .y    = y
    , .gen  = entry
    , .env  = site
    , .rep  = rep
    )
bam_cite(
    x     = g
  , style = "bibtex"
  )
#> @Article{Crossa2011,
#>   title = {Bayesian Estimation of the Additive Main Effects and Multiplicative Interaction Model},
#>   author = {Jose Crossa and Sergio Perez-Elizalde and Diego Jarquin and Jose Miguel Cotes and Kert Viele and Genzhou Liu and Paul L. Cornelius},
#>   journal = {Crop Science},
#>   year = {2011},
#>   volume = {51},
#>   number = {4},
#>   pages = {1458--1469},
#>   doi = {10.2135/cropsci2010.06.0343},
#> }
```
