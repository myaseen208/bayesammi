### Name: biplots
### Title: Biplots
### Aliases: biplots biplots.default

### ** Examples

## Not run: 
##D data(Maiz)
##D 
##D fm1 <-
##D   bayes_ammi(
##D     .data = Maiz,
##D     .y = y,
##D     .gen = entry,
##D     .env = site,
##D     .rep = rep,
##D     .nIter = 200
##D   )
##D 
##D library(ggplot2)
##D 
##D output_05 <- biplots(model = fm1, plot_stable = TRUE, plot_unstable = TRUE, pb = 0.05)
##D output_05
##D 
##D output_95 <- biplots(model = fm1, plot_stable = TRUE, plot_unstable = TRUE, pb = 0.95)
##D output_95
## End(Not run)




