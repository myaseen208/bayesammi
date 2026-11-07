# Bayesian AMMI and GGE analysis in a web browser

`bam_app()` starts the application of bayesammi: load a trial (the
example [`Maiz`](https://myaseen208.com/bayesammi/reference/Maiz.md) or
a CSV file with one row per plot), fit the Bayesian AMMI or GGE model,
check convergence, draw the biplot with HPD regions, and turn the
posterior into decisions (which genotype wins where, probabilities of
superior performance and stability, ranks, Bayesian stability,
predictions of missing cells), and compare models. Every page shows an
interpretation computed from the results
([`bam_interpret()`](https://myaseen208.com/bayesammi/reference/bam_interpret.md)),
the output as R prints it, formatted tables, plots and the R code that
repeats the analysis, each with a Copy button.

## Usage

``` r
bam_app(launch.browser = interactive(), ...)

bam_app_object()
```

## Arguments

- launch.browser:

  passed to
  [`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html): `TRUE`
  (default in an interactive session) opens the application in the web
  browser.

- ...:

  further arguments to
  [`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html), such
  as `port`.

## Value

`bam_app()` is called for its side effect and returns the value of
[`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html) when the
application is closed; `bam_app_object()` returns the application object
(for [`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html),
[`shiny::testServer()`](https://rdrr.io/pkg/shiny/man/testServer.html)
or deployment on a server).

## Details

The application needs the suggested packages shiny, bslib and DT. Plots
are sent to the browser as WebP images when magick is installed (smaller
than PNG; PNG otherwise) and can be downloaded as PNG, PDF or WebP. Text
is set in IBM Plex Sans and output in IBM Plex Mono when the browser can
reach Google Fonts, in the system fonts otherwise.

## See also

[`bayes_ammi`](https://myaseen208.com/bayesammi/reference/bayes_ammi.md),
[`bam_interpret`](https://myaseen208.com/bayesammi/reference/bam_interpret.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
if (interactive()) {
  bam_app(
      launch.browser = TRUE
    )
}
```
