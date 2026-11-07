#!/usr/bin/env Rscript
# inst/notation/Notations.tex from inst/notation/notation.csv, so the two never drift apart (the document step of
# tools/check-all.R runs it). Symbols are split at top-level commas; the argument, when given, follows the meaning.
pkg <- read.dcf("DESCRIPTION", fields = "Package")[1L, 1L]
csv <- "inst/notation/notation.csv"
if (!file.exists(csv)) quit(save = "no")
d <- utils::read.csv(csv, colClasses = "character", encoding = "UTF-8")
split_top <- function(s) {
  ch <- strsplit(s, "")[[1L]]
  depth <- cumsum(ch == "{") - cumsum(ch == "}")
  cut <- which(ch == "," & depth == 0L)
  trimws(substring(s, c(1L, cut + 1L), c(cut - 1L, nchar(s))))
}
tt <- function(x) gsub("_", "\\\\_", x)
rows <- vapply(seq_len(nrow(d)), function(i) {
  sym <- paste0("$", split_top(d$symbol[i]), "$", collapse = ", ")
  m <- gsub("%", "\\\\%", d$meaning[i])
  if (nzchar(d$argument[i])) m <- sprintf("%s (\\texttt{%s})", m, tt(d$argument[i]))
  sprintf("%s & %s\\\\", sym, m)
}, "")
out <- c(
    sprintf("%% Notation of %s: the appendix used in articles and the reference for the help pages and vignettes.", pkg)
  , "% Rule: scalars upper case in normal face; indices lower-case italic; vectors lower case in bold; matrices upper case"
  , "% in bold; Greek letters follow the same rule (\\boldsymbol for bold Greek); sets calligraphic; operators upright."
  , "% Generated from notation.csv by tools/notation-tex.R: edit the CSV, not this file."
  , "\\begin{longtable}{>{\\raggedright\\arraybackslash}p{0.28\\textwidth}>{\\raggedright\\arraybackslash}p{0.66\\textwidth}}"
  , sprintf("\\caption{Notation of \\texttt{%s}.}\\\\", pkg), "\\toprule", "Symbol & Meaning\\\\", "\\midrule", "\\endhead"
  , rows, "\\bottomrule", "\\end{longtable}")
writeLines(out, "inst/notation/Notations.tex", useBytes = TRUE)
