# Writes inst/CITATION (the package and the method of every analysis) and inst/REFERENCES.bib (every reference) from
# the registry in R/references.R. Run after a change of the registry: Rscript tools/write-citation.R
pkgload::load_all(".", quiet = TRUE, export_all = TRUE)
refs <- .bam_refs_list()
# the methods cited by citation(): every analysis of the package
cite_keys <- c("Crossa2011", "PerezElizalde2012", "Dias2022", "FinlayWilkinson1963", "EberhartRussell1966",
               "Shukla1972", "LianCampos2016", "Oliveira2018", "EdwardsJannink2006", "Eskridge1990", "Wricke1962",
               "LinBinns1988", "Annicchiarico1992", "Purchase2000", "Sneller1997", "Zobel1994", "Zali2012",
               "Olivoto2019")
q <- function(x) paste0("\"", gsub("\"", "\\\\\"", gsub("\\\\", "\\\\\\\\", x)), "\"")
# a long string as paste() of pieces of at most 90 characters
qlong <- function(x, indent) {
  w <- strwrap(x, width = 90L)
  if (length(w) == 1L) return(q(x))
  pieces <- vapply(w, q, "")
  paste0("paste(", paste(pieces, collapse = paste0(",\n", strrep(" ", indent + 6L))), ")")
}
person_code <- function(p) {
  paste0("person(", if (length(p$given)) paste0("c(", paste(q(p$given), collapse = ", "), ")") else "NULL", ", ",
         q(p$family), ")")
}
entry_code <- function(r, header = NULL) {
  b <- unclass(r$bib)[[1L]]
  fields <- setdiff(names(b), c("author", "title"))
  au <- paste0("c(", paste(vapply(b$author, person_code, ""), collapse = paste0(",\n", strrep(" ", 18L))), ")")
  lines <- c(sprintf("    bibtype     = %s", q(attr(b, "bibtype"))), sprintf("  , key         = %s", q(r$key)),
             if (!is.null(header)) sprintf("  , header      = %s", q(header)),
             sprintf("  , title       = %s", qlong(b$title, 16L)), sprintf("  , author      = %s", au),
             vapply(fields, function(f) sprintf("  , %-11s = %s", f, q(b[[f]])), ""),
             sprintf("  , textVersion = %s", qlong(r$text, 16L)))
  paste0("bibentry(\n", paste(lines, collapse = "\n"), "\n  )")
}
old <- readLines("inst/CITATION")
cut <- grep("^bibentry\\(", old)[2L]
head <- old[seq_len(cut - 1L)]
body <- vapply(seq_along(cite_keys), function(i) entry_code(refs[[cite_keys[i]]],
              if (i == 1L) "and for the methods (bam_cite() gives the references of a single analysis):"), "")
writeLines(c(head, body), "inst/CITATION")
bib <- do.call(c, lapply(refs, `[[`, "bib"))
writeLines(c("% References of the methods of bayesammi (written by tools/write-citation.R from R/references.R)", "",
             format(bib, style = "Bibtex"), ""), "inst/REFERENCES.bib")
cat("inst/CITATION:", length(cite_keys) + 1L, "entries; inst/REFERENCES.bib:", length(refs), "entries\n")
