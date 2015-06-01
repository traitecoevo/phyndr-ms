## Ported from callr for the meantime:
run_system <- function(command, args, env=character(), max_lines=20,
                       p=0.8) {
  res <- suppressWarnings(system2(command, args,
                                  env=env, stdout=TRUE, stderr=TRUE))
  ok <- attr(res, "status")
  if (!is.null(ok) && ok != 0) {
    max_nc <- getOption("warning.length")

    cmd <- paste(c(env, shQuote(command), args), collapse = " ")
    msg <- sprintf("Running command:\n  %s\nhad status %d", cmd, ok)
    errmsg <- attr(cmd, "errmsg")
    if (!is.null(errmsg)) {
      msg <- c(msg, sprintf("%s\nerrmsg: %s", errmsg))
    }
    sep <- paste(rep("-", getOption("width")), collapse="")

    ## Truncate message:
    if (length(res) > max_lines) {
      n <- ceiling(max_lines * p)
      res <- c(head(res, ceiling(max_lines - n)),
               sprintf("[[... %d lines dropped ...]]", length(res) - max_lines),
               tail(res, ceiling(n)))
    }

    ## compute the number of characters so far, including three new lines:
    nc <- (nchar(msg) + nchar(sep) * 2) + 3
    i <- max(1, which(cumsum(rev(nchar(res) + 1L)) < (max_nc - nc)))
    res <- res[(length(res) - i + 1L):length(res)]
    msg <- c(msg, "Program output:", sep, res, sep)
    stop(paste(msg, collapse="\n"))
  }
  invisible(res)
}

## Ported from traitecoevo/plant_paper until included in remake:
latex_build <- function(filename, bibliography=NULL,
                        chdir=TRUE, interaction="nonstopmode",
                        max_attempts=5L, clean=FALSE, engine="pdflatex") {
  if (chdir && dirname(filename) != "") {
    owd <- setwd(dirname(filename))
    on.exit(setwd(owd))
    filename <- basename(filename)
  }

  res <- run_latex(filename, interaction, engine)
  if(engine=="xelatex") {
      res <- run_latex(filename, interaction, engine)    
  }
  if (!is.null(bibliography)) {
    run_bibtex(filename)
    res <- run_latex(filename, interaction, engine)
  }

  pat <- c("Rerun to get cross-references right", # labels
           "Rerun to get citations correct",      # bibtex
           "Rerun to get outlines right")         # tikz
  isin <- function(p, x) {
    any(grepl(p, x))
  }
  for (i in seq_len(max_attempts)) {
    if (any(vapply(pat, isin, logical(1), res))) {
      res <- run_latex(filename, interaction, engine)
    } else {
      break
    }
  }

  if (clean) {
    latex_clean(filename)
  }

  invisible(NULL)
}

latex_clean <- function(filename) {
  filebase <- sub(".tex$", "", filename)
  exts <- c(".log", ".aux", ".bbl", ".blg", ".fls", ".out", ".snm",
            ".nav", ".tdo", ".toc")
  aux <- paste0(filebase, exts)
  file.remove(aux[file.exists(aux)])
}

run_latex <- function(filename, interaction="nonstopmode", engine="pdflatex") {
  args <- c(paste0("-interaction=", interaction),
            "-halt-on-error",
            filename)
  run_system(Sys_which(engine), args)
}

run_bibtex <- function(filename) {
  run_system(Sys_which("bibtex"), sub(".tex$", "", filename))
}

Sys_which <- function(x) {
  ret <- Sys.which(x)
  if (ret == "") {
    stop(sprintf("%s not found in $PATH", x))
  }
  ret
}

