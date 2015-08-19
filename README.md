## A simple approach for maximizing the overlap of phylogenetic and comparative data

This repo contains tex files for our manuscript and R code to reproduce all example analyses. The paper describes two packages [phyndr](https://github.com/richfitz/phyndr) and [taxonlookup](https://github.com/wcornwell/taxonlookup) (the source code for the packages are in their own repos; links above).

The analyses contained in the paper can be reproduced using the [remake](https://github.com/richfitz/remake) R package.

A vignette, which works through the basic functionality of both phyndr and taxonlookup, is available [here](vignette/vignette.md).

```
## Install devtools if you do not have it
## install.packages("devtools")

## Install remake
devtools::install_github("richfitz/remake")

## Re-run the analysis with remake
remake::install_missing_packages()
remake::make()
```

This project is a collaboration between [Matthew Pennell](https://mwpennell.github.io), [Rich FitzJohn](http://richfitz.github.io/), and [Will Cornwell](http://willcornwell.org/).
