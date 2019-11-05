# powerScPop (working title)

R package for the manuscript "Design of single cell transcriptomics experiments for differential expression and eQTL analysis between samples"

## Installation

You will need the latest version of devtools. First install the release version:

```R
install.packages("devtools")
```

then update it to the latest developement version:

```R
devtools::install_github("hadley/devtools")
```

Finally you can install the latest development version of QTLnetwork from github with:

```R
devtools::install_github("KatharinaSchmid/powerScPop")
```

## Shiny app

Run the Shiny app of the package with the following commands:

```{R}
library(powerScPop)
runShiny()
```

## Reproduction of paper plots

All plots shown in the manuscript were created with the package. To reproduce them, the code can be found in the corresponding vignette called 'reproduce-paper-plots'. It is not created by default when downloading the package, instead the specific parameter need to be set:

```{R}
devtools::install_github("KatharinaSchmid/powerScPop",build_vignettes=TRUE)
browseVignettes("powerScPop")
```
