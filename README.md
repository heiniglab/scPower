# scPower

R package for the manuscript "Design and power analysis for multi-sample single cell genomics experiments"

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
devtools::install_github("KatharinaSchmid/scPower")
```

## Shiny app

Run the Shiny app of the package with the following commands:

```{R}
library(scPower)
runShiny()
```

The shiny app is also available at the webpage: scpower.helmholtz-muenchen.de

## Vignettes for introduction and for reproduction of paper plots

The package provides two vignettes for the user. The first one, called "introduction-scPower", is a user manual describing all functions with small toy examples. The second one, called "reproduce-paper-plots", shows how all plots in the manuscript were created with the package. The vignettes are not incorporated by default when downloading the package, instead the specific parameter need to be set:

```{R}
devtools::install_github("KatharinaSchmid/scPower",build_vignettes=TRUE)
browseVignettes("scPower")
```
