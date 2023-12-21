# scPower

scPower is a R package for design and power analysis of cell type specific interindividual DE and eQTL studies using single cell RNA-seq. It enables the user to calculate the power for a given experimental setup and to choose for a restricted budget the optimal combination of experimental parameters which maximizes the power. Necessary experimental priors, e.g. effect sizes and expression distributions, can be taken from example data sets, saved in the package, or estimated from new data sets. The tool was evaluated with data from different tissues and single cell technologies, based on UMI counts and read counts. 

The calculation can also performed using a graphical interface of a shiny app or at our webpage [scpower](http://scpower.helmholtz-muenchen.de), which is however temporarily not available ([see section below](#update-web-server-for-scpower-temporarily-down)).

A short tutorial is given in the vignette [introduction-scPower](vignettes/introduction-scPower.pdf).

A detailed description of all methods and citation of all used tools and data sets can be found in the associated paper 

Schmid, K. T., Höllbacher, B., Cruceanu, C., Boettcher, A., Lickert, H., Binder, E. B., Theis, F. J., & Heinig, M. (2021). scPower accelerates and optimizes the design of multi-sample single cell transcriptomic studies. Nature Communications. https://doi.org/10.1038/s41467-021-26779-7

An explanation how the plots in the paper were generated is shown in the second vignette [reproduce-paper-plots](vignettes/reproduce-paper-plots.pdf) and an example how to combine our model with more complex designs is shown in the third vignette [extension-complex-design](vignettes/extension-complex-design.pdf).

## A new and improved web server for scPower is now available

A new version of the scPower webserver https://scpower.helmholtz-muenchen.de/ is now featuring parameters derived atlas scale data sets, enabling power analysis across 66 tissues, 891 cell types, profiled with 6 different technology platforms.

## Installation

You will need the latest version of devtools. First install the release version:

```R
install.packages("devtools")
```

then update it to the latest developement version:

```R
devtools::install_github("hadley/devtools")
```

Finally you can install the latest development version of scPower from github with:

```R
devtools::install_github("heiniglab/scPower")
```

If you have problems installing the package, please try to install necessary packages yourself from CRAN.

```R
#CRAN packages
install.packages(c("pwr","MKmisc","reshape2","HardyWeinberg","plotly", "shiny"))
```

## Shiny app

Run the Shiny app of the package with the following commands:

```{R}
library(scPower)
runShiny()
```

The shiny app is also available at the webpage [scpower](http://scpower.helmholtz-muenchen.de) .

## Vignettes for introduction and for reproduction of paper plots

The package provides two vignettes for the user. The first one, called [introduction-scPower](vignettes/introduction-scPower.pdf), is a user manual describing all functions with small toy examples. The second one, called [reproduce-paper-plots](vignettes/reproduce-paper-plots.pdf), shows how all plots in the manuscript were created with the package. The vignettes can be built when downloading the package with setting the option "build_vignettes=TRUE", however, this will take several minutes, as the vignettes are quite detailed. Alternatively, the pdf files are directly uploaded and can be accessed in github.

```{R}
devtools::install_github("heiniglab/scPower",build_vignettes=TRUE)
browseVignettes("scPower")
```
Be aware that additional libraries are required for building the vignettes:

```R
#CRAN packages
install.packages(c("knitr","rmarkdown","ggpubr","data.table", "viridis", "RColorBrewer","gridExtra","dplyr"))
```
