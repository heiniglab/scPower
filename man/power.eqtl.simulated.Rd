% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/power.R
\name{power.eqtl.simulated}
\alias{power.eqtl.simulated}
\title{Power calculation for an eQTL gene using simulations}
\usage{
power.eqtl.simulated(
  count.mean,
  heritability,
  sig.level,
  nSamples,
  rep.times = 100
)
}
\arguments{
\item{count.mean}{Expression mean in the pseudobulk}

\item{heritability}{Heritability of the trait}

\item{sig.level}{Significane threshold}

\item{nSamples}{Sample size}

\item{rep.times}{Number of repetitions used for the}
}
\value{
Power to detect the eQTL gene
}
\description{
The eQTL power is independent of the mean except for very small mean values,
where small effect sizes might not be detectable due to the discrete counts
In these cases, the power can instead be simulated.
}
