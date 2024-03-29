% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expression_fit.R
\name{estimate.exp.prob.count.param}
\alias{estimate.exp.prob.count.param}
\title{Estimate gene expression probability based on experimental parameters - variant for Smart-seq or when calculating
directly based on UMI counts (without read UMI fit)}
\usage{
estimate.exp.prob.count.param(
  nSamples,
  nCellsCt,
  meanCellCounts,
  gamma.mixed.fits,
  ct,
  disp.fun.param,
  min.counts = 3,
  perc.indiv.expr = 0.5,
  cutoffVersion = "absolute",
  nGenes = 21000,
  samplingMethod = "quantiles",
  countMethod = "UMI",
  gene_length_prior = NULL
)
}
\arguments{
\item{nSamples}{Sample size}

\item{nCellsCt}{Mean number of cells per individual and cell type}

\item{gamma.mixed.fits}{Data frame with gamma mixed fit parameters for each cell type
(required columns: parameter, ct (cell type), intercept, meanUMI (slope))}

\item{ct}{Cell type of interest (name from the gamma mixed models)}

\item{disp.fun.param}{Data frame with parameters to fit the dispersion parameter dependent on the mean
(required columns: ct (cell type), asymptDisp, extraPois (both from taken from DEseq))}

\item{min.counts}{Expression cutoff in one individual: if cutoffVersion=absolute,
more than this number of UMI counts for each gene per individual and
cell type is required; if cutoffVersion=percentage, more than this percentage
of cells need to have a count value large than 0}

\item{perc.indiv.expr}{Expression cutoff on the population level: if number < 1, percentage of
individuals that need to have this gene expressed to define it as globally expressed;
if number >=1 absolute number of individuals that need to have this gene expressed}

\item{cutoffVersion}{Either "absolute" or "percentage" leading to different
interpretations of min.counts (see description above)}

\item{nGenes}{Number of genes to simulate (should match the number of genes used for the fitting)}

\item{samplingMethod}{Approach to sample the gene mean values (either taking quantiles or random sampling)}

\item{countMethod}{Specify if it is a UMI or read based method (by "UMI" or "read"). For a read based method,
the dispersion function is fitted linear, for a UMI based method constant.}

\item{gene_length_prior}{Vector with gene length, ordered by expression rank 
(the first entry for the highest expressed gene); only required for Smart-seq (countMethod="read")}

\item{meanCellCount}{Either mean UMI counts per cell for UMI-based methods or mean read counts per cell for read-based methods}
}
\value{
Vector with expression probabilities for each gene
}
\description{
This is an adaption of function estimate.exp.prob.param were directly the mean UMI parameter can be set in the
variable meanCellCount (no read.umi.fit and read depths required). This version can also be used to model read based single cell methods,
such as Smart-seq, by setting the parameter meanCellCounts to their read depth
}
