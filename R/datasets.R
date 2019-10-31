#' eQTL priors
#'
#' Data frame with eQTL priors (beta values and gene expression ranks) for
#' calculation of the eQTL power. Taken from different celltype sorted bulk
#' eQTL studies.
#'
#' @docType data
#'
#' @usage data(eQTLRefStudy)
#'
#' @format Data frame with prior information of eQTL genes
#'
#' @keywords datasets
#'
#' @examples
#' data(eQTLRefStudy)
#' head(eqtl.ref.study)
"eqtl.ref.study"

#' DE priors
#'
#' Data frame with DE priors (fold changes and gene expression ranks) for
#' calculation of the DE power. Taken from different celltype sorted bulk
#' DE studies.
#'
#' @docType data
#'
#' @usage data(DERefStudy)
#'
#' @format Data frame with prior information of DE genes
#'
#' @keywords datasets
#'
#' @examples
#' data(DERefStudy)
#' head(de.ref.study)
"de.ref.study"

#' Parameters for read umi fit
#'
#' Parameters to estimate the mean UMI counts per cell with the mean read depth,
#' using a logarithmic function, fitted from our example PBMC data set
#'
#' @docType data
#'
#' @usage data(readDepthUmiFit)
#'
#' @format Data frame with parameter of read-UMI fit
#'
#' @keywords datasets
#'
#' @examples
#' data(readDepthUmiFit)
#' head(read.umi.fit)
"read.umi.fit"

#' Parameters for gamma umi fit
#'
#' Parameters to estimate parameter for the gamma mixed distribution of the
#' gene means with the mean UMI counts per cell using a linear function,
#' fitted from our example PBMC data set
#'
#' @docType data
#'
#' @usage data(gammaUmiFits)
#'
#' @format Data frame with parameter of UMI-gamma parameter fit
#'
#' @keywords datasets
#'
#' @examples
#' data(gammaUmiFits)
#' head(gamma.mixed.fits)
"gamma.mixed.fits"

#' Parameters for mean dispersion fit
#'
#' Parameters to estimate the dispersion parameter of each gene with the
#' mean of each gene, using the mean-dispersion function of DESeq
#'
#' @docType data
#'
#' @usage data(dispFunParam)
#'
#' @format Data frame with parameter of mean-dispersion function
#'
#' @keywords datasets
#'
#' @examples
#' data(dispFunParam)
#' head(disp.fun.param)
"disp.fun.param"
