##################################
# Description of all data sets
################################

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
#' DE studies. Gene length is required for estimation of power of read count
#' based single cell technologies and reported for part of the studies.
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
#' fitted from example PBMC 10X data set
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

#' Parameters for gamma umi fit - Drop-seq
#'
#' Parameters to estimate parameter for the gamma mixed distribution of the
#' gene means with the mean UMI counts per cell using a linear function,
#' fitted from example lung Drop-seq data set
#'
#' @docType data
#'
#' @usage data(gammaMixedFitsDrop)
#'
#' @format Data frame with parameter of UMI-gamma parameter fit
#'
#' @keywords datasets
#'
#' @examples
#' data(gammaMixedFitsDrop)
#' head(gamma.mixed.fits.drop)
"gamma.mixed.fits.drop"

#' Parameters for gamma read fit - Smart-seq2 data
#'
#' Parameters to estimate parameter for the gamma mixed distribution of the
#' gene means with the mean read counts per cell using a linear function,
#' fitted from example pancreas Smart-seq2 data
#'
#' @docType data
#'
#' @usage data(gammaMixedFitsSmart)
#'
#' @format Data frame with parameter of read-gamma parameter fit
#'
#' @keywords datasets
#'
#' @examples
#' data(gammaMixedFitsSmart)
#' head(gamma.mixed.fits.smart)
"gamma.mixed.fits.smart"

#' Probability parameter for the gamma umi fit
#'
#' Probability parameter of the first component of the gamma mixed model,
#' in contrast to the other parameters independent of the mean UMI counts per cell,
#' fitted from example PBMC 10X data set (mean value per cell type)
#'
#' @docType data
#'
#' @usage data(gammaProbs)
#'
#' @format Data frame with probability parameter of UMI-gamma parameter fit
#'
#' @keywords datasets
#'
#' @examples
#' data(gammaProbs)
#' head(gamma.probs)
"gamma.probs"

#' Probability parameter for the gamma umi fit - Drop-seq data
#'
#' Probability parameter of the first component of the gamma mixed model,
#' in contrast to the other parameters independent of the mean UMI counts per cell,
#' fitted from example PBMC 10X data set (mean value per cell type)
#'
#' @docType data
#'
#' @usage data(gammaProbsDrop)
#'
#' @format Data frame with probability parameter of UMI-gamma parameter fit
#'
#' @keywords datasets
#'
#' @examples
#' data(gammaProbsDrop)
#' head(gamma.probs.drop)
"gamma.probs.drop"

#' Probability parameter for the gamma read fit - Smart-seq2 data
#'
#' Probability parameter of the first component of the gamma mixed model,
#' in contrast to the other parameters independent of the mean read counts per cell,
#' fitted from pancreas Smart-seq2 data (mean value per cell type)
#'
#' @docType data
#'
#' @usage data(gammaProbs)
#'
#' @format Data frame with probability parameter of read-gamma parameter fit
#'
#' @keywords datasets
#'
#' @examples
#' data(gammaProbsSmart)
#' head(gamma.probs.smart)
"gamma.probs.smart"

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

#' Parameters for mean dispersion fit - Drop-seq data
#'
#' Parameters to estimate the dispersion parameter of each gene with the
#' mean of each gene, using the mean-dispersion function of DESeq
#'
#' @docType data
#'
#' @usage data(dispFunParamDrop)
#'
#' @format Data frame with parameter of mean-dispersion function
#'
#' @keywords datasets
#'
#' @examples
#' data(dispFunParamDrop)
#' head(disp.fun.param.drop)
"disp.fun.param.drop"

#' Parameters for mean dispersion fit - Smart-seq2 data
#'
#' Parameters to estimate the dispersion parameter of each gene
#' linear dependent on the read depth for the mean-dispersion function of DESeq
#'
#' @docType data
#'
#' @usage data(dispFunParamSmart)
#'
#' @format Data frame with parameter of mean-dispersion function
#'
#' @keywords datasets
#'
#' @examples
#' data(dispFunParamSmart)
#' head(disp.fun.param.smart)
"disp.fun.param.smart"

#' Example count matrices to test the expression curve estimation
#'
#' Small subset from the 10X PBMC data to visualize expression curve
#' estimation in the introduction vignette
#'
#' @docType data
#'
#' @usage data(countMatrixExample)
#'
#' @format List with four count matrices (complete and subsampled to 75%, 50% and 25%)
#'
#' @keywords datasets
#'
#' @examples
#' data(countMatrixExample)
#' head(count.matrix.example)
"count.matrix.example"
