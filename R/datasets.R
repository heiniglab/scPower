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

#' Precalcuated size estimates for eQTL power simulation
#'
#' Precalculated size estimates for eQTL power simulation
#' (using the function estimate.size.simulation) to speed simulation
#'
#' @docType data
#'
#' @usage data(sizeEstimates)
#'
#' @format Data frame with precalculated size estimates (depending on Rsq, af and mean)
#'
#' @keywords datasets
#'
#' @examples
#' data(sizeEstimates)
#' head(size.estimates)
"size.estimates"

#' Precalcuated p-values from eQTL simulation
#'
#' Precalculated p-values from eQTL simulation to speed calculation time
#'
#' @docType data
#'
#' @usage data(simEqtlPvals)
#'
#' @format Matrix of with 100 p-values (columns) and rownames in the format
#' "mean_Rsq_sampleSize" showing the corresponding parameter combination with
#' which the p-values were generated
#'
#' @keywords datasets
#'
#' @examples
#' data(simEqtlPvals)
#' head(sim.eqtl.pvals)
"sim.eqtl.pvals"

#' Example count matrices to test the expression curve estimation
#'
#' Small subset from the 10X PBMC data to visualize expression curve
#' estimation in the introduction vignette
#'
#' @docType data
#'
#' @usage data(countMatrixExample)
#'
#' @format List with four count matrices (complete and subsampled to 75\%, 50\% and 25\%)
#'
#' @keywords datasets
#'
#' @examples
#' data(countMatrixExample)
#' head(count.matrix.example)
"count.matrix.example"

#' Example annotation data frame matching the example count matrix to 
#' test the expression curve estimation
#'
#' Contains three columns necessary for fitting the expression priors 
#' (cell barcode, cell type and donor)
#'
#' @docType data
#'
#' @usage data(annotDfExample)
#'
#' @format Data frame with three columns (cell, cell_type, donor)
#'
#' @keywords datasets
#'
#' @examples
#' data(annotDfExample)
#' head(annot.df)
"annot.df"

#' Example Smart-seq2 count matrix for the tutorial
#'
#' Contains a read count matrix of genes times cells
#'
#' @docType data
#'
#' @usage data(smartseqExample)
#'
#' @format Read count matrix of genes times cells
#'
#' @keywords datasets
#'
#' @examples
#' data(smartseqExample)
#' dim(counts_smartseq)
"counts_smartseq"

#' Example Smart-seq2 annotation data frame (matching the Smart-seq2 count matrix)
#'
#' Contains an annotation data frame with cell type and donor information per cell
#'
#' @docType data
#'
#' @usage data(smartseqExample)
#'
#' @format Annotation data frame with cell type and donor information per cell
#'
#' @keywords datasets
#'
#' @examples
#' data(smartseqExample)
#' head(annot_smartseq)
"annot_smartseq"

#' Example gene length data frame
#'
#' Contains gene length (only exonic regions) in bp for the genes from the Smart-seq2 example data frame;
#' information obtained from UCSC (gene annotation GRCh37hg19)
#'
#' @docType data
#'
#' @usage data(smartseqExample)
#'
#' @format Gene length data frame with gene symbol and gene length in bp (exonic regions)
#'
#' @keywords datasets
#'
#' @examples
#' data(smartseqExample)
#' head(gene_length)
"gene_length"

#' Observed gene counts from multiple data sets to reproduce figures curves in the vignette
#' reproduce-figure-plots faster
#'
#' @docType data
#'
#' @usage data(precalculatedObservedGeneCounts)
#'
#' @format Data frame with observed gene counts from multiple data sets
#'
#' @keywords datasets
#'
#' @examples
#' data(precalculatedObservedGeneCounts)
#' head(observed.gene.counts)
"observed.gene.counts"

#' Optimal parameters for increasing budget and simulated and observed priors from different
#' single cell technologies to reproduce figures curves in the vignette reproduce-figure-plots faster
#'
#' @docType data
#'
#' @usage data(precalculatedBudgetOptim)
#'
#' @format Data frame with optimal parameter for increasing budget and different priors
#'
#' @keywords datasets
#'
#' @examples
#' data(precalculatedBudgetOptim)
#' head(budget.optimization)
"budget.optimization"

#' Results of DE power generated by simulations with powsimR
#'
#' @docType data
#'
#' @usage data(simulatedDEValues)
#'
#' @format Data frame with simulated DE power
#'
#' @keywords datasets
#'
#' @examples
#' data(simulatedDEValues)
#' head(simulated.DE.values)
"simulated.DE.values"

#' Results of eQTL power generated by self-implemented simulations
#' (see vignette reproduce-paper-plots)
#'
#' @docType data
#'
#' @usage data(simulatedEQTLValues)
#'
#' @format Data frame with simulated eQTL power
#'
#' @keywords datasets
#'
#' @examples
#' data(simulatedEQTLValues)
#' head(simulated.eQTL.values)
"simulated.eQTL.values"
